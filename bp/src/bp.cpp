#include "bp/bp.hpp"
#include "bp_internal.hpp"
#include "benchmark_internal.hpp"
#include "hp_transcript_internal.hpp"

#include <mcl/fp.hpp>
#include <algorithm>
#include <chrono>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string_view>

namespace bp {
namespace {

using Clock = std::chrono::steady_clock;
using BenchmarkProfile = benchmark_internal::Profile;
thread_local BenchmarkProfile* active_profile = nullptr;

double elapsed_ms(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

constexpr std::string_view kProtocol = "BP-IPA/BN254-G2/v1";
constexpr std::string_view kAlphaDomain = "BP-IPA-ALPHA-v1";
constexpr std::string_view kRoundDomain = "BP-IPA-ROUND-v1";
constexpr std::array<std::uint8_t, 8> kProofMagic{
    'B','P','I','P','A','G','2','1'};
constexpr std::uint16_t kProofVersion = 1;


constexpr Digest kRejectionLimit = {
    0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,
    0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,
    0xfd,0xbd,0x00,0x00,0x00,0x00,0x00,0x63,
    0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};

Bytes u64be(std::uint64_t value) {
  Bytes out;
  out.reserve(8);
  for (int shift = 56; shift >= 0; shift -= 8)
    out.push_back(static_cast<std::uint8_t>(value >> shift));
  return out;
}

void append_raw(Bytes& out, std::span<const std::uint8_t> field) {
  out.insert(out.end(), field.begin(), field.end());
}

void frame(Bytes& out, std::span<const std::uint8_t> field) {
  append_raw(out, u64be(field.size()));
  append_raw(out, field);
}

void frame(Bytes& out, std::string_view field) {
  frame(out, {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

template<class T>
Bytes encode(const T& value) {
  Bytes out(512);
  const auto written = value.serialize(out.data(), out.size());
  if (written == 0) throw std::runtime_error("MCL serialization failed");
  out.resize(written);
  return out;
}

bool valid_group(const Group& p, bool nonidentity) {
  return p.isValid() && p.isValidOrder() && (!nonidentity || !p.isZero());
}

Group add(const Group& a, const Group& b) {
  Group out;
  Group::add(out, a, b);
  return out;
}

Group mul(const Group& point, const Scalar& scalar) {
  Group out;
  Group::mul(out, point, scalar);
  return out;
}

Scalar add(const Scalar& a, const Scalar& b) {
  Scalar out;
  Scalar::add(out, a, b);
  return out;
}

Scalar mul(const Scalar& a, const Scalar& b) {
  Scalar out;
  Scalar::mul(out, a, b);
  return out;
}

Scalar inverse(const Scalar& value) {
  if (value.isZero()) throw std::invalid_argument("zero challenge");
  Scalar out;
  Scalar::inv(out, value);
  return out;
}

Scalar inner_product(std::span<const Scalar> a, std::span<const Scalar> b) {
  if (a.size() != b.size()) throw std::invalid_argument("inner-product length mismatch");
  Scalar out;
  out.clear();
  for (std::size_t i = 0; i < a.size(); ++i)
    out = add(out, mul(a[i], b[i]));
  return out;
}

Group msm(std::span<const Group> points, std::span<const Scalar> scalars) {
  if (points.size() != scalars.size()) throw std::invalid_argument("MSM length mismatch");
  const auto msm_start = Clock::now();
  Group out;
  out.clear();
  if (!points.empty()) {
    const auto copy_start = Clock::now();
    std::vector<Group> work(points.begin(), points.end());
    if (active_profile)
      active_profile->msm_input_copy_ms += elapsed_ms(copy_start);
    Group::mulVec(out, work.data(), scalars.data(), work.size());
  }
  if (active_profile)
    active_profile->proof_message_msm_ms += elapsed_ms(msm_start);
  return out;
}

bool power_of_two(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

std::size_t log2_exact(std::size_t n) {
  if (!power_of_two(n)) throw std::invalid_argument("dimension must be a power of two");
  std::size_t d = 0;
  while (n > 1) { n >>= 1; ++d; }
  return d;
}

Digest public_params_digest(const PublicParams& pp) {
  Bytes input;
  frame(input, "BP-IPA-PUBLIC-PARAMS-v1");
  frame(input, pp.group_identifier);
  frame(input, pp.scalar_modulus);
  frame(input, pp.hash_suite_identifier);
  frame(input, pp.transcript_domain);
  frame(input, u64be(pp.n));
  frame(input, u64be(pp.d));
  for (const auto& p : pp.G) frame(input, serialize(p));
  for (const auto& p : pp.H) frame(input, serialize(p));
  frame(input, serialize(pp.K));
  return sha256(input);
}

class Transcript {
 public:
  Transcript(const PublicParams& pp, const Group& Z,
             const Digest* cached_public_parameter_digest = nullptr) {
    if (pp.transcript_domain == kHpTranscriptDomain) {
      hp_.emplace(pp, Z);
      return;
    }
    Bytes input;
    frame(input, kProtocol);
    frame(input, pp.transcript_domain);
    frame(input, pp.group_identifier);
    frame(input, pp.scalar_modulus);
    frame(input, pp.hash_suite_identifier);
    frame(input, cached_public_parameter_digest
                     ? *cached_public_parameter_digest
                     : public_params_digest(pp));
    frame(input, serialize(Z));
    state_ = sha256(input);
  }

  Scalar round(std::size_t k, const Group& A, const Group& B) {
    if (hp_) return hp_->round(k, A, B);
    Bytes message;
    frame(message, kRoundDomain);
    frame(message, state_);
    frame(message, u64be(k));
    frame(message, serialize(A));
    frame(message, serialize(B));
    const Digest round_state = sha256(message);
    Scalar alpha = hash_to_nonzero(round_state);
    Bytes next;
    frame(next, kProtocol);
    frame(next, round_state);
    frame(next, serialize(alpha));
    state_ = sha256(next);
    return alpha;
  }

 private:
  static Scalar hash_to_nonzero(const Digest& input) {
    for (std::uint64_t counter = 0;; ++counter) {
      Bytes candidate;
      frame(candidate, kAlphaDomain);
      frame(candidate, input);
      frame(candidate, u64be(counter));
      const Digest hash = sha256(candidate);
      if (!std::lexicographical_compare(
              hash.begin(), hash.end(), kRejectionLimit.begin(), kRejectionLimit.end()))
        continue;
      Scalar out;
      out.setBigEndianMod(hash.data(), hash.size());
      if (!out.isZero()) return out;
    }
  }

  Digest state_{};
  std::optional<hp_internal::HpBpTranscript> hp_;
};

Group setup_base(std::string_view domain, std::span<const std::uint8_t> seed,
                 std::uint64_t index, bool has_index) {
  for (std::uint64_t counter = 0;; ++counter) {
    Bytes input;
    frame(input, kProtocol);
    frame(input, domain);
    frame(input, seed);
    if (has_index) frame(input, u64be(index));
    frame(input, u64be(counter));
    Group out;
    mcl::bn::hashAndMapToG2(out, input.data(), input.size());
    if (valid_group(out, true)) return out;
  }
}

std::vector<Scalar> fold_scalars(std::span<const Scalar> left,
                                 std::span<const Scalar> right,
                                 const Scalar& factor) {
  std::vector<Scalar> out;
  out.reserve(left.size());
  for (std::size_t i = 0; i < left.size(); ++i)
    out.push_back(add(left[i], mul(factor, right[i])));
  return out;
}

std::vector<Group> fold_groups(std::span<const Group> left,
                               std::span<const Group> right,
                               const Scalar& factor) {
  std::vector<Group> out;
  out.reserve(left.size());
  for (std::size_t i = 0; i < left.size(); ++i)
    out.push_back(add(left[i], mul(right[i], factor)));
  return out;
}

Group relation(std::span<const Group> G, std::span<const Group> H,
               const Group& K, std::span<const Scalar> x,
               std::span<const Scalar> y) {
  return add(add(msm(G, x), msm(H, y)), mul(K, inner_product(x, y)));
}

bool canonical_scalar(std::span<const std::uint8_t> bytes, Scalar& out) {
  Scalar candidate;
  if (bytes.size() != scalar_bytes() ||
      candidate.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      serialize(candidate) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = candidate;
  return true;
}

bool canonical_group(std::span<const std::uint8_t> bytes, Group& out) {
  Group candidate;
  if (bytes.size() != group_bytes() ||
      candidate.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      !valid_group(candidate, false) ||
      serialize(candidate) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = candidate;
  return true;
}

std::uint64_t read_u64(std::span<const std::uint8_t> bytes, std::size_t pos) {
  std::uint64_t out = 0;
  for (std::size_t i = 0; i < 8; ++i) out = (out << 8) | bytes[pos + i];
  return out;
}

bool valid_typed_proof_inputs(const PublicParams& pp, const Group& Z,
                              const Proof& proof) {
  if (!valid_group(Z, false) || proof.rounds.size() != pp.d) return false;
  for (const auto& round : proof.rounds)
    if (!valid_group(round.A, false) || !valid_group(round.B, false)) return false;
  return true;
}

bool verify_core(const PublicParams& pp, const Group& Z, const Proof& proof,
                 Trace* trace) {
  if (trace) { trace->challenges.clear(); trace->round_invariants.clear(); }

  const auto copy_start = Clock::now();
  std::vector<Group> G(pp.G), H(pp.H);
  Group current = Z;
  if (active_profile)
    active_profile->initial_allocation_copy_ms += elapsed_ms(copy_start);

  const auto transcript_start = Clock::now();
  Transcript transcript(pp, Z);
  if (active_profile)
    active_profile->transcript_challenge_ms += elapsed_ms(transcript_start);

  for (std::size_t round = 0; round < pp.d; ++round) {
    const std::size_t k = pp.d - round;
    const auto& message = proof.rounds[round];

    const auto challenge_start = Clock::now();
    const Scalar alpha = transcript.round(k, message.A, message.B);
    const Scalar alpha_inv = inverse(alpha);
    if (active_profile)
      active_profile->transcript_challenge_ms += elapsed_ms(challenge_start);

    const std::size_t half = G.size() / 2;
    const auto generator_start = Clock::now();
    auto Gn = fold_groups(std::span(G).first(half), std::span(G).subspan(half), alpha);
    auto Hn = fold_groups(std::span(H).first(half), std::span(H).subspan(half), alpha_inv);
    if (active_profile)
      active_profile->generator_folding_ms += elapsed_ms(generator_start);

    const auto update_start = Clock::now();
    current = add(add(current, mul(message.A, alpha_inv)), mul(message.B, alpha));
    if (active_profile)
      active_profile->verifier_z_update_ms += elapsed_ms(update_start);
    if (trace) trace->challenges.push_back(alpha);
    G = std::move(Gn);
    H = std::move(Hn);
  }

  const auto terminal_start = Clock::now();
  const Group expected = add(add(mul(G[0], proof.x_final), mul(H[0], proof.y_final)),
                             mul(pp.K, mul(proof.x_final, proof.y_final)));
  const bool accepted = current == expected;
  if (active_profile)
    active_profile->terminal_verification_ms += elapsed_ms(terminal_start);
  return accepted;
}

}

void initialize() {
  static std::once_flag once;
  std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}

Bytes serialize(const Scalar& value) { return encode(value); }
Bytes serialize(const Group& value) { return encode(value); }

Digest sha256(std::span<const std::uint8_t> input) {
  if (input.size() > std::numeric_limits<std::uint32_t>::max())
    throw std::length_error("SHA-256 input too large");
  Digest out{};
  const auto written = mcl::fp::sha256(
      out.data(), static_cast<std::uint32_t>(out.size()), input.data(),
      static_cast<std::uint32_t>(input.size()));
  if (written != out.size()) throw std::runtime_error("SHA-256 failed");
  return out;
}

std::size_t scalar_bytes() {
  initialize();
  Scalar one;
  one = 1;
  return serialize(one).size();
}

std::size_t group_bytes() {
  initialize();
  static const std::size_t size = [] {
    Group p;
    mcl::bn::hashAndMapToG2(p, "bp-size", 7);
    return serialize(p).size();
  }();
  return size;
}

PublicParams Setup(std::size_t n, std::span<const std::uint8_t> public_seed) {
  initialize();
  const std::size_t d = log2_exact(n);
  PublicParams pp;
  pp.n = n;
  pp.d = d;
  pp.G.reserve(n);
  pp.H.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    pp.G.push_back(setup_base("BP-IPA-G-v1", public_seed, i, true));
    pp.H.push_back(setup_base("BP-IPA-H-v1", public_seed, i, true));
  }
  pp.K = setup_base("BP-IPA-K-v1", public_seed, 0, false);
  return pp;
}

bool validate_public_params(const PublicParams& pp) noexcept {
  try {
    initialize();
    const bool standard_transcript =
        pp.transcript_domain == kTranscriptDomain && !pp.transcript_crs_digest;
    const bool hp_transcript =
        pp.transcript_domain == kHpTranscriptDomain && pp.transcript_crs_digest.has_value();
    if (!power_of_two(pp.n) || pp.d != log2_exact(pp.n) ||
        pp.G.size() != pp.n || pp.H.size() != pp.n ||
        pp.group_identifier != kGroupIdentifier ||
        pp.scalar_modulus != kScalarModulus ||
        pp.hash_suite_identifier != kHashSuiteIdentifier ||
        (!standard_transcript && !hp_transcript))
      return false;
    for (const auto& p : pp.G) if (!valid_group(p, true)) return false;
    for (const auto& p : pp.H) if (!valid_group(p, true)) return false;
    return valid_group(pp.K, true);
  } catch (...) { return false; }
}

Group commit(const PublicParams& pp, std::span<const Scalar> x,
             std::span<const Scalar> y) {
  if (!validate_public_params(pp)) throw std::invalid_argument("invalid public parameters");
  if (x.size() != pp.n || y.size() != pp.n)
    throw std::invalid_argument("invalid witness dimensions");
  return relation(pp.G, pp.H, pp.K, x, y);
}

Proof Prove(const PublicParams& pp, const Group& Z,
            std::span<const Scalar> x, std::span<const Scalar> y,
            Trace* trace) {
  const auto validation_start = Clock::now();
  const bool valid_pp = validate_public_params(pp);
  const bool valid_statement = valid_group(Z, false);
  const bool valid_dimensions = x.size() == pp.n && y.size() == pp.n;
  if (active_profile)
    active_profile->public_parameter_validation_ms += elapsed_ms(validation_start);
  if (!valid_pp) throw std::invalid_argument("invalid public parameters");
  if (!valid_statement) throw std::invalid_argument("invalid statement");
  if (!valid_dimensions)
    throw std::invalid_argument("invalid witness dimensions");
  if (trace) { trace->challenges.clear(); trace->round_invariants.clear(); }

  const auto copy_start = Clock::now();
  std::vector<Group> G(pp.G), H(pp.H);
  std::vector<Scalar> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  if (active_profile)
    active_profile->initial_allocation_copy_ms += elapsed_ms(copy_start);

  const auto transcript_start = Clock::now();
  Transcript transcript(pp, Z);
  if (active_profile)
    active_profile->transcript_challenge_ms += elapsed_ms(transcript_start);
  Proof proof;
  proof.rounds.reserve(pp.d);

  for (std::size_t round = 0; round < pp.d; ++round) {
    const std::size_t k = pp.d - round;
    const std::size_t half = xv.size() / 2;
    const auto GL = std::span(G).first(half), GR = std::span(G).subspan(half);
    const auto HL = std::span(H).first(half), HR = std::span(H).subspan(half);
    const auto xL = std::span(xv).first(half), xR = std::span(xv).subspan(half);
    const auto yL = std::span(yv).first(half), yR = std::span(yv).subspan(half);

    const Group A = add(add(msm(GL, xR), msm(HR, yL)),
                        mul(pp.K, inner_product(xR, yL)));
    const Group B = add(add(msm(GR, xL), msm(HL, yR)),
                        mul(pp.K, inner_product(xL, yR)));

    const auto challenge_start = Clock::now();
    const Scalar alpha = transcript.round(k, A, B);
    const Scalar alpha_inv = inverse(alpha);
    if (active_profile)
      active_profile->transcript_challenge_ms += elapsed_ms(challenge_start);

    const auto witness_start = Clock::now();
    auto xn = fold_scalars(xL, xR, alpha_inv);
    auto yn = fold_scalars(yL, yR, alpha);
    if (active_profile)
      active_profile->witness_folding_ms += elapsed_ms(witness_start);

    const auto generator_start = Clock::now();
    auto Gn = fold_groups(GL, GR, alpha);
    auto Hn = fold_groups(HL, HR, alpha_inv);
    if (active_profile)
      active_profile->generator_folding_ms += elapsed_ms(generator_start);
    if (trace) {
      trace->challenges.push_back(alpha);
      const Group before = relation(G, H, pp.K, xv, yv);
      const Group rhs = add(add(before, mul(A, alpha_inv)), mul(B, alpha));
      const Group after = relation(Gn, Hn, pp.K, xn, yn);
      trace->round_invariants.push_back(after == rhs);
    }
    proof.rounds.push_back({A, B});
    xv = std::move(xn); yv = std::move(yn);
    G = std::move(Gn); H = std::move(Hn);
  }
  proof.x_final = xv[0];
  proof.y_final = yv[0];
  return proof;
}

bool Verify(const PublicParams& pp, const Group& Z, const Proof& proof,
            Trace* trace) noexcept {
  try {
    const auto validation_start = Clock::now();
    const bool valid_pp = validate_public_params(pp);
    if (active_profile)
      active_profile->public_parameter_validation_ms += elapsed_ms(validation_start);
    const auto proof_validation_start = Clock::now();
    const bool valid_proof = valid_pp && valid_typed_proof_inputs(pp, Z, proof);
    if (active_profile)
      active_profile->proof_parsing_validation_ms += elapsed_ms(proof_validation_start);
    if (!valid_proof) return false;
    return verify_core(pp, Z, proof, trace);
  } catch (...) { return false; }
}

std::size_t proof_payload_bytes(std::size_t n) {
  const auto d = log2_exact(n);
  return 2 * d * group_bytes() + 2 * scalar_bytes();
}

std::size_t proof_wire_bytes(std::size_t n) {
  return 26 + proof_payload_bytes(n);
}

Bytes serialize_proof(const PublicParams& pp, const Proof& proof) {
  if (!validate_public_params(pp) || proof.rounds.size() != pp.d)
    throw std::invalid_argument("invalid proof shape");
  for (const auto& round : proof.rounds)
    if (!valid_group(round.A, false) || !valid_group(round.B, false))
      throw std::invalid_argument("invalid proof group element");
  Bytes out;
  out.reserve(proof_wire_bytes(pp.n));
  append_raw(out, kProofMagic);
  out.push_back(static_cast<std::uint8_t>(kProofVersion >> 8));
  out.push_back(static_cast<std::uint8_t>(kProofVersion));
  append_raw(out, u64be(pp.n));
  append_raw(out, u64be(pp.d));
  for (const auto& round : proof.rounds) {
    append_raw(out, serialize(round.A));
    append_raw(out, serialize(round.B));
  }
  append_raw(out, serialize(proof.x_final));
  append_raw(out, serialize(proof.y_final));
  return out;
}

bool deserialize_proof(const PublicParams& pp, std::span<const std::uint8_t> bytes,
                       Proof& proof) noexcept {
  try {
    const auto pp_validation_start = Clock::now();
    const bool valid_pp = validate_public_params(pp);
    if (active_profile)
      active_profile->public_parameter_validation_ms += elapsed_ms(pp_validation_start);
    const auto proof_start = Clock::now();
    if (!valid_pp || bytes.size() != proof_wire_bytes(pp.n) ||
        !std::equal(kProofMagic.begin(), kProofMagic.end(), bytes.begin()) ||
        bytes[8] != static_cast<std::uint8_t>(kProofVersion >> 8) ||
        bytes[9] != static_cast<std::uint8_t>(kProofVersion) ||
        read_u64(bytes, 10) != pp.n || read_u64(bytes, 18) != pp.d)
      return false;
    const std::size_t gb = group_bytes(), sb = scalar_bytes();
    std::size_t pos = 26;
    Proof candidate;
    candidate.rounds.resize(pp.d);
    for (auto& round : candidate.rounds) {
      if (!canonical_group(bytes.subspan(pos, gb), round.A)) return false;
      pos += gb;
      if (!canonical_group(bytes.subspan(pos, gb), round.B)) return false;
      pos += gb;
    }
    if (!canonical_scalar(bytes.subspan(pos, sb), candidate.x_final)) return false;
    pos += sb;
    if (!canonical_scalar(bytes.subspan(pos, sb), candidate.y_final)) return false;
    proof = std::move(candidate);
    if (active_profile)
      active_profile->proof_parsing_validation_ms += elapsed_ms(proof_start);
    return true;
  } catch (...) { return false; }
}

bool VerifySerialized(const PublicParams& pp, const Group& Z,
                      std::span<const std::uint8_t> bytes) noexcept {
  Proof proof;
  return deserialize_proof(pp, bytes, proof) && Verify(pp, Z, proof);
}

namespace internal {

Digest bp_public_parameter_digest(const PublicParams& pp) {
  return public_params_digest(pp);
}

bool replay_bp_challenges_prevalidated(
    const PublicParams& pp, const Digest& parameter_digest,
    const Group& Z, const Proof& proof,
    std::vector<Scalar>& challenges) noexcept {
  try {
    challenges.clear();
    if (proof.rounds.size() != pp.d) return false;
    Transcript transcript(pp, Z, &parameter_digest);
    challenges.reserve(pp.d);
    for (std::size_t round = 0; round < pp.d; ++round)
      challenges.push_back(transcript.round(
          pp.d - round, proof.rounds[round].A, proof.rounds[round].B));
    return true;
  } catch (...) {
    challenges.clear();
    return false;
  }
}

bool replay_bp_challenges_prevalidated(
    const PublicParams& pp, const Group& Z, const Proof& proof,
    std::vector<Scalar>& challenges) noexcept {
  try {
    challenges.clear();
    if (proof.rounds.size() != pp.d) return false;
    Transcript transcript(pp, Z);
    challenges.reserve(pp.d);
    for (std::size_t round = 0; round < pp.d; ++round)
      challenges.push_back(transcript.round(
          pp.d - round, proof.rounds[round].A, proof.rounds[round].B));
    return true;
  } catch (...) {
    challenges.clear();
    return false;
  }
}

bool replay_bp_challenges(const PublicParams& pp, const Group& Z,
                          const Proof& proof,
                          std::vector<Scalar>& challenges) noexcept {
  try {
    challenges.clear();
    if (!validate_public_params(pp) ||
        !valid_typed_proof_inputs(pp, Z, proof))
      return false;
    return replay_bp_challenges_prevalidated(
        pp, Z, proof, challenges);
  } catch (...) {
    challenges.clear();
    return false;
  }
}

}

namespace benchmark_internal {

Profile& Profile::operator+=(const Profile& o) {
  total_ms += o.total_ms;
  public_parameter_validation_ms += o.public_parameter_validation_ms;
  proof_parsing_validation_ms += o.proof_parsing_validation_ms;
  initial_allocation_copy_ms += o.initial_allocation_copy_ms;
  msm_input_copy_ms += o.msm_input_copy_ms;
  proof_message_msm_ms += o.proof_message_msm_ms;
  generator_folding_ms += o.generator_folding_ms;
  witness_folding_ms += o.witness_folding_ms;
  transcript_challenge_ms += o.transcript_challenge_ms;
  verifier_z_update_ms += o.verifier_z_update_ms;
  terminal_verification_ms += o.terminal_verification_ms;
  return *this;
}

Profile& Profile::operator/=(double d) {
  total_ms /= d;
  public_parameter_validation_ms /= d;
  proof_parsing_validation_ms /= d;
  initial_allocation_copy_ms /= d;
  msm_input_copy_ms /= d;
  proof_message_msm_ms /= d;
  generator_folding_ms /= d;
  witness_folding_ms /= d;
  transcript_challenge_ms /= d;
  verifier_z_update_ms /= d;
  terminal_verification_ms /= d;
  return *this;
}

double Profile::accounted_ms() const {
  return public_parameter_validation_ms + proof_parsing_validation_ms +
         initial_allocation_copy_ms + proof_message_msm_ms +
         generator_folding_ms + witness_folding_ms +
         transcript_challenge_ms + verifier_z_update_ms +
         terminal_verification_ms;
}

double Profile::residual_overhead_ms() const {
  return std::max(0.0, total_ms - accounted_ms());
}

Proof ProveProfiled(const PublicParams& pp, const Group& Z,
                    std::span<const Scalar> x, std::span<const Scalar> y,
                    Profile& profile) {
  profile = {};
  active_profile = &profile;
  const auto start = Clock::now();
  try {
    Proof proof = Prove(pp, Z, x, y);
    profile.total_ms = elapsed_ms(start);
    active_profile = nullptr;
    return proof;
  } catch (...) {
    active_profile = nullptr;
    throw;
  }
}

bool VerifySerializedProfiled(const PublicParams& pp, const Group& Z,
                              std::span<const std::uint8_t> proof_bytes,
                              Profile& profile) {
  profile = {};
  active_profile = &profile;
  const auto start = Clock::now();
  const bool accepted = VerifySerialized(pp, Z, proof_bytes);
  profile.total_ms = elapsed_ms(start);
  active_profile = nullptr;
  return accepted;
}

bool PrevalidateVerification(const PublicParams& pp, const Group& Z,
                             const Proof& proof,
                             PrevalidatedVerification& output) noexcept {
  output = {};
  try {
    if (!validate_public_params(pp) || !valid_typed_proof_inputs(pp, Z, proof))
      return false;
    output = {&pp, &Z, &proof};
    return true;
  } catch (...) { return false; }
}

bool VerifyPrevalidatedProfiled(const PrevalidatedVerification& inputs,
                                Profile& profile) noexcept {
  profile = {};
  if (!inputs.pp || !inputs.Z || !inputs.proof) return false;
  active_profile = &profile;
  const auto start = Clock::now();
  try {
    const bool accepted = verify_core(*inputs.pp, *inputs.Z, *inputs.proof, nullptr);
    profile.total_ms = elapsed_ms(start);
    active_profile = nullptr;
    return accepted;
  } catch (...) {
    profile.total_ms = elapsed_ms(start);
    active_profile = nullptr;
    return false;
  }
}

}

}

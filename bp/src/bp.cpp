#include "bp/bp.hpp"
#include "bp_internal.hpp"
#include "benchmark_internal.hpp"
#include "internal/bp_transcript.hpp"
#include "internal/protocol_utils.hpp"
#include "internal/protocol_validation.hpp"

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

constexpr std::array<std::uint8_t, 8> kProofMagic{
    'B','P','I','P','A','G','2','1'};
constexpr std::uint16_t kProofVersion = 1;

using internal::u64be;

void append_raw(Bytes& out, std::span<const std::uint8_t> field) {
  internal::append_raw(out, field);
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

std::size_t log2_exact(std::size_t n) {
  return internal::exact_log2(n);
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
  return internal::canonical_scalar(bytes, out);
}

bool canonical_group(std::span<const std::uint8_t> bytes, Group& out) {
  return internal::canonical_group(bytes, out);
}

std::uint64_t read_u64(std::span<const std::uint8_t> bytes, std::size_t pos) {
  std::uint64_t out = 0;
  if (!internal::read_u64(bytes, pos, out))
    throw std::out_of_range("truncated u64");
  return out;
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
  internal::BpTranscript transcript(pp, Z);
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
    pp.G.push_back(internal::derive_bp_generator(
        "BP-IPA-G-v1", public_seed, i, true));
    pp.H.push_back(internal::derive_bp_generator(
        "BP-IPA-H-v1", public_seed, i, true));
  }
  pp.K = internal::derive_bp_generator(
      "BP-IPA-K-v1", public_seed, 0, false);
  return pp;
}

bool validate_public_params(const PublicParams& pp) noexcept {
  return internal::validate_bp_public_params(pp);
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
  internal::BpTranscript transcript(pp, Z);
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
    const bool valid_proof =
        valid_pp && internal::validate_bp_proof_shape(pp, Z, proof);
    if (active_profile)
      active_profile->proof_parsing_validation_ms += elapsed_ms(proof_validation_start);
    if (!valid_proof) return false;
    return verify_core(pp, Z, proof, trace);
  } catch (...) { return false; }
}

std::size_t proof_payload_bytes(std::size_t n) {
  const auto d = log2_exact(n);
  std::size_t rounds = 0, finals = 0, total = 0;
  if (!internal::checked_mul(d, 2, rounds) ||
      !internal::checked_mul(rounds, group_bytes(), rounds) ||
      !internal::checked_mul(2, scalar_bytes(), finals) ||
      !internal::checked_add(rounds, finals, total))
    throw std::overflow_error("proof size overflow");
  return total;
}

std::size_t proof_wire_bytes(std::size_t n) {
  std::size_t total = 0;
  if (!internal::checked_add(26, proof_payload_bytes(n), total))
    throw std::overflow_error("proof wire size overflow");
  return total;
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
  if (out.size() != proof_wire_bytes(pp.n))
    throw std::runtime_error("unexpected proof encoding size");
  return out;
}

bool deserialize_proof(const PublicParams& pp, std::span<const std::uint8_t> bytes,
                       Proof& proof) noexcept {
  proof = {};
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
    pos += sb;
    if (pos != bytes.size()) return false;
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
  return bp_transcript_parameter_digest(pp);
}

bool replay_bp_challenges_prevalidated(
    const PublicParams& pp, const Digest& parameter_digest,
    const Group& Z, const Proof& proof,
    std::vector<Scalar>& challenges) noexcept {
  try {
    challenges.clear();
    if (proof.rounds.size() != pp.d) return false;
    BpTranscript transcript(pp, Z, &parameter_digest);
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
    BpTranscript transcript(pp, Z);
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
        !internal::validate_bp_proof_shape(pp, Z, proof))
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
    if (!validate_public_params(pp) ||
        !internal::validate_bp_proof_shape(pp, Z, proof))
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

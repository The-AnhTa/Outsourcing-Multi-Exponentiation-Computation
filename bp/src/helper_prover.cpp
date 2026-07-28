#include "bp/helper_prover.hpp"
#include "hp_protocol_internal.hpp"
#include "hp_transcript_internal.hpp"
#include "hp_vme_internal.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <limits>
#include <random>
#include <stdexcept>
#include <string_view>

namespace bp {
namespace {

using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(Clock::now() - start).count();
}

Bytes u64be(std::uint64_t value) {
  Bytes out;
  for (int shift = 56; shift >= 0; shift -= 8)
    out.push_back(static_cast<std::uint8_t>(value >> shift));
  return out;
}

void raw(Bytes& out, std::span<const std::uint8_t> field) {
  out.insert(out.end(), field.begin(), field.end());
}

void frame(Bytes& out, std::span<const std::uint8_t> field) {
  raw(out, u64be(field.size()));
  raw(out, field);
}

void frame(Bytes& out, std::string_view field) {
  frame(out, {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

bool power_of_two(std::size_t n) { return n && !(n & (n - 1)); }

std::size_t exact_log2(std::size_t n) {
  if (!power_of_two(n)) throw std::invalid_argument("n must be a power of two");
  std::size_t d = 0;
  while (n > 1) { n >>= 1; ++d; }
  return d;
}

Scalar fadd(const Scalar& a, const Scalar& b) {
  Scalar out; Scalar::add(out, a, b); return out;
}
Scalar fmul(const Scalar& a, const Scalar& b) {
  Scalar out; Scalar::mul(out, a, b); return out;
}
Scalar finv(const Scalar& a) {
  if (a.isZero()) throw std::invalid_argument("zero scalar");
  Scalar out; Scalar::inv(out, a); return out;
}
Scalar fneg(const Scalar& a) {
  Scalar out; Scalar::neg(out, a); return out;
}
Group gadd(const Group& a, const Group& b) {
  Group out; Group::add(out, a, b); return out;
}
Group gmul(const Group& p, const Scalar& x) {
  Group out; Group::mul(out, p, x); return out;
}
Scalar inner(std::span<const Scalar> a, std::span<const Scalar> b) {
  if (a.size() != b.size()) throw std::invalid_argument("inner length");
  Scalar out; out.clear();
  for (std::size_t i = 0; i < a.size(); ++i)
    out = fadd(out, fmul(a[i], b[i]));
  return out;
}
Group msm(std::span<const Group> points, std::span<const Scalar> scalars) {
  if (points.size() != scalars.size()) throw std::invalid_argument("MSM length");
  Group out; out.clear();
  if (!points.empty()) {
    std::vector<Group> work(points.begin(), points.end());
    Group::mulVec(out, work.data(), scalars.data(), work.size());
  }
  return out;
}
bool valid_group(const Group& p, bool nonidentity = false) {
  return p.isValid() && p.isValidOrder() && (!nonidentity || !p.isZero());
}

Group derive_g2(std::string_view domain, std::span<const std::uint8_t> seed,
                std::size_t index, bool indexed) {
  Bytes input;
  frame(input, domain);
  frame(input, seed);
  if (indexed) frame(input, u64be(index));
  Group out;
  mcl::bn::hashAndMapToG2(out, input.data(), input.size());
  if (!valid_group(out, true)) throw std::runtime_error("invalid derived G2");
  return out;
}

Digest hp_crs_digest(const HpPublicParams& pp) {
  Bytes input;
  frame(input, "BPVME/HP/CRS/v1");
  frame(input, pp.bp.group_identifier);
  frame(input, pp.bp.scalar_modulus);
  frame(input, u64be(pp.bp.n));
  frame(input, u64be(pp.N));
  for (const auto& p : pp.bp.G) frame(input, serialize(p));
  for (const auto& p : pp.bp.H) frame(input, serialize(p));
  frame(input, serialize(pp.bp.K));
  frame(input, pp.vme.digest);
  frame(input, pp.bp.hash_suite_identifier);
  return sha256(input);
}

Digest z_gamma_digest(std::span<const Scalar> z) {
  Bytes input;
  frame(input, "BPVME/HP/Z-GAMMA/v1");
  for (const auto& value : z) frame(input, serialize(value));
  return sha256(input);
}

Digest application_context(const HpPublicParams& pp, const Group& Z,
                           std::span<const RoundMessage> rounds,
                           const Scalar& x0, const Scalar& y0,
                           const Scalar& gamma,
                           std::span<const Scalar> z_gamma) {
  Bytes input;
  frame(input, "BPVME/HP/VME-CONTEXT/v1");
  frame(input, pp.crs_digest);
  frame(input, u64be(pp.bp.n));
  frame(input, serialize(Z));
  for (const auto& round : rounds) {
    frame(input, serialize(round.A));
    frame(input, serialize(round.B));
  }
  frame(input, serialize(x0));
  frame(input, serialize(y0));
  frame(input, serialize(gamma));
  frame(input, z_gamma_digest(z_gamma));
  return sha256(input);
}

bool fast_hp_binding(const HpPublicParams& pp, const Digest& precomp_crs,
                     const VmePrecomputation& precomp) {
  try {
    if (!power_of_two(pp.bp.n) || pp.bp.d != exact_log2(pp.bp.n) ||
        pp.bp.n > std::numeric_limits<std::size_t>::max() / 2 ||
        pp.N != 2 * pp.bp.n || pp.vme.dimension != pp.N ||
        pp.vme.log_dimension != pp.bp.d + 1 ||
        pp.bp.G.size() != pp.bp.n || pp.bp.H.size() != pp.bp.n ||
        pp.vme.auxiliary_G.size() != pp.N ||
        pp.bp.transcript_domain != kHpTranscriptDomain ||
        !pp.bp.transcript_crs_digest ||
        *pp.bp.transcript_crs_digest != pp.crs_digest ||
        precomp_crs != pp.crs_digest ||
        pp.vme.fixed_P.size() != pp.N)
      return false;
    for (std::size_t i = 0; i < pp.bp.n; ++i)
      if (pp.vme.fixed_P[i] != pp.bp.G[i] ||
          pp.vme.fixed_P[pp.bp.n + i] != pp.bp.H[i])
        return false;
    return hp_crs_digest(pp) == pp.crs_digest &&
           hp_internal::validate_vme_precomputation_binding(pp.vme, precomp);
  } catch (...) { return false; }
}

struct GeneratedRounds {
  std::vector<RoundMessage> rounds;
  std::vector<Scalar> alphas;
  Scalar x0, y0;
};

std::vector<Scalar> fold_scalars(std::span<const Scalar> left,
                                 std::span<const Scalar> right,
                                 const Scalar& factor) {
  std::vector<Scalar> out;
  out.reserve(left.size());
  for (std::size_t i = 0; i < left.size(); ++i)
    out.push_back(fadd(left[i], fmul(factor, right[i])));
  return out;
}

std::vector<Group> fold_groups(std::span<const Group> left,
                               std::span<const Group> right,
                               const Scalar& factor) {
  std::vector<Group> out;
  out.reserve(left.size());
  for (std::size_t i = 0; i < left.size(); ++i)
    out.push_back(gadd(left[i], gmul(right[i], factor)));
  return out;
}

GeneratedRounds compute_rounds(const HpPublicParams& pp, const Group& Z,
                               std::span<const Scalar> x,
                               std::span<const Scalar> y) {
  std::vector<Group> G(pp.bp.G), H(pp.bp.H);
  std::vector<Scalar> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  hp_internal::HpBpTranscript transcript(pp.bp, Z);
  GeneratedRounds out;
  out.rounds.reserve(pp.bp.d); out.alphas.reserve(pp.bp.d);
  for (std::size_t round_index = 0; round_index < pp.bp.d; ++round_index) {
    const std::size_t k = pp.bp.d - round_index, half = xv.size() / 2;
    const auto GL = std::span(G).first(half), GR = std::span(G).subspan(half);
    const auto HL = std::span(H).first(half), HR = std::span(H).subspan(half);
    const auto xL = std::span(xv).first(half), xR = std::span(xv).subspan(half);
    const auto yL = std::span(yv).first(half), yR = std::span(yv).subspan(half);
    const Scalar cA = inner(xR, yL), cB = inner(xL, yR);
    const Group A = gadd(gadd(msm(GL, xR), msm(HR, yL)), gmul(pp.bp.K, cA));
    const Group B = gadd(gadd(msm(GR, xL), msm(HL, yR)), gmul(pp.bp.K, cB));
    const Scalar alpha = transcript.round(k, A, B), alpha_inv = finv(alpha);
    out.rounds.push_back({A, B}); out.alphas.push_back(alpha);
    xv = fold_scalars(xL, xR, alpha_inv);
    yv = fold_scalars(yL, yR, alpha);
    G = fold_groups(GL, GR, alpha);
    H = fold_groups(HL, HR, alpha_inv);
  }
  out.x0 = xv[0]; out.y0 = yv[0];
  return out;
}

struct ReplayResult {
  std::vector<Scalar> alphas;
  Scalar x0, y0, gamma;
};

ReplayResult replay_rounds(const HpPublicParams& pp, const Group& Z,
                           std::span<const Scalar> x,
                           std::span<const Scalar> y,
                           std::span<const RoundMessage> rounds) {
  if (rounds.size() != pp.bp.d) throw std::invalid_argument("round count");
  hp_internal::HpBpTranscript transcript(pp.bp, Z);
  std::vector<Scalar> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  ReplayResult out;
  out.alphas.reserve(pp.bp.d);
  for (std::size_t round_index = 0; round_index < pp.bp.d; ++round_index) {
    const std::size_t k = pp.bp.d - round_index, half = xv.size() / 2;
    const Scalar alpha =
        transcript.round(k, rounds[round_index].A, rounds[round_index].B);
    const Scalar alpha_inv = finv(alpha);
    out.alphas.push_back(alpha);
    const auto xL = std::span(xv).first(half), xR = std::span(xv).subspan(half);
    const auto yL = std::span(yv).first(half), yR = std::span(yv).subspan(half);
    xv = fold_scalars(xL, xR, alpha_inv);
    yv = fold_scalars(yL, yR, alpha);
  }
  out.x0 = xv[0]; out.y0 = yv[0];
  out.gamma = transcript.outer_gamma(out.x0, out.y0);
  return out;
}

hp_internal::HpAggregate build_aggregate(
    const HpPublicParams& pp, std::span<const Scalar> x,
    std::span<const Scalar> y, std::span<const RoundMessage> rounds,
    std::span<const Scalar> alphas, const Scalar& gamma) {
  if (rounds.size() != pp.bp.d || alphas.size() != pp.bp.d || gamma.isZero())
    throw std::invalid_argument("aggregate shape");
  hp_internal::HpAggregate out;
  out.z_gamma.resize(pp.N);
  for (auto& value : out.z_gamma) value.clear();
  out.c_gamma.clear();
  out.X_gamma.clear();
  std::vector<Scalar> xv(x.begin(), x.end()), yv(y.begin(), y.end());
  std::vector<Scalar> rhoG(1), rhoH(1);
  rhoG[0] = 1; rhoH[0] = 1;
  Scalar gamma_power = gamma;

  for (std::size_t round_index = 0; round_index < pp.bp.d; ++round_index) {
    const std::size_t m = xv.size(), half = m / 2;
    const auto xL = std::span(xv).first(half), xR = std::span(xv).subspan(half);
    const auto yL = std::span(yv).first(half), yR = std::span(yv).subspan(half);
    const Scalar weightA = gamma_power;
    const Scalar weightB = fmul(gamma_power, gamma);
    for (std::size_t r = 0; r < rhoG.size(); ++r) {
      const std::size_t gbase = r * m;
      const std::size_t hbase = pp.bp.n + r * m;
      for (std::size_t i = 0; i < half; ++i) {
        out.z_gamma[gbase + i] = fadd(
            out.z_gamma[gbase + i], fmul(weightA, fmul(rhoG[r], xR[i])));
        out.z_gamma[gbase + half + i] = fadd(
            out.z_gamma[gbase + half + i],
            fmul(weightB, fmul(rhoG[r], xL[i])));
        out.z_gamma[hbase + half + i] = fadd(
            out.z_gamma[hbase + half + i],
            fmul(weightA, fmul(rhoH[r], yL[i])));
        out.z_gamma[hbase + i] = fadd(
            out.z_gamma[hbase + i], fmul(weightB, fmul(rhoH[r], yR[i])));
      }
    }
    const Scalar cA = inner(xR, yL), cB = inner(xL, yR);
    out.c_gamma = fadd(out.c_gamma,
        fadd(fmul(weightA, cA), fmul(weightB, cB)));
    out.X_gamma = gadd(out.X_gamma,
        gadd(gmul(rounds[round_index].A, weightA),
             gmul(rounds[round_index].B, weightB)));
    const Scalar alpha = alphas[round_index], alpha_inv = finv(alpha);
    xv = fold_scalars(xL, xR, alpha_inv);
    yv = fold_scalars(yL, yR, alpha);
    std::vector<Scalar> nextG, nextH;
    nextG.reserve(rhoG.size() * 2); nextH.reserve(rhoH.size() * 2);
    for (const auto& value : rhoG) {
      nextG.push_back(value); nextG.push_back(fmul(value, alpha));
    }
    for (const auto& value : rhoH) {
      nextH.push_back(value); nextH.push_back(fmul(value, alpha_inv));
    }
    rhoG = std::move(nextG); rhoH = std::move(nextH);
    gamma_power = fmul(weightB, gamma);
  }
  out.X_gamma = gadd(out.X_gamma, gmul(pp.bp.K, fneg(out.c_gamma)));
  out.x0 = xv[0]; out.y0 = yv[0];
  return out;
}

bool canonical_group(std::span<const std::uint8_t> bytes, Group& out) {
  Group value;
  if (bytes.size() != group_bytes() ||
      value.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
      !valid_group(value) || serialize(value) != Bytes(bytes.begin(), bytes.end()))
    return false;
  out = value;
  return true;
}

std::uint64_t read_u64(std::span<const std::uint8_t> bytes, std::size_t pos) {
  std::uint64_t out = 0;
  for (std::size_t i = 0; i < 8; ++i) out = (out << 8) | bytes[pos + i];
  return out;
}

constexpr std::array<std::uint8_t, 8> kHpMagic{
    'B','P','H','P','G','2','0','1'};
constexpr std::uint16_t kHpVersion = 1;
constexpr std::size_t kHpHeader = 27;

}

HpSetupResult setup_hp(std::size_t n, std::span<const std::uint8_t> seed,
                       HpSetupTimings* timings) {
  initialize();
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  const std::size_t d = exact_log2(n);
  if (n > std::numeric_limits<std::size_t>::max() / 2)
    throw std::overflow_error("N overflow");
  HpSetupResult out;
  out.pp.bp.n = n; out.pp.bp.d = d;
  out.pp.bp.transcript_domain = kHpTranscriptDomain;
  out.pp.bp.G.reserve(n); out.pp.bp.H.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    out.pp.bp.G.push_back(derive_g2("BPVME/HP/BP-G/v1", seed, i, true));
    out.pp.bp.H.push_back(derive_g2("BPVME/HP/BP-H/v1", seed, i, true));
  }
  out.pp.bp.K = derive_g2("BPVME/HP/BP-K/v1", seed, 0, false);
  out.pp.N = 2 * n;
  std::vector<Group> P;
  P.reserve(out.pp.N);
  P.insert(P.end(), out.pp.bp.G.begin(), out.pp.bp.G.end());
  P.insert(P.end(), out.pp.bp.H.begin(), out.pp.bp.H.end());
  Bytes vme_seed;
  frame(vme_seed, "BPVME/HP/VME/v1");
  frame(vme_seed, seed);
  out.pp.vme = hp_internal::setup_vme(P, vme_seed);
  out.pp.crs_digest = hp_crs_digest(out.pp);
  out.pp.bp.transcript_crs_digest = out.pp.crs_digest;
  const auto precompute_start = Clock::now();
  const VmePrecomputation pre = hp_internal::precompute_vme(out.pp.vme);
  out.helper_precomp = {pre, out.pp.crs_digest};
  out.client_precomp = {pre, out.pp.crs_digest};
  if (timings) {
    timings->precomputation_ms = milliseconds(precompute_start);
    timings->setup_ms = milliseconds(total_start) - timings->precomputation_ms;
  }
  return out;
}

HpProof prove_hp(const HpPublicParams& pp, const HpHelperPrecomp& pre,
                 const Group& Z, std::span<const Scalar> x,
                 std::span<const Scalar> y, HpProveTimings* timings) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  if (!fast_hp_binding(pp, pre.crs_digest, pre.vme) || !valid_group(Z) ||
      x.size() != pp.bp.n || y.size() != pp.bp.n)
    throw std::invalid_argument("invalid helper input");
  if (pp.bp.d == 0) {
    if (timings) timings->total_ms = milliseconds(total_start);
    return {};
  }
  const auto rounds_start = Clock::now();
  const GeneratedRounds generated = compute_rounds(pp, Z, x, y);
  if (timings) timings->bp_round_generation_ms = milliseconds(rounds_start);
  hp_internal::HpBpTranscript transcript(pp.bp, Z);
  for (std::size_t i = 0; i < generated.rounds.size(); ++i)
    (void)transcript.round(pp.bp.d - i, generated.rounds[i].A,
                           generated.rounds[i].B);
  const Scalar gamma = transcript.outer_gamma(generated.x0, generated.y0);
  const auto aggregate_start = Clock::now();
  const auto aggregate = build_aggregate(
      pp, x, y, generated.rounds, generated.alphas, gamma);
  const Digest context = application_context(
      pp, Z, generated.rounds, aggregate.x0, aggregate.y0, gamma,
      aggregate.z_gamma);
  if (timings) timings->outer_aggregation_ms = milliseconds(aggregate_start);
  const auto vme_start = Clock::now();
  VmeProof vme = hp_internal::prove_vme(
      pp.vme, pre.vme, context, aggregate.X_gamma, aggregate.z_gamma);
  if (timings) {
    timings->vme_prove_ms = milliseconds(vme_start);
    timings->total_ms = milliseconds(total_start);
  }
  return {generated.rounds, std::move(vme)};
}

static std::optional<Proof> verify_hp_impl(
    const HpPublicParams& pp, const HpClientPrecomp& pre, const Group& Z,
    std::span<const Scalar> x, std::span<const Scalar> y,
    const HpProof& helper, HpVerifyTimings* timings,
    bool proof_prevalidated) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  const auto validation_start = Clock::now();
  if (!fast_hp_binding(pp, pre.crs_digest, pre.vme) || !valid_group(Z) ||
      x.size() != pp.bp.n || y.size() != pp.bp.n ||
      helper.rounds.size() != pp.bp.d) {
    if (timings) {
      timings->parse_and_validation_ms = milliseconds(validation_start);
      timings->total_ms = milliseconds(total_start);
    }
    return std::nullopt;
  }
  if (!proof_prevalidated)
    for (const auto& round : helper.rounds)
      if (!valid_group(round.A) || !valid_group(round.B)) return std::nullopt;
  if (pp.bp.d == 0) {
    if (helper.vme_proof) return std::nullopt;
    Proof result{{}, x[0], y[0]};
    if (timings) {
      timings->parse_and_validation_ms = milliseconds(validation_start);
      timings->total_ms = milliseconds(total_start);
    }
    return result;
  }
  if (!helper.vme_proof ||
      (!proof_prevalidated &&
       !hp_internal::validate_vme_proof(pp.vme, *helper.vme_proof)))
    return std::nullopt;
  if (timings)
    timings->parse_and_validation_ms = milliseconds(validation_start);
  const auto replay_start = Clock::now();
  ReplayResult replay;
  try { replay = replay_rounds(pp, Z, x, y, helper.rounds); }
  catch (...) { return std::nullopt; }
  if (timings)
    timings->bp_transcript_replay_ms = milliseconds(replay_start);
  const auto aggregate_start = Clock::now();
  hp_internal::HpAggregate aggregate;
  Digest context;
  try {
    aggregate = build_aggregate(
        pp, x, y, helper.rounds, replay.alphas, replay.gamma);
    if (aggregate.x0 != replay.x0 || aggregate.y0 != replay.y0)
      return std::nullopt;
    context = application_context(pp, Z, helper.rounds, replay.x0, replay.y0,
                                  replay.gamma, aggregate.z_gamma);
  } catch (...) { return std::nullopt; }
  if (timings)
    timings->outer_aggregation_ms = milliseconds(aggregate_start);
  const auto vme_start = Clock::now();
  hp_internal::VmeVerificationStats vme_stats;
  const bool accepted = hp_internal::verify_vme_prevalidated_optimized(
      pp.vme, pre.vme, context, aggregate.X_gamma, aggregate.z_gamma,
      *helper.vme_proof, &vme_stats);
  if (timings) {
    timings->vme_verify_ms = milliseconds(vme_start);
    timings->vme_transcript_ms = vme_stats.transcript_ms;
    timings->vme_batch_inversion_ms = vme_stats.batch_inversion_ms;
    timings->vme_symbolic_recurrence_ms = vme_stats.recurrence_ms;
    timings->vme_terminal_assembly_ms = vme_stats.terminal_ms;
    timings->vme_gt_multiexp_ms = vme_stats.gt_multiexp_ms;
    timings->vme_multi_pairing_ms = vme_stats.multi_pairing_ms;
    timings->total_ms = milliseconds(total_start);
  }
  if (!accepted) return std::nullopt;
  return Proof{helper.rounds, replay.x0, replay.y0};
}

std::optional<Proof> verify_hp(const HpPublicParams& pp,
                               const HpClientPrecomp& pre, const Group& Z,
                               std::span<const Scalar> x,
                               std::span<const Scalar> y,
                               const HpProof& helper,
                               HpVerifyTimings* timings) {
  return verify_hp_impl(pp, pre, Z, x, y, helper, timings, false);
}

std::optional<Proof> verify_hp_serialized(
    const HpPublicParams& pp, const HpClientPrecomp& pre, const Group& Z,
    std::span<const Scalar> x, std::span<const Scalar> y,
    std::span<const std::uint8_t> bytes, HpVerifyTimings* timings) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  const auto parse_start = Clock::now();
  HpProof proof;
  if (!deserialize_hp_proof(pp, bytes, proof)) {
    if (timings) {
      timings->parse_and_validation_ms = milliseconds(parse_start);
      timings->total_ms = milliseconds(total_start);
    }
    return std::nullopt;
  }
  const double parse_ms = milliseconds(parse_start);
  HpVerifyTimings inner_timings;
  auto result =
      verify_hp_impl(pp, pre, Z, x, y, proof, &inner_timings, true);
  if (timings) {
    *timings = inner_timings;
    timings->parse_and_validation_ms += parse_ms;
    timings->total_ms = milliseconds(total_start);
  }
  return result;
}

HpInstance generate_hp_instance(const HpPublicParams& pp) {
  if (!power_of_two(pp.bp.n)) throw std::invalid_argument("invalid HP parameters");
  HpInstance out;
  out.x.resize(pp.bp.n); out.y.resize(pp.bp.n);
  for (std::size_t i = 0; i < pp.bp.n; ++i) {
    out.x[i].setByCSPRNG();
    out.y[i].setByCSPRNG();
  }
  out.Z = commit(pp.bp, out.x, out.y);
  return out;
}

std::size_t vme_proof_payload_bytes(std::size_t n) {
  const std::size_t d = exact_log2(n), D = d + 1;
  if (d == 0) return 0;
  return (4 * (D - 1) + 7 * D) * hp_internal::gt_bytes() +
         2 * hp_internal::g1_bytes() + group_bytes();
}

std::size_t hp_proof_payload_bytes(std::size_t n) {
  const std::size_t d = exact_log2(n);
  return 2 * d * group_bytes() + (d == 0 ? 0 : vme_proof_payload_bytes(n));
}

std::size_t hp_proof_wire_bytes(std::size_t n) {
  return kHpHeader + hp_proof_payload_bytes(n);
}

Bytes serialize_hp_proof(const HpPublicParams& pp, const HpProof& proof) {
  if (proof.rounds.size() != pp.bp.d ||
      (pp.bp.d == 0 ? proof.vme_proof.has_value() : !proof.vme_proof.has_value()))
    throw std::invalid_argument("HP proof shape");
  for (const auto& round : proof.rounds)
    if (!valid_group(round.A) || !valid_group(round.B))
      throw std::invalid_argument("invalid HP round element");
  if (proof.vme_proof &&
      !hp_internal::validate_vme_proof(pp.vme, *proof.vme_proof))
    throw std::invalid_argument("invalid VME proof");
  Bytes out;
  out.reserve(hp_proof_wire_bytes(pp.bp.n));
  raw(out, kHpMagic);
  out.push_back(static_cast<std::uint8_t>(kHpVersion >> 8));
  out.push_back(static_cast<std::uint8_t>(kHpVersion));
  raw(out, u64be(pp.bp.n)); raw(out, u64be(pp.bp.d));
  out.push_back(proof.vme_proof ? 1 : 0);
  for (const auto& round : proof.rounds) {
    raw(out, serialize(round.A)); raw(out, serialize(round.B));
  }
  if (proof.vme_proof) {
    const auto& v = *proof.vme_proof;
    for (const auto& c : v.rexp_claims)
      for (const GT* field : {&c.E, &c.F, &c.TL, &c.TR})
        raw(out, hp_internal::serialize_gt(*field));
    raw(out, hp_internal::serialize_g1(v.R));
    for (const auto& f : v.dory_folds)
      for (const GT* field : {&f.D1L, &f.D1R, &f.D2L, &f.D2R, &f.W1, &f.W2})
        raw(out, hp_internal::serialize_gt(*field));
    for (const auto& U : v.batch_U) raw(out, hp_internal::serialize_gt(U));
    raw(out, hp_internal::serialize_g1(v.phi_final));
    raw(out, serialize(v.theta_final));
  }
  return out;
}

bool deserialize_hp_proof(const HpPublicParams& pp,
                          std::span<const std::uint8_t> bytes,
                          HpProof& proof) noexcept {
  try {
    if (bytes.size() != hp_proof_wire_bytes(pp.bp.n) ||
        !std::equal(kHpMagic.begin(), kHpMagic.end(), bytes.begin()) ||
        bytes[8] != static_cast<std::uint8_t>(kHpVersion >> 8) ||
        bytes[9] != static_cast<std::uint8_t>(kHpVersion) ||
        read_u64(bytes, 10) != pp.bp.n || read_u64(bytes, 18) != pp.bp.d ||
        bytes[26] != static_cast<std::uint8_t>(pp.bp.d == 0 ? 0 : 1))
      return false;
    std::size_t pos = kHpHeader;
    HpProof candidate;
    candidate.rounds.resize(pp.bp.d);
    for (auto& round : candidate.rounds) {
      if (!canonical_group(bytes.subspan(pos, group_bytes()), round.A)) return false;
      pos += group_bytes();
      if (!canonical_group(bytes.subspan(pos, group_bytes()), round.B)) return false;
      pos += group_bytes();
    }
    if (pp.bp.d > 0) {
      VmeProof v;
      const std::size_t D = pp.bp.d + 1, gtb = hp_internal::gt_bytes();
      v.rexp_claims.resize(D - 1);
      for (auto& c : v.rexp_claims)
        for (GT* field : {&c.E, &c.F, &c.TL, &c.TR}) {
          if (!hp_internal::deserialize_gt(bytes.subspan(pos, gtb), *field))
            return false;
          pos += gtb;
        }
      if (!hp_internal::deserialize_g1(
              bytes.subspan(pos, hp_internal::g1_bytes()), v.R))
        return false;
      pos += hp_internal::g1_bytes();
      v.dory_folds.resize(D);
      for (auto& f : v.dory_folds)
        for (GT* field : {&f.D1L, &f.D1R, &f.D2L, &f.D2R, &f.W1, &f.W2}) {
          if (!hp_internal::deserialize_gt(bytes.subspan(pos, gtb), *field))
            return false;
          pos += gtb;
        }
      v.batch_U.resize(D);
      for (auto& U : v.batch_U) {
        if (!hp_internal::deserialize_gt(bytes.subspan(pos, gtb), U)) return false;
        pos += gtb;
      }
      if (!hp_internal::deserialize_g1(
              bytes.subspan(pos, hp_internal::g1_bytes()), v.phi_final))
        return false;
      pos += hp_internal::g1_bytes();
      if (!canonical_group(bytes.subspan(pos, group_bytes()), v.theta_final))
        return false;
      pos += group_bytes();
      candidate.vme_proof = std::move(v);
    }
    if (pos != bytes.size()) return false;
    proof = std::move(candidate);
    return true;
  } catch (...) { return false; }
}

namespace hp_internal {

HpRoundTrace reconstruct_hp_trace(const HpPublicParams& pp, const Group& Z,
                                  std::span<const Scalar> x,
                                  std::span<const Scalar> y,
                                  std::span<const RoundMessage> rounds) {
  const ReplayResult replay = replay_rounds(pp, Z, x, y, rounds);
  HpRoundTrace trace;
  trace.alphas = replay.alphas;
  trace.x0 = replay.x0; trace.y0 = replay.y0; trace.gamma = replay.gamma;
  trace.aggregate =
      build_aggregate(pp, x, y, rounds, replay.alphas, replay.gamma);
  trace.application_context = application_context(
      pp, Z, rounds, replay.x0, replay.y0, replay.gamma,
      trace.aggregate.z_gamma);
  return trace;
}

Group direct_hp_aggregate_msm_for_test(const HpPublicParams& pp,
                                       const HpAggregate& aggregate) {
  return msm(pp.vme.fixed_P, aggregate.z_gamma);
}

bool verify_hp_reference_for_test(const HpPublicParams& pp,
                                  const HpClientPrecomp& pre, const Group& Z,
                                  std::span<const Scalar> x,
                                  std::span<const Scalar> y,
                                  const HpProof& proof) {
  try {
    if (!proof.vme_proof || proof.rounds.size() != pp.bp.d) return false;
    const auto trace = reconstruct_hp_trace(pp, Z, x, y, proof.rounds);
    return verify_vme_reference(
        pp.vme, pre.vme, trace.application_context, trace.aggregate.X_gamma,
        trace.aggregate.z_gamma, *proof.vme_proof);
  } catch (...) { return false; }
}

bool verify_hp_wrong_round_encoding_for_test(
    const HpPublicParams& pp, const HpClientPrecomp& pre, const Group& Z,
    std::span<const Scalar> x, std::span<const Scalar> y,
    const HpProof& proof) {
  try {
    if (!proof.vme_proof || proof.rounds.size() != pp.bp.d) return false;
    HpBpTranscript transcript(pp.bp, Z);
    std::vector<Scalar> xv(x.begin(), x.end()), yv(y.begin(), y.end()), alphas;
    for (std::size_t round_index = 0; round_index < pp.bp.d; ++round_index) {
      const std::size_t k = pp.bp.d - round_index;
      Bytes little_endian(8);
      for (std::size_t i = 0; i < 8; ++i)
        little_endian[i] = static_cast<std::uint8_t>(k >> (8 * i));
      const Scalar alpha = transcript.round_with_index_bytes_for_test(
          little_endian, proof.rounds[round_index].A,
          proof.rounds[round_index].B);
      const Scalar alpha_inv = finv(alpha);
      const std::size_t half = xv.size() / 2;
      const auto xL = std::span(xv).first(half), xR = std::span(xv).subspan(half);
      const auto yL = std::span(yv).first(half), yR = std::span(yv).subspan(half);
      xv = fold_scalars(xL, xR, alpha_inv);
      yv = fold_scalars(yL, yR, alpha);
      alphas.push_back(alpha);
    }
    const Scalar gamma = transcript.outer_gamma(xv[0], yv[0]);
    const auto aggregate =
        build_aggregate(pp, x, y, proof.rounds, alphas, gamma);
    const Digest context = application_context(
        pp, Z, proof.rounds, aggregate.x0, aggregate.y0, gamma,
        aggregate.z_gamma);
    return verify_vme_optimized(
        pp.vme, pre.vme, context, aggregate.X_gamma, aggregate.z_gamma,
        *proof.vme_proof);
  } catch (...) { return false; }
}

}

}

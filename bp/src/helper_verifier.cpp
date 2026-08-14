#include "bp/helper_verifier.hpp"

#include "bp_internal.hpp"
#include "hp_vme_internal.hpp"
#include "hv_internal.hpp"
#include "internal/protocol_utils.hpp"
#include "internal/protocol_validation.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <limits>
#include <random>
#include <stdexcept>
#include <string_view>

namespace bp {

namespace internal {

struct HvPreparedStatementCacheAccess {
  static bool load(const HvPreparedStatementCache& cache,
                   const HvPublicParams& pp, const Digest& context,
                   hv_internal::PreparedStatement& prepared) {
    if (!cache.ready_ || cache.crs_digest_ != pp.crs_digest ||
        cache.statement_context_ != context || cache.z_v_.size() != pp.N ||
        !valid_group(cache.X0_))
      return false;
    prepared.z_v = cache.z_v_;
    prepared.X0 = cache.X0_;
    prepared.context = cache.statement_context_;
    return true;
  }

  static void store(HvPreparedStatementCache& cache,
                    const HvPublicParams& pp,
                    const hv_internal::PreparedStatement& prepared) {
    cache.ready_ = true;
    cache.crs_digest_ = pp.crs_digest;
    cache.statement_context_ = prepared.context;
    cache.z_v_ = prepared.z_v;
    cache.X0_ = prepared.X0;
  }
};

}  

namespace {

using Clock = std::chrono::steady_clock;

double milliseconds(Clock::time_point start) {
  return std::chrono::duration<double, std::milli>(
      Clock::now() - start).count();
}

using internal::exact_log2;
using internal::fmul;
using internal::fneg;
using internal::frame;
using internal::gadd;
using internal::gmul;
using internal::msm;
using internal::power_of_two;
using internal::u64be;
using internal::valid_group;

void raw(Bytes& out, std::span<const std::uint8_t> value) {
  internal::append_raw(out, value);
}

Scalar fone() { Scalar out; out = 1; return out; }

Digest bp_parameter_digest(const PublicParams& pp) {
  Bytes input;
  frame(input, "BPVME/HV/BP-PARAMETERS/v1");
  frame(input, pp.group_identifier);
  frame(input, pp.scalar_modulus);
  frame(input, pp.hash_suite_identifier);
  frame(input, pp.transcript_domain);
  frame(input, u64be(pp.n));
  for (const auto& point : pp.G) frame(input, serialize(point));
  for (const auto& point : pp.H) frame(input, serialize(point));
  frame(input, serialize(pp.K));
  return sha256(input);
}

Digest p_digest(const HvPublicParams& pp) {
  Bytes input;
  frame(input, "BPVME/HV/P/v1");
  frame(input, u64be(pp.N));
  for (const auto& point : pp.bp.G) frame(input, serialize(point));
  for (const auto& point : pp.bp.H) frame(input, serialize(point));
  return sha256(input);
}

Digest hv_crs_digest(const HvPublicParams& pp) {
  Bytes input;
  frame(input, "BPVME/HV/CRS/v1");
  frame(input, kGroupIdentifier);
  frame(input, kScalarModulus);
  frame(input, kHashSuiteIdentifier);
  frame(input, kHvContextDomain);
  frame(input, kHvVmeDomain);
  frame(input, u64be(pp.bp.n));
  frame(input, u64be(pp.bp.d));
  frame(input, u64be(pp.N));
  for (const auto& point : pp.bp.G) frame(input, serialize(point));
  for (const auto& point : pp.bp.H) frame(input, serialize(point));
  frame(input, serialize(pp.bp.K));
  frame(input, pp.vme.digest);
  frame(input, pp.bp_parameter_digest);
  frame(input, pp.bp_transcript_parameter_digest);
  frame(input, pp.P_digest);
  return sha256(input);
}

Digest hv_context(const HvPublicParams& pp, const HvStatement& statement) {
  Bytes input;
  frame(input, kHvContextDomain);
  frame(input, pp.crs_digest);
  frame(input, u64be(pp.N));
  frame(input, serialize(statement.Z));
  frame(input, u64be(statement.bulletproof.rounds.size()));
  for (const auto& round : statement.bulletproof.rounds) {
    frame(input, serialize(round.A));
    frame(input, serialize(round.B));
  }
  frame(input, serialize(statement.bulletproof.x_final));
  frame(input, serialize(statement.bulletproof.y_final));
  return sha256(input);
}

bool valid_statement(const HvPublicParams& pp, const HvStatement& statement) {
  return internal::validate_hv_statement_shape(pp, statement);
}

bool fast_binding(const HvPublicParams& pp, const VmePrecomputation& vme,
                  const Digest& precomp_crs) {
  try {
    if (!internal::validate_hv_public_params_shape(pp) ||
        pp.bp_parameter_digest != bp_parameter_digest(pp.bp) ||
        pp.P_digest != p_digest(pp) ||
        precomp_crs != pp.crs_digest)
      return false;
    return hv_crs_digest(pp) == pp.crs_digest &&
           hp_internal::validate_vme_precomputation_binding(pp.vme, vme);
  } catch (...) { return false; }
}

bool online_binding(const HvPublicParams& pp,
                    const HvVerifierPrecomputation& pre) {
  return internal::validate_hv_public_params_shape(pp) &&
         pre.crs_digest == pp.crs_digest &&
         pre.bp_parameter_digest == pp.bp_parameter_digest &&
         pre.bp_transcript_parameter_digest ==
             pp.bp_transcript_parameter_digest &&
         pre.P_digest == pp.P_digest &&
         pre.vme.pairing_x.size() == pp.vme.log_dimension + 1 &&
         pre.vme.delta1R.size() == pp.vme.log_dimension &&
         pre.vme.delta2R.size() == pp.vme.log_dimension;
}

std::vector<Scalar> batch_inverse(std::span<const Scalar> values) {
  if (values.empty()) return {};
  std::vector<Scalar> prefix(values.size()), out(values.size());
  Scalar product = fone();
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i].isZero()) throw std::invalid_argument("zero challenge");
    prefix[i] = product;
    product = fmul(product, values[i]);
  }
  Scalar inverse;
  Scalar::inv(inverse, product);
  for (std::size_t i = values.size(); i-- > 0;) {
    out[i] = fmul(inverse, prefix[i]);
    inverse = fmul(inverse, values[i]);
  }
  return out;
}

bool prepare(const HvPublicParams& pp, const HvStatement& statement,
             hv_internal::PreparedStatement& out,
             HvVerifyTimings* timings = nullptr,
             const Digest* known_context = nullptr) {
  out = {};
  const auto replay_start = Clock::now();
  std::vector<Scalar> proof_order;
  if (!internal::replay_bp_challenges_prevalidated(
          pp.bp, pp.bp_transcript_parameter_digest, statement.Z,
          statement.bulletproof, proof_order))
    return false;
  out.alphas.resize(pp.bp.d);
  for (std::size_t k = 1; k <= pp.bp.d; ++k)
    out.alphas[k - 1] = proof_order[pp.bp.d - k];
  if (timings) timings->bp_transcript_ms = milliseconds(replay_start);

  const auto inversion_start = Clock::now();
  out.alpha_inverses = batch_inverse(out.alphas);
  if (timings)
    timings->batch_inversion_ms = milliseconds(inversion_start);

  const auto weight_start = Clock::now();
  out.g = {fone()};
  out.h = {fone()};
  for (std::size_t k = 0; k < pp.bp.d; ++k) {
    const std::size_t old = out.g.size();
    std::vector<Scalar> next_g(2 * old), next_h(2 * old);
    for (std::size_t i = 0; i < old; ++i) {
      next_g[i] = out.g[i];
      next_g[old + i] = fmul(out.alphas[k], out.g[i]);
      next_h[i] = out.h[i];
      next_h[old + i] = fmul(out.alpha_inverses[k], out.h[i]);
    }
    out.g = std::move(next_g);
    out.h = std::move(next_h);
  }
  if (out.g.size() != pp.bp.n || out.h.size() != pp.bp.n) return false;
  if (timings)
    timings->weight_generation_ms = milliseconds(weight_start);

  const auto z_start = Clock::now();
  out.z_v.resize(pp.N);
  for (std::size_t i = 0; i < pp.bp.n; ++i) {
    out.z_v[i] =
        fmul(statement.bulletproof.x_final, out.g[i]);
    out.z_v[pp.bp.n + i] =
        fmul(statement.bulletproof.y_final, out.h[i]);
  }
  if (timings) timings->z_vector_ms = milliseconds(z_start);

  const auto x0_start = Clock::now();
  std::vector<Group> bases;
  std::vector<Scalar> scalars;
  bases.reserve(2 * pp.bp.d + 2);
  scalars.reserve(2 * pp.bp.d + 2);
  bases.push_back(statement.Z);
  scalars.push_back(fone());
  bases.push_back(pp.bp.K);
  scalars.push_back(fneg(fmul(statement.bulletproof.x_final,
                             statement.bulletproof.y_final)));
  for (std::size_t proof_index = 0; proof_index < pp.bp.d; ++proof_index) {
    const std::size_t k = pp.bp.d - proof_index;
    bases.push_back(statement.bulletproof.rounds[proof_index].A);
    scalars.push_back(out.alpha_inverses[k - 1]);
    bases.push_back(statement.bulletproof.rounds[proof_index].B);
    scalars.push_back(out.alphas[k - 1]);
  }
  out.X0 = msm(bases, scalars);
  if (timings) timings->x0_msm_ms = milliseconds(x0_start);

  const auto context_start = Clock::now();
  out.context = known_context ? *known_context : hv_context(pp, statement);
  if (timings && !known_context)
    timings->context_digest_ms = milliseconds(context_start);
  return true;
}

bool canonical_group(std::span<const std::uint8_t> bytes, Group& out) {
  return internal::canonical_group(bytes, out);
}

bool canonical_scalar(std::span<const std::uint8_t> bytes, Scalar& out) {
  return internal::canonical_scalar(bytes, out);
}

std::uint64_t read_u64(std::span<const std::uint8_t> bytes, std::size_t pos) {
  std::uint64_t out = 0;
  if (!internal::read_u64(bytes, pos, out))
    throw std::out_of_range("truncated u64");
  return out;
}

constexpr std::array<std::uint8_t, 8> kHvMagic{
    'B','P','H','V','G','2','0','1'};
constexpr std::uint16_t kHvVersion = 1;
constexpr std::size_t kHvHeader = 26;
constexpr std::array<std::uint8_t, 8> kHvStatementMagic{
    'B','P','H','V','S','T','0','1'};
constexpr std::size_t kHvStatementHeader = 42;
constexpr std::array<std::uint8_t, 8> kBpProofMagic{
    'B','P','I','P','A','G','2','1'};

bool verify_impl(const HvPublicParams& pp,
                 const HvVerifierPrecomputation& pre,
                 const HvStatement& statement, const HvProof& proof,
                 HvVerifyTimings* timings, bool proof_prevalidated,
                 bool reference, bool statement_prevalidated = false,
                 HvPreparedStatementCache* cache = nullptr) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  try {
    const auto validation_start = Clock::now();
    if (!online_binding(pp, pre) ||
        (!statement_prevalidated && !valid_statement(pp, statement)) ||
        (!proof_prevalidated &&
         !internal::validate_hv_proof_shape(pp, proof, true))) {
      if (timings) {
        timings->validation_ms = milliseconds(validation_start);
        timings->total_ms = milliseconds(total_start);
      }
      return false;
    }
    if (timings)
      timings->validation_ms = milliseconds(validation_start);
    hv_internal::PreparedStatement prepared;
    if (cache) {
      const auto context_start = Clock::now();
      const Digest current_context = hv_context(pp, statement);
      if (timings)
        timings->context_digest_ms = milliseconds(context_start);
      if (!internal::HvPreparedStatementCacheAccess::load(
              *cache, pp, current_context, prepared)) {
        if (!prepare(pp, statement, prepared, timings, &current_context))
          return false;
        internal::HvPreparedStatementCacheAccess::store(
            *cache, pp, prepared);
      }
    } else if (!prepare(pp, statement, prepared, timings)) {
      return false;
    }
    const auto vme_start = Clock::now();
    hp_internal::VmeVerificationStats vme_stats;
    const bool accepted = reference
        ? hp_internal::verify_vme_reference(
              pp.vme, pre.vme, prepared.context, prepared.X0,
              prepared.z_v, proof.vme_proof, &vme_stats)
        : hp_internal::verify_vme_prevalidated_optimized(
              pp.vme, pre.vme, prepared.context, prepared.X0,
              prepared.z_v, proof.vme_proof, &vme_stats);
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
    return accepted;
  } catch (...) {
    if (timings) timings->total_ms = milliseconds(total_start);
    return false;
  }
}

}

HvSetupResult setup_helper_verifier(
    std::size_t n, std::span<const std::uint8_t> seed,
    HvSetupTimings* timings) {
  initialize();
  if (timings) *timings = {};
  const auto setup_start = Clock::now();
  if (!power_of_two(n) ||
      n > std::numeric_limits<std::size_t>::max() / 2)
    throw std::invalid_argument("n must be a power of two");
  HvSetupResult out;
  out.pp.bp = Setup(n, seed);
  out.pp.N = 2 * n;
  std::vector<Group> P;
  P.reserve(out.pp.N);
  P.insert(P.end(), out.pp.bp.G.begin(), out.pp.bp.G.end());
  P.insert(P.end(), out.pp.bp.H.begin(), out.pp.bp.H.end());
  Bytes vme_seed;
  frame(vme_seed, "BPVME/HV/VME-SEED/v1");
  frame(vme_seed, seed);
  frame(vme_seed, bp_parameter_digest(out.pp.bp));
  out.pp.vme = hp_internal::setup_vme(P, vme_seed, kHvVmeDomain);
  out.pp.bp_parameter_digest = bp_parameter_digest(out.pp.bp);
  out.pp.bp_transcript_parameter_digest =
      internal::bp_public_parameter_digest(out.pp.bp);
  out.pp.P_digest = p_digest(out.pp);
  out.pp.crs_digest = hv_crs_digest(out.pp);
  const auto precompute_start = Clock::now();
  const VmePrecomputation vme = hp_internal::precompute_vme(out.pp.vme);
  out.helper_precomp = {
      vme, out.pp.crs_digest, out.pp.bp_transcript_parameter_digest};
  out.verifier_precomp = {
      vme, out.pp.crs_digest, out.pp.bp_parameter_digest,
      out.pp.bp_transcript_parameter_digest, out.pp.P_digest};
  if (timings) {
    timings->verifier_precomputation_ms = milliseconds(precompute_start);
    timings->setup_ms =
        milliseconds(setup_start) - timings->verifier_precomputation_ms;
  }
  return out;
}

std::optional<HvProof> prove_helper_verifier(
    const HvPublicParams& pp, const HvHelperPrecomputation& pre,
    const HvStatement& statement, HvProveTimings* timings) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  try {
    if (!fast_binding(pp, pre.vme, pre.crs_digest) ||
        pre.bp_transcript_parameter_digest !=
            pp.bp_transcript_parameter_digest ||
        !valid_statement(pp, statement))
      return std::nullopt;
    const auto prepare_start = Clock::now();
    hv_internal::PreparedStatement prepared;
    if (!prepare(pp, statement, prepared)) return std::nullopt;
    if (timings)
      timings->statement_preparation_ms = milliseconds(prepare_start);
    const auto relation_start = Clock::now();
    if (msm(pp.vme.fixed_P, prepared.z_v) != prepared.X0)
      return std::nullopt;
    if (timings)
      timings->relation_msm_ms = milliseconds(relation_start);
    const auto vme_start = Clock::now();
    HvProof proof{hp_internal::prove_vme(
        pp.vme, pre.vme, prepared.context, prepared.X0, prepared.z_v)};
    if (timings) {
      timings->vme_prove_ms = milliseconds(vme_start);
      timings->total_ms = milliseconds(total_start);
    }
    return proof;
  } catch (...) { return std::nullopt; }
}

bool verify_helper_verifier(
    const HvPublicParams& pp, const HvVerifierPrecomputation& pre,
    const HvStatement& statement, const HvProof& proof,
    HvVerifyTimings* timings) {
  return verify_impl(pp, pre, statement, proof, timings, false, false);
}

bool verify_helper_verifier_cached(
    const HvPublicParams& pp, const HvVerifierPrecomputation& pre,
    const HvStatement& statement, const HvProof& proof,
    HvPreparedStatementCache& cache, HvVerifyTimings* timings) {
  return verify_impl(
      pp, pre, statement, proof, timings, false, false, false, &cache);
}

bool verify_helper_verifier_serialized(
    const HvPublicParams& pp, const HvVerifierPrecomputation& pre,
    const HvStatement& statement, std::span<const std::uint8_t> bytes,
    HvVerifyTimings* timings) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  const auto validation_start = Clock::now();
  HvProof proof;
  if (!deserialize_hv_proof(pp, bytes, proof)) {
    if (timings) {
      timings->validation_ms = milliseconds(validation_start);
      timings->total_ms = milliseconds(total_start);
    }
    return false;
  }
  const double parse_ms = milliseconds(validation_start);
  HvVerifyTimings inner;
  const bool accepted =
      verify_impl(pp, pre, statement, proof, &inner, true, false);
  if (timings) {
    *timings = inner;
    timings->validation_ms += parse_ms;
    timings->total_ms = milliseconds(total_start);
  }
  return accepted;
}

bool verify_helper_verifier_serialized(
    const HvPublicParams& pp, const HvVerifierPrecomputation& pre,
    std::span<const std::uint8_t> statement_bytes,
    std::span<const std::uint8_t> proof_bytes, HvVerifyTimings* timings) {
  if (timings) *timings = {};
  const auto total_start = Clock::now();
  const auto validation_start = Clock::now();
  HvStatement statement;
  HvProof proof;
  if (!deserialize_hv_statement(pp, statement_bytes, statement) ||
      !deserialize_hv_proof(pp, proof_bytes, proof)) {
    if (timings) {
      timings->validation_ms = milliseconds(validation_start);
      timings->total_ms = milliseconds(total_start);
    }
    return false;
  }
  const double parse_ms = milliseconds(validation_start);
  HvVerifyTimings inner;
  const bool accepted = verify_impl(
      pp, pre, statement, proof, &inner, true, false, true);
  if (timings) {
    *timings = inner;
    timings->validation_ms += parse_ms;
    timings->total_ms = milliseconds(total_start);
  }
  return accepted;
}

HvInstance generate_helper_verifier_instance(const HvPublicParams& pp) {
  if (!validate_public_params(pp.bp))
    throw std::invalid_argument("invalid helper-verifier parameters");
  HvInstance out;
  out.x.resize(pp.bp.n);
  out.y.resize(pp.bp.n);
  for (std::size_t i = 0; i < pp.bp.n; ++i) {
    out.x[i].setByCSPRNG();
    out.y[i].setByCSPRNG();
  }
  out.statement.Z = commit(pp.bp, out.x, out.y);
  out.statement.bulletproof =
      Prove(pp.bp, out.statement.Z, out.x, out.y);
  return out;
}

std::size_t hv_proof_payload_bytes(std::size_t n) {
  const std::size_t D = exact_log2(n) + 1;
  std::size_t gt_count = 0, gt_total = 0, g1_total = 0, total = 0;
  if (!internal::checked_mul(11, D, gt_count) || gt_count < 4 ||
      !internal::checked_mul(gt_count - 4, hp_internal::gt_bytes(), gt_total) ||
      !internal::checked_mul(2, hp_internal::g1_bytes(), g1_total) ||
      !internal::checked_add(gt_total, g1_total, total) ||
      !internal::checked_add(total, group_bytes(), total))
    throw std::overflow_error("HV proof size overflow");
  return total;
}

std::size_t hv_proof_wire_bytes(std::size_t n) {
  std::size_t total = 0;
  if (!internal::checked_add(kHvHeader, hv_proof_payload_bytes(n), total))
    throw std::overflow_error("HV proof wire size overflow");
  return total;
}

Bytes serialize_hv_proof(const HvPublicParams& pp, const HvProof& proof) {
  if (!internal::validate_hv_public_params_shape(pp) ||
      !internal::validate_hv_proof_shape(pp, proof, false))
    throw std::invalid_argument("invalid helper-verifier proof");
  Bytes out;
  out.reserve(hv_proof_wire_bytes(pp.bp.n));
  raw(out, kHvMagic);
  out.push_back(static_cast<std::uint8_t>(kHvVersion >> 8));
  out.push_back(static_cast<std::uint8_t>(kHvVersion));
  raw(out, u64be(pp.N));
  raw(out, u64be(pp.vme.log_dimension));
  for (const auto& claims : proof.vme_proof.rexp_claims)
    for (const GT* field : {&claims.E, &claims.F, &claims.TL, &claims.TR})
      raw(out, hp_internal::serialize_gt(*field));
  raw(out, hp_internal::serialize_g1(proof.vme_proof.R));
  for (const auto& fold : proof.vme_proof.dory_folds)
    for (const GT* field :
         {&fold.D1L, &fold.D1R, &fold.D2L, &fold.D2R,
          &fold.W1, &fold.W2})
      raw(out, hp_internal::serialize_gt(*field));
  for (const auto& U : proof.vme_proof.batch_U)
    raw(out, hp_internal::serialize_gt(U));
  raw(out, hp_internal::serialize_g1(proof.vme_proof.phi_final));
  raw(out, serialize(proof.vme_proof.theta_final));
  if (out.size() != hv_proof_wire_bytes(pp.bp.n))
    throw std::runtime_error("unexpected HV proof encoding size");
  return out;
}

bool deserialize_hv_proof(const HvPublicParams& pp,
                          std::span<const std::uint8_t> bytes,
                          HvProof& proof) noexcept {
  proof = {};
  try {
    if (!internal::validate_hv_public_params_shape(pp) ||
        bytes.size() != hv_proof_wire_bytes(pp.bp.n) ||
        !std::equal(kHvMagic.begin(), kHvMagic.end(), bytes.begin()) ||
        bytes[8] != static_cast<std::uint8_t>(kHvVersion >> 8) ||
        bytes[9] != static_cast<std::uint8_t>(kHvVersion) ||
        read_u64(bytes, 10) != pp.N ||
        read_u64(bytes, 18) != pp.vme.log_dimension)
      return false;
    const std::size_t D = pp.vme.log_dimension;
    const std::size_t gtb = hp_internal::gt_bytes();
    std::size_t pos = kHvHeader;
    HvProof candidate;
    candidate.vme_proof.rexp_claims.resize(D - 1);
    for (auto& claims : candidate.vme_proof.rexp_claims)
      for (GT* field : {&claims.E, &claims.F, &claims.TL, &claims.TR}) {
        if (!hp_internal::deserialize_gt_canonical_unchecked(
                bytes.subspan(pos, gtb), *field))
          return false;
        pos += gtb;
      }
    if (!hp_internal::deserialize_g1(
            bytes.subspan(pos, hp_internal::g1_bytes()),
            candidate.vme_proof.R))
      return false;
    pos += hp_internal::g1_bytes();
    candidate.vme_proof.dory_folds.resize(D);
    for (auto& fold : candidate.vme_proof.dory_folds)
      for (GT* field :
           {&fold.D1L, &fold.D1R, &fold.D2L, &fold.D2R,
            &fold.W1, &fold.W2}) {
        if (!hp_internal::deserialize_gt_canonical_unchecked(
                bytes.subspan(pos, gtb), *field))
          return false;
        pos += gtb;
      }
    candidate.vme_proof.batch_U.resize(D);
    for (auto& U : candidate.vme_proof.batch_U) {
      if (!hp_internal::deserialize_gt_canonical_unchecked(
              bytes.subspan(pos, gtb), U))
        return false;
      pos += gtb;
    }
    if (!hp_internal::deserialize_g1(
            bytes.subspan(pos, hp_internal::g1_bytes()),
            candidate.vme_proof.phi_final))
      return false;
    pos += hp_internal::g1_bytes();
    if (!canonical_group(bytes.subspan(pos, group_bytes()),
                         candidate.vme_proof.theta_final))
      return false;
    pos += group_bytes();
    if (pos != bytes.size()) return false;
    if (!hp_internal::validate_vme_proof_batched(
            pp.vme, candidate.vme_proof))
      return false;
    proof = std::move(candidate);
    return true;
  } catch (...) { return false; }
}

std::size_t hv_statement_wire_bytes(std::size_t n) {
  std::size_t total = 0;
  if (!internal::checked_add(kHvStatementHeader, group_bytes(), total) ||
      !internal::checked_add(total, proof_wire_bytes(n), total))
    throw std::overflow_error("HV statement wire size overflow");
  return total;
}

Bytes serialize_hv_statement(const HvPublicParams& pp,
                             const HvStatement& statement) {
  if (!internal::validate_hv_public_params_shape(pp) ||
      !valid_statement(pp, statement))
    throw std::invalid_argument("invalid helper-verifier statement");
  Bytes out;
  out.reserve(hv_statement_wire_bytes(pp.bp.n));
  raw(out, kHvStatementMagic);
  out.push_back(static_cast<std::uint8_t>(kHvVersion >> 8));
  out.push_back(static_cast<std::uint8_t>(kHvVersion));
  raw(out, pp.crs_digest);
  raw(out, serialize(statement.Z));
  raw(out, serialize_proof(pp.bp, statement.bulletproof));
  if (out.size() != hv_statement_wire_bytes(pp.bp.n))
    throw std::runtime_error("unexpected HV statement encoding size");
  return out;
}

bool deserialize_hv_statement(const HvPublicParams& pp,
                              std::span<const std::uint8_t> bytes,
                              HvStatement& statement) noexcept {
  statement = {};
  try {
    if (!internal::validate_hv_public_params_shape(pp) ||
        bytes.size() != hv_statement_wire_bytes(pp.bp.n) ||
        !std::equal(kHvStatementMagic.begin(), kHvStatementMagic.end(),
                    bytes.begin()) ||
        bytes[8] != static_cast<std::uint8_t>(kHvVersion >> 8) ||
        bytes[9] != static_cast<std::uint8_t>(kHvVersion) ||
        !std::equal(pp.crs_digest.begin(), pp.crs_digest.end(),
                    bytes.begin() + 10))
      return false;
    HvStatement candidate;
    std::size_t pos = kHvStatementHeader;
    if (!canonical_group(bytes.subspan(pos, group_bytes()), candidate.Z))
      return false;
    pos += group_bytes();
    const auto bp_wire = bytes.subspan(pos);
    if (bp_wire.size() != proof_wire_bytes(pp.bp.n) ||
        !std::equal(kBpProofMagic.begin(), kBpProofMagic.end(),
                    bp_wire.begin()) ||
        bp_wire[8] != 0 || bp_wire[9] != 1 ||
        read_u64(bp_wire, 10) != pp.bp.n ||
        read_u64(bp_wire, 18) != pp.bp.d)
      return false;
    std::size_t bp_pos = 26;
    candidate.bulletproof.rounds.resize(pp.bp.d);
    for (auto& round : candidate.bulletproof.rounds) {
      if (!canonical_group(
              bp_wire.subspan(bp_pos, group_bytes()), round.A))
        return false;
      bp_pos += group_bytes();
      if (!canonical_group(
              bp_wire.subspan(bp_pos, group_bytes()), round.B))
        return false;
      bp_pos += group_bytes();
    }
    if (!canonical_scalar(
            bp_wire.subspan(bp_pos, scalar_bytes()),
            candidate.bulletproof.x_final))
      return false;
    bp_pos += scalar_bytes();
    if (!canonical_scalar(
            bp_wire.subspan(bp_pos, scalar_bytes()),
            candidate.bulletproof.y_final))
      return false;
    bp_pos += scalar_bytes();
    if (bp_pos != bp_wire.size()) return false;
    statement = std::move(candidate);
    return true;
  } catch (...) { return false; }
}

namespace hv_internal {

bool prepare_statement_for_test(const HvPublicParams& pp,
                                const HvStatement& statement,
                                PreparedStatement& out) {
  try {
    return valid_statement(pp, statement) && prepare(pp, statement, out);
  } catch (...) { return false; }
}

bool verify_helper_verifier_reference_for_test(
    const HvPublicParams& pp, const HvVerifierPrecomputation& pre,
    const HvStatement& statement, const HvProof& proof) {
  return verify_impl(pp, pre, statement, proof, nullptr, false, true);
}

Group direct_vme_relation_for_test(const HvPublicParams& pp,
                                   const PreparedStatement& prepared) {
  return msm(pp.vme.fixed_P, prepared.z_v);
}

}
}

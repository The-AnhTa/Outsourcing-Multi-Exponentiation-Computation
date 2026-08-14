#include "internal/protocol_validation.hpp"

#include "hp_vme_internal.hpp"
#include "internal/protocol_utils.hpp"

#include <limits>

namespace bp::internal {

bool validate_bp_public_params(const PublicParams& pp) noexcept {
  try {
    initialize();
    const bool standard_transcript =
        pp.transcript_domain == kTranscriptDomain && !pp.transcript_crs_digest;
    const bool hp_transcript =
        pp.transcript_domain == kHpTranscriptDomain && pp.transcript_crs_digest;
    if (!power_of_two(pp.n) || pp.d != exact_log2(pp.n) ||
        pp.G.size() != pp.n || pp.H.size() != pp.n ||
        pp.group_identifier != kGroupIdentifier ||
        pp.scalar_modulus != kScalarModulus ||
        pp.hash_suite_identifier != kHashSuiteIdentifier ||
        (!standard_transcript && !hp_transcript))
      return false;
    for (const auto& point : pp.G)
      if (!valid_group(point, true)) return false;
    for (const auto& point : pp.H)
      if (!valid_group(point, true)) return false;
    return valid_group(pp.K, true);
  } catch (...) {
    return false;
  }
}

bool validate_bp_proof_shape(const PublicParams& pp, const Group& statement,
                             const Proof& proof) noexcept {
  if (!valid_group(statement) || proof.rounds.size() != pp.d) return false;
  for (const auto& round : proof.rounds)
    if (!valid_group(round.A) || !valid_group(round.B)) return false;
  return true;
}

bool validate_hp_public_params_shape(const HpPublicParams& pp) noexcept {
  try {
    std::size_t expected_dimension = 0;
    if (!validate_bp_public_params(pp.bp) ||
        !checked_mul(pp.bp.n, 2, expected_dimension) ||
        pp.N != expected_dimension || pp.vme.dimension != pp.N ||
        pp.vme.log_dimension != pp.bp.d + 1 ||
        pp.bp.transcript_domain != kHpTranscriptDomain ||
        !pp.bp.transcript_crs_digest ||
        *pp.bp.transcript_crs_digest != pp.crs_digest ||
        pp.vme.fixed_P.size() != pp.N ||
        !hp_internal::validate_vme_params(pp.vme))
      return false;
    for (std::size_t i = 0; i < pp.bp.n; ++i)
      if (pp.vme.fixed_P[i] != pp.bp.G[i] ||
          pp.vme.fixed_P[pp.bp.n + i] != pp.bp.H[i])
        return false;
    return true;
  } catch (...) {
    return false;
  }
}

bool validate_hp_proof_shape(const HpPublicParams& pp,
                             const HpProof& proof) noexcept {
  if (proof.rounds.size() != pp.bp.d ||
      (pp.bp.d == 0 ? proof.vme_proof.has_value()
                    : !proof.vme_proof.has_value()))
    return false;
  for (const auto& round : proof.rounds)
    if (!valid_group(round.A) || !valid_group(round.B)) return false;
  return !proof.vme_proof ||
         hp_internal::validate_vme_proof(pp.vme, *proof.vme_proof);
}

bool validate_hv_public_params_shape(const HvPublicParams& pp) noexcept {
  try {
    std::size_t expected_dimension = 0;
    if (!validate_bp_public_params(pp.bp) ||
        !checked_mul(pp.bp.n, 2, expected_dimension) ||
        pp.N != expected_dimension || pp.vme.dimension != pp.N ||
        pp.vme.log_dimension != pp.bp.d + 1 ||
        pp.vme.transcript_domain != kHvVmeDomain ||
        pp.vme.fixed_P.size() != pp.N ||
        !hp_internal::validate_vme_params(pp.vme))
      return false;
    for (std::size_t i = 0; i < pp.bp.n; ++i)
      if (pp.vme.fixed_P[i] != pp.bp.G[i] ||
          pp.vme.fixed_P[pp.bp.n + i] != pp.bp.H[i])
        return false;
    return true;
  } catch (...) {
    return false;
  }
}

bool validate_hv_statement_shape(const HvPublicParams& pp,
                                 const HvStatement& statement) noexcept {
  return validate_bp_proof_shape(pp.bp, statement.Z, statement.bulletproof);
}

bool validate_hv_proof_shape(const HvPublicParams& pp,
                             const HvProof& proof,
                             bool batched_gt_validation) noexcept {
  return batched_gt_validation
             ? hp_internal::validate_vme_proof_batched(pp.vme, proof.vme_proof)
             : hp_internal::validate_vme_proof(pp.vme, proof.vme_proof);
}

}  // namespace bp::internal

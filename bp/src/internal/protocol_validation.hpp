#pragma once

#include "bp/helper_verifier.hpp"

namespace bp::internal {

bool validate_bp_public_params(const PublicParams& pp) noexcept;
bool validate_bp_proof_shape(const PublicParams& pp, const Group& statement,
                             const Proof& proof) noexcept;
bool validate_hp_public_params_shape(const HpPublicParams& pp) noexcept;
bool validate_hp_proof_shape(const HpPublicParams& pp,
                             const HpProof& proof) noexcept;
bool validate_hv_public_params_shape(const HvPublicParams& pp) noexcept;
bool validate_hv_statement_shape(const HvPublicParams& pp,
                                 const HvStatement& statement) noexcept;
bool validate_hv_proof_shape(const HvPublicParams& pp,
                             const HvProof& proof,
                             bool batched_gt_validation) noexcept;

}  // namespace bp::internal

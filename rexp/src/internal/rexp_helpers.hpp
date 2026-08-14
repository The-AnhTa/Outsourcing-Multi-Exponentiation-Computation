#pragma once

#include "rexp/rexp.hpp"

namespace rexp::internal {

bool validate_rexp_proof_gt(
    const RexpProof& proof, RexpProofValidationMetrics* metrics);

inline void require_bound(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement) {
    if (statement.crsDigest() != params.digest()
        || statement.H().size() != params.n()) {
        throw std::invalid_argument(
            "prepared statement belongs to another CRS");
    }
}

inline bool rexp_proof_shape(
    const PreparedPublicParameters& params, const RexpProof& proof) {
    if (proof.doryProofs.size() != params.d()
        || proof.dynamicRoundMessages.size()
            != (params.d() ? params.d() - 1 : 0)
        || !proof.R.isValid() || !proof.R.isValidOrder()) {
        return false;
    }
    for (std::size_t round = 0; round < params.d(); ++round) {
        if (proof.doryProofs[round].rounds.size()
            != params.d() - round - 1) {
            return false;
        }
    }
    return true;
}

} 

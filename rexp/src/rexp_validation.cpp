#include "rexp/verify.hpp"

#include "internal/rexp_helpers.hpp"

#include <chrono>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rexp {
namespace internal {

void validate_rexp_proof_shape(const RexpProof& proof, std::size_t d) {
    if (d >= std::numeric_limits<std::size_t>::digits
        || proof.doryProofs.size() != d
        || proof.dynamicRoundMessages.size() != (d ? d - 1 : 0)
        || !proof.R.isValid() || !proof.R.isValidOrder()) {
        throw std::invalid_argument("invalid proof shape");
    }
    for (std::size_t round = 0; round < d; ++round) {
        const DoryProof& dory = proof.doryProofs[round];
        if (dory.rounds.size() != d - round - 1
            || !dory.PhiFinal.isValid() || !dory.PhiFinal.isValidOrder()
            || !dory.ThetaFinal.isValid() || !dory.ThetaFinal.isValidOrder()) {
            throw std::invalid_argument("invalid embedded Dory proof");
        }
    }
}

bool validate_rexp_proof_gt(
    const RexpProof& proof, RexpProofValidationMetrics* metrics) {
    using Clock = std::chrono::steady_clock;
    const auto start = metrics ? Clock::now() : Clock::time_point{};
    std::size_t checked = 0;
    const auto check = [&](const GT& value) {
        ++checked;
        if (!mcl::bn::isValidGT(value)) {
            throw std::invalid_argument("proof GT is outside target subgroup");
        }
    };
    for (const auto& message : proof.dynamicRoundMessages) {
        check(message.E);
        check(message.F);
        check(message.TL);
        check(message.TR);
    }
    for (const auto& dory : proof.doryProofs) {
        for (const auto& round : dory.rounds) {
            check(round.D1L);
            check(round.D1R);
            check(round.D2L);
            check(round.D2R);
            check(round.W1);
            check(round.W2);
        }
    }
    if (metrics) {
        metrics->gt_elements_checked = checked;
        metrics->gt_subgroup_validation_ms =
            std::chrono::duration<double, std::milli>(Clock::now() - start)
                .count();
    }
    return true;
}

} 

bool IsValidGTSubgroup(const GT& value) {
    return mcl::bn::isValidGT(value);
}

ValidatedRexpProof ValidateRexpProof(
    const RexpProof& proof,
    std::size_t d,
    RexpProofValidationMetrics* metrics) {
    internal::validate_rexp_proof_shape(proof, d);
    internal::validate_rexp_proof_gt(proof, metrics);
    ValidatedRexpProof result;
    result.proof_ = proof;
    result.d_ = d;
    return result;
}

ValidatedRexpProof ValidateRexpProof(
    RexpProof&& proof,
    std::size_t d,
    RexpProofValidationMetrics* metrics) {
    internal::validate_rexp_proof_shape(proof, d);
    internal::validate_rexp_proof_gt(proof, metrics);
    ValidatedRexpProof result;
    result.proof_ = std::move(proof);
    result.d_ = d;
    return result;
}

} 

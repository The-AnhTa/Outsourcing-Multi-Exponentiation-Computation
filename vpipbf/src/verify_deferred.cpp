#include "vpip_bf/verify_deferred.hpp"
#include "vpip_bf/verify_combined.hpp"
#include "vpip_bf/verify_online.hpp"
#include "vpip_bf/group_utils.hpp"
#include "internal/verification.hpp"

namespace vpip_bf {

std::optional<ValidatedVerificationInputs> PrevalidateVerificationInputs(
    const VpipBfCRS& crs, const VpipBfPrecomputation& precomputation,
    const VpipBfStatement& statement, const VpipBfProof& proof) {
    if (!validate_verification_inputs(crs, precomputation, statement, proof))
        return std::nullopt;
    return ValidatedVerificationInputs(crs, precomputation, statement, proof);
}

bool verify_online(const ValidatedVerificationInputs& inputs) {
    return verify_core_symbolic_unchecked(
        inputs.crs(), inputs.precomputation(),
        inputs.statement(), inputs.proof()).accepted;
}

bool verify_online_with_timing(const ValidatedVerificationInputs& inputs,
                               OnlineTimingBreakdown& timing) {
    timing = {};
    return verify_core_symbolic_unchecked(
        inputs.crs(), inputs.precomputation(), inputs.statement(),
        inputs.proof(), &timing).accepted;
}

DeferredVerificationTrace verify_deferred_with_trace(
    const VpipBfCRS& crs, const VpipBfPrecomputation& precomputation,
    const VpipBfStatement& statement, const VpipBfProof& proof) {
    DeferredVerificationTrace trace;
    const auto inputs = PrevalidateVerificationInputs(
        crs, precomputation, statement, proof);
    if (!inputs) return trace;
    const auto reference = verify_core_unchecked(
        inputs->crs(), inputs->precomputation(),
        inputs->statement(), inputs->proof());
    const auto optimized = verify_core_symbolic_unchecked(
        inputs->crs(), inputs->precomputation(),
        inputs->statement(), inputs->proof());
    const bool engines_agree = reference.rho == optimized.rho
        && reference.beta == optimized.beta
        && reference.alpha == optimized.alpha
        && reference.gamma == optimized.gamma
        && reference.epsilon == optimized.epsilon
        && reference.dory_residual == optimized.dory_residual
        && reference.rexp_residual == optimized.rexp_residual;
    trace.accepted = reference.accepted && optimized.accepted
        && engines_agree;
    trace.dory_accepted = reference.dory_accepted
        && optimized.dory_accepted;
    trace.rexp_accepted = reference.rexp_accepted
        && optimized.rexp_accepted;
    trace.rho = optimized.rho;
    trace.beta = optimized.beta;
    trace.alpha = optimized.alpha;
    trace.gamma = optimized.gamma;
    trace.epsilon = optimized.epsilon;
    trace.dory_residual_reference = reference.dory_residual;
    trace.dory_residual_pippenger = optimized.dory_residual;
    trace.rexp_residual_reference = reference.rexp_residual;
    trace.rexp_residual_pippenger = optimized.rexp_residual;
    return trace;
}

bool verify_deferred(const VpipBfCRS& crs,
                     const VpipBfPrecomputation& precomputation,
                     const VpipBfStatement& statement,
                     const VpipBfProof& proof) {
    return verify_deferred_with_trace(
        crs, precomputation, statement, proof).accepted;
}

CombinedVerificationTrace verify_deferred_combined_with_trace(
    const VpipBfCRS& crs, const VpipBfPrecomputation& precomputation,
    const VpipBfStatement& statement, const VpipBfProof& proof) {
    CombinedVerificationTrace trace;
    const auto inputs = PrevalidateVerificationInputs(
        crs, precomputation, statement, proof);
    if (!inputs) return trace;
    const auto optimized = verify_core_symbolic_unchecked(
        inputs->crs(), inputs->precomputation(),
        inputs->statement(), inputs->proof());
    const auto reference = verify_core_unchecked(
        inputs->crs(), inputs->precomputation(),
        inputs->statement(), inputs->proof());
    const bool engines_agree = optimized.rho == reference.rho
        && optimized.beta == reference.beta
        && optimized.alpha == reference.alpha
        && optimized.gamma == reference.gamma
        && optimized.epsilon == reference.epsilon
        && optimized.dory_residual == reference.dory_residual
        && optimized.rexp_residual == reference.rexp_residual;
    trace.accepted = optimized.accepted && reference.accepted
        && engines_agree;
    trace.rho = optimized.rho;
    trace.beta = optimized.beta;
    trace.alpha = optimized.alpha;
    trace.gamma = optimized.gamma;
    trace.epsilon = optimized.epsilon;
    trace.evaluated_dory_residual = optimized.dory_residual;
    trace.evaluated_rexp_residual = optimized.rexp_residual;
    trace.evaluated_combined_residual = gt_mul(
        optimized.dory_residual, optimized.rexp_residual);
    trace.evaluated_combined_reference = gt_mul(
        reference.dory_residual, reference.rexp_residual);
    return trace;
}

bool verify_deferred_combined(
    const VpipBfCRS& crs, const VpipBfPrecomputation& precomputation,
    const VpipBfStatement& statement, const VpipBfProof& proof) {
    return verify_deferred_combined_with_trace(
        crs, precomputation, statement, proof).accepted;
}

CombinedVerificationTrace verify_online_with_trace(
    const ValidatedVerificationInputs& inputs) {
    CombinedVerificationTrace trace;
    const auto result = verify_core_symbolic_unchecked(
        inputs.crs(), inputs.precomputation(),
        inputs.statement(), inputs.proof());
    trace.accepted = result.accepted;
    trace.rho = result.rho;
    trace.beta = result.beta;
    trace.alpha = result.alpha;
    trace.gamma = result.gamma;
    trace.epsilon = result.epsilon;
    trace.evaluated_dory_residual = result.dory_residual;
    trace.evaluated_rexp_residual = result.rexp_residual;
    trace.evaluated_combined_residual = gt_mul(
        result.dory_residual, result.rexp_residual);
    return trace;
}

} 

#include "vme_ibf/verify_deferred.hpp"

#include "vme_ibf/verify_combined.hpp"
#include "vme_ibf/verify_online.hpp"
#include "internal/crypto.hpp"
#include "internal/core_context.hpp"
#include "internal/protocol.hpp"
#include "internal/verification.hpp"

namespace vme_ibf {
namespace {

struct Resolvers {
    GtAtomResolver gt;
    PairingAtomResolver pairing;
};

Resolvers resolvers(const VmeIbfCRS& crs,
                    const VmeIbfPrecomputation& precomputation,
                    const VmeIbfStatement& statement,
                    const VmeIbfProof& proof,
                    const internal::VerificationEquations& equations) {
    const auto* crs_ptr = &crs;
    const auto* precomputation_ptr = &precomputation;
    const auto* statement_ptr = &statement;
    const auto* proof_ptr = &proof;
    const auto* equations_ptr = &equations;
    return {
        [precomputation_ptr, proof_ptr](GtAtomId id) -> const GT& {
            return internal::resolve_gt_atom(
                id, *precomputation_ptr, *proof_ptr);
        },
        [crs_ptr, statement_ptr, proof_ptr, equations_ptr](PairingAtomId id) {
            return internal::resolve_pairing_atom(
                id, *crs_ptr, *statement_ptr, *proof_ptr, *equations_ptr);
        }};
}

DeferredVerificationTrace evaluate_deferred(
    const VmeIbfCRS& crs,
    const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement,
    const VmeIbfProof& proof,
    bool audit_reference) {
    DeferredVerificationTrace trace;
    internal::VerificationEquations equations;
    if (!internal::build_verification_equations(
            crs, precomputation, statement, proof, equations)) return trace;
    trace.rho = equations.challenges.rho;
    trace.beta = equations.challenges.beta;
    trace.alpha = equations.challenges.alpha;
    trace.gamma = equations.challenges.gamma;
    trace.epsilon = equations.challenges.epsilon;
    trace.dory_gt_terms_before_normalize = equations.dory.gt_terms.size();
    trace.rexp_gt_terms_before_normalize = equations.rexp.gt_terms.size();
    normalize(equations.dory);
    normalize(equations.rexp);
    trace.dory_gt_terms_after_normalize = equations.dory.gt_terms.size();
    trace.rexp_gt_terms_after_normalize = equations.rexp.gt_terms.size();
    trace.dory_pairing_terms = equations.dory.pairing_terms.size();
    trace.rexp_pairing_terms = equations.rexp.pairing_terms.size();
    const auto resolver = resolvers(
        crs, precomputation, statement, proof, equations);
    if (audit_reference) {
        trace.dory_residual_reference = evaluate_symbolic_expression(
            equations.dory, resolver.gt, resolver.pairing,
            MultiexpMode::Reference);
        trace.rexp_residual_reference = evaluate_symbolic_expression(
            equations.rexp, resolver.gt, resolver.pairing,
            MultiexpMode::Reference);
    }
    trace.dory_residual_pippenger = evaluate_symbolic_expression(
        equations.dory, resolver.gt, resolver.pairing, MultiexpMode::Pippenger);
    trace.rexp_residual_pippenger = evaluate_symbolic_expression(
        equations.rexp, resolver.gt, resolver.pairing, MultiexpMode::Pippenger);
    GT one;
    one.setOne();
    trace.dory_accepted = trace.dory_residual_pippenger == one;
    trace.rexp_accepted = trace.rexp_residual_pippenger == one;
    trace.accepted = trace.dory_accepted && trace.rexp_accepted;
    trace.gt_multiexp_calls = 2;
    trace.pairing_product_calls = 2;
    return trace;
}

CombinedVerificationTrace evaluate_combined(
    const VmeIbfCRS& crs,
    const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement,
    const VmeIbfProof& proof,
    bool audit_reference,
    const Digest* start_state = nullptr) {
    CombinedVerificationTrace trace;
    const auto total_start = internal::Clock::now();
    internal::VerificationEquations equations;
    if (!internal::build_verification_equations(
            crs, precomputation, statement, proof, equations,
            start_state)) return trace;
    const auto resolver = resolvers(
        crs, precomputation, statement, proof, equations);
    Transcript transcript = Transcript::resume(equations.challenges.after_epsilon);
    const Fr eta = internal::derive_eta(transcript, crs.d);
    auto combined = equations.dory;
    multiply_in_place(combined, powered(equations.rexp, eta));
    trace.gt_terms_before_normalize = combined.gt_terms.size();
    trace.pairing_terms_before_normalize = combined.pairing_terms.size();
    const auto before_normalize = internal::Clock::now();
    normalize(combined);
    const auto after_normalize = internal::Clock::now();
    trace.gt_terms_after_normalize = combined.gt_terms.size();
    trace.pairing_terms_after_normalize = combined.pairing_terms.size();
    trace.rho = equations.challenges.rho;
    trace.beta = equations.challenges.beta;
    trace.alpha = equations.challenges.alpha;
    trace.gamma = equations.challenges.gamma;
    trace.epsilon = equations.challenges.epsilon;
    trace.eta = eta;
    trace.final_transcript_digest = transcript.digest();
    if (audit_reference) {
        trace.evaluated_dory_residual = evaluate_symbolic_expression(
            equations.dory, resolver.gt, resolver.pairing,
            MultiexpMode::Reference);
        trace.evaluated_rexp_residual = evaluate_symbolic_expression(
            equations.rexp, resolver.gt, resolver.pairing,
            MultiexpMode::Reference);
        trace.evaluated_combined_reference = evaluate_symbolic_expression(
            combined, resolver.gt, resolver.pairing, MultiexpMode::Reference);
    }
    PairingProductStats stats;
    trace.evaluated_combined_residual = evaluate_symbolic_expression(
        combined, resolver.gt, resolver.pairing, MultiexpMode::Pippenger,
        nullptr, &stats);
    const auto after_evaluation = internal::Clock::now();
    GT one;
    one.setOne();
    trace.accepted = trace.evaluated_combined_residual == one;
    trace.coalesced_pairing_terms = stats.coalesced_pairing_terms;
    trace.miller_loop_batches = stats.miller_loop_batches;
    trace.miller_loop_terms = stats.miller_loop_terms;
    trace.final_exponentiations = stats.final_exponentiations;
    trace.gt_msm_ms = stats.gt_msm_ms;
    trace.multi_pairing_ms = stats.multi_pairing_ms;
    trace.gt_multiexp_calls = 1;
    trace.pairing_product_calls = 1;
    trace.intermediate_gt_exponentiations = 0;
    trace.underlying_pairings = combined.pairing_terms.size();
    trace.transcript_ms = equations.timings.transcript_ms;
    trace.batch_inversion_ms = equations.timings.batch_inversion_ms;
    trace.recurrence_ms = equations.timings.recurrence_ms;
    trace.terminal_assembly_ms = equations.timings.terminal_assembly_ms;
    trace.combined_normalize_ms = internal::milliseconds(
        before_normalize, after_normalize);
    trace.gt_msm_pairing_ms = internal::milliseconds(
        after_normalize, after_evaluation);
    trace.total_ms = internal::milliseconds(total_start, internal::Clock::now());
    return trace;
}

} 

DeferredVerificationTrace verify_deferred_with_trace(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof) {
    if (!internal::validate_verification_objects(
            crs, precomputation, statement, proof, true)) return {};
    return evaluate_deferred(crs, precomputation, statement, proof, true);
}

bool verify_deferred(const VmeIbfCRS& crs,
                     const VmeIbfPrecomputation& precomputation,
                     const VmeIbfStatement& statement,
                     const VmeIbfProof& proof) {
    if (!internal::validate_verification_objects(
            crs, precomputation, statement, proof, true)) return false;
    return evaluate_deferred(crs, precomputation, statement, proof, false).accepted;
}

CombinedVerificationTrace verify_deferred_combined_with_trace(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof) {
    if (!internal::validate_verification_objects(
            crs, precomputation, statement, proof, true)) return {};
    return evaluate_combined(crs, precomputation, statement, proof, true);
}

std::optional<ValidatedVerificationInputs> validate_verification_inputs(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof) {
    if (!internal::validate_verification_objects(
            crs, precomputation, statement, proof, true)) return std::nullopt;
    return ValidatedVerificationInputs(crs, precomputation, statement, proof);
}

bool verify_online(const ValidatedVerificationInputs& inputs) {
    return evaluate_combined(inputs.crs(), inputs.precomputation(),
                             inputs.statement(), inputs.proof(), false).accepted;
}

CombinedVerificationTrace verify_online_with_trace(
    const ValidatedVerificationInputs& inputs) {
    return evaluate_combined(inputs.crs(), inputs.precomputation(),
                             inputs.statement(), inputs.proof(), false);
}

bool verify_deferred_combined(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof) {
    const auto inputs = validate_verification_inputs(
        crs, precomputation, statement, proof);
    return inputs && verify_online(*inputs);
}

namespace internal {

CombinedVerificationTrace verify_with_transcript_unchecked(
    const VmeIbfCRS& crs,
    const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement,
    const VmeIbfProof& proof,
    const Digest& transcript_state,
    bool audit_reference) {
    return evaluate_combined(
        crs, precomputation, statement, proof,
        audit_reference, &transcript_state);
}

} 

} 

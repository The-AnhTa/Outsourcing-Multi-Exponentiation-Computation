#include "rexp/verify.hpp"

#include "internal/crypto.hpp"
#include "internal/rexp_helpers.hpp"
#include "internal/rexp_transcript.hpp"

#include <chrono>

namespace rexp {
namespace {

using Clock = std::chrono::steady_clock;

bool verify_core(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const RexpProof& proof,
    bool reference,
    bool already_validated,
    RexpVerifyMetrics* output) {
    const auto start = output ? Clock::now() : Clock::time_point{};
    if (output) *output = RexpVerifyMetrics{};
    try {
        internal::require_bound(params, statement);
        if (!internal::rexp_proof_shape(params, proof)) return false;
        if (!already_validated) {
            internal::validate_rexp_proof_gt(proof, nullptr);
        }
        Digest transcript =
            internal::rexp_initial_transcript(params, statement);
        GT d1 = statement.D1Initial();
#ifdef REXP_ENABLE_PROFILING
        if (output) output->per_dory.reserve(params.d());
#endif
        for (std::size_t round = 0; round < params.d(); ++round) {
            const std::size_t dimension = params.n() >> round;
            const std::size_t half = dimension / 2;
            const RexpRoundMessage message = round
                ? proof.dynamicRoundMessages[round - 1]
                : internal::rexp_initial_round_message(statement);
            const Digest round_digest = internal::rexp_absorb_round(
                transcript, round, dimension, message);
            const Fr rho = ChallengeNonzeroFr(
                round_digest, "REXP-G1-RHO-V1", round);
            const Fr rho_inverse = internal::inverse(rho);
            const DoryStatement dory_statement {
                internal::gt_mul(
                    internal::gt_mul(d1, internal::gt_pow(message.E, rho)),
                    internal::gt_pow(message.F, rho_inverse)),
                internal::gt_mul(
                    message.TL, internal::gt_pow(message.TR, rho)),
                internal::gt_mul(
                    params.X()[round + 1],
                    internal::gt_pow(
                        params.Delta2R()[round], rho_inverse))};
            const Digest dory_input = internal::rexp_enter_dory(
                round_digest, round, half, dory_statement);
            Digest dory_end;
            if (reference) {
                if (!VerifyEmbeddedReference(
                        params.levelCRS(round + 1),
                        params.levelPrecomputation(round + 1),
                        dory_statement, proof.doryProofs[round], dory_input,
                        &dory_end)) {
                    return false;
                }
            } else {
                VerifyMetrics metrics;
                const auto& level_crs = params.levelCRS(round + 1);
                const auto& level_precomp =
                    params.levelPrecomputation(round + 1);
                if (!VerifyEmbeddedDeferred(
                        level_crs, level_precomp, dory_statement,
                        proof.doryProofs[round], dory_input, &dory_end,
                        output ? &metrics : nullptr)) {
                    return false;
                }
                if (output) {
                    ++output->dory_verifications;
                    output->gt_multiexponentiations += metrics.gt_multiexp_calls;
                    output->dory_terminal_pairings += metrics.terminal_pairings;
                    output->actual_gt_bases += metrics.actual_gt_bases;
                    output->transcript_ms += metrics.transcript_ms;
                    output->gt_multiexp_ms += metrics.gt_multiexp_ms;
                    output->dory_pairing_ms += metrics.terminal_pairing_ms;
#ifdef REXP_ENABLE_PROFILING
                    output->per_dory.push_back({
                        round, level_crs.d,
                        metrics.symbolic_atom_insertions,
                        metrics.actual_gt_bases,
                        metrics.coalesced_duplicate_bases,
                        metrics.zero_coefficients_removed,
                        metrics.identity_bases_removed,
                        metrics.gt_multiexp_calls,
                        metrics.terminal_pairings});
#endif
                }
            }
            transcript = internal::rexp_leave_dory(dory_end, round);
            d1 = dory_statement.D1;
        }
        (void)internal::rexp_absorb_final(transcript, proof.R);
        const auto pairing_start = output ? Clock::now() : Clock::time_point{};
        if (params.Lambda()[0].isZero()) return false;
        GT pairing;
        mcl::bn::pairing(pairing, proof.R, params.Lambda()[0]);
        if (output) {
            output->final_pairings = 1;
            output->final_pairing_ms =
                std::chrono::duration<double, std::milli>(
                    Clock::now() - pairing_start).count();
        }
        const bool accepted = pairing == d1;
        if (output) {
            output->total_ms = std::chrono::duration<double, std::milli>(
                Clock::now() - start).count();
        }
        return accepted;
    } catch (...) {
        return false;
    }
}

} // namespace

bool VerifyPrepared(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const RexpProof& proof,
    RexpVerifyMetrics* metrics) {
    return verify_core(params, statement, proof, false, false, metrics);
}

bool VerifyValidatedProof(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const ValidatedRexpProof& proof,
    RexpVerifyMetrics* metrics) {
    if (proof.d() != params.d()) return false;
    return verify_core(
        params, statement, proof.proof(), false, true, metrics);
}

bool VerifyOptimized(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const ValidatedRexpProof& proof) {
    return VerifyValidatedProof(params, statement, proof);
}

bool VerifyReference(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const RexpProof& proof) {
    return verify_core(params, statement, proof, true, false, nullptr);
}

bool Verify(
    const RawRexpCRS& crs,
    const RawRexpStatement& statement,
    const RexpProof& proof,
    RexpVerifyMetrics* metrics) {
    try {
        const PreparedPublicParameters params = PreparePublicParameters(crs);
        const PreparedStatement prepared = PrepareStatement(params, statement);
        return VerifyPrepared(params, prepared, proof, metrics);
    } catch (...) {
        return false;
    }
}

} // namespace rexp

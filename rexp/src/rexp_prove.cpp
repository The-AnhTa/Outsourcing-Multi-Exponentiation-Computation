#include "rexp/prove.hpp"

#include "internal/crypto.hpp"
#include "internal/rexp_helpers.hpp"
#include "internal/rexp_transcript.hpp"

#include <stdexcept>
#include <utility>

namespace rexp {

RexpProof Prove(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement,
    const RexpProverInput& input) {
    internal::require_bound(params, statement);
    if (input.H != statement.H()) {
        throw std::invalid_argument("prover input differs from H");
    }
    RexpProof proof;
    Digest transcript = internal::rexp_initial_transcript(params, statement);
    std::vector<G1> points = statement.H();
    GT d1 = statement.D1Initial();
    for (std::size_t round = 0; round < params.d(); ++round) {
        const std::size_t dimension = params.n() >> round;
        const std::size_t half = dimension / 2;
        RexpRoundMessage message;
        if (round == 0) {
            message = internal::rexp_initial_round_message(statement);
        } else {
            message = {
                internal::pairing_product(
                    points, half, params.Lambda(), 0, half),
                internal::pairing_product(
                    points, 0, params.Lambda(), half, half),
                internal::pairing_product(
                    points, 0, params.Lambda(), 0, half),
                internal::pairing_product(
                    points, half, params.Lambda(), 0, half)};
            proof.dynamicRoundMessages.push_back(message);
        }
        const Digest round_digest = internal::rexp_absorb_round(
            transcript, round, dimension, message);
        const Fr rho = ChallengeNonzeroFr(
            round_digest, "REXP-G1-RHO-V1", round);
        const Fr rho_inverse = internal::inverse(rho);
        std::vector<G1> next(half);
        std::vector<G2> theta(half);
        for (std::size_t i = 0; i < half; ++i) {
            next[i] = internal::g1_add(
                points[i], internal::g1_mul(points[half + i], rho));
            theta[i] = internal::g2_add(
                params.Lambda()[i],
                internal::g2_mul(params.Lambda()[half + i], rho_inverse));
        }
        const DoryStatement dory_statement {
            internal::gt_mul(
                internal::gt_mul(d1, internal::gt_pow(message.E, rho)),
                internal::gt_pow(message.F, rho_inverse)),
            internal::gt_mul(
                message.TL, internal::gt_pow(message.TR, rho)),
            internal::gt_mul(
                params.X()[round + 1],
                internal::gt_pow(params.Delta2R()[round], rho_inverse))};
        const Digest dory_input = internal::rexp_enter_dory(
            round_digest, round, half, dory_statement);
        Digest dory_end;
        proof.doryProofs.push_back(ProveEmbedded(
            params.levelCRS(round + 1), {next, theta}, dory_input,
            &dory_end));
        transcript = internal::rexp_leave_dory(dory_end, round);
        d1 = dory_statement.D1;
        points = std::move(next);
    }
    proof.R = points[0];
    (void)internal::rexp_absorb_final(transcript, proof.R);
    return proof;
}

} // namespace rexp

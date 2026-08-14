#include "vpip_bf/phase1.hpp"
#include "internal/protocol.hpp"

#include <stdexcept>

namespace vpip_bf {

std::vector<Fr> tensor_vector(std::span<const Fr> challenges) {
    std::vector<Fr> weights(1);
    weights[0] = 1;
    for (std::size_t j = challenges.size(); j > 0; --j) {
        const Fr& challenge = challenges[j - 1];
        std::vector<Fr> next;
        next.reserve(weights.size() * 2);
        for (const auto& weight : weights) {
            next.push_back(weight);
            Fr scaled;
            Fr::mul(scaled, weight, challenge);
            next.push_back(scaled);
        }
        weights = std::move(next);
    }
    return weights;
}

Digest compute_statement_digest(const VpipBfCRS& crs,
                                const VpipBfPrecomputation&,
                                const VpipBfStatementInput& input,
                                const GT& commitment) {
    Bytes bytes;
    append_frame(bytes, "vpipbf/statement/v1");
    append_frame(bytes, crs.digest);
    append_frame(bytes, input.digest);
    for (const auto& point : input.X) append_frame(bytes, serialize(point));
    for (const auto& point : crs.H) append_frame(bytes, serialize(point));
    append_frame(bytes, serialize(crs.Lprime));
    append_frame(bytes, serialize(commitment));
    return sha256(bytes);
}

Phase1Result prove_phase1(const VpipBfCRS& crs,
                          const VpipBfPrecomputation& precomputation,
                          const VpipBfStatementInput& input) {
    if (!validate_crs(crs) || !audit_precomputation(crs, precomputation)
        || !validate_statement_input(crs, input))
        throw std::invalid_argument("malformed Phase-I input");

    Phase1Result result;
    result.statement.X = input.X;
    result.statement.C = pairing_product(input.X, crs.H);
    result.statement.digest = compute_statement_digest(
        crs, precomputation, input, result.statement.C);
    Transcript transcript(result.statement.digest);
    result.transcript_start = transcript.digest();
    result.rho.resize(crs.d);
    result.fresh.resize(crs.d + 1);
    std::vector<G1> folded_g = crs.G;
    GT outer = precomputation.pairing_x[0];

    for (std::size_t round = 0; round < crs.d; ++round) {
        const std::size_t dimension = crs.n >> round;
        const std::size_t half = dimension / 2;
        const std::size_t level = round + 1;
        RexpClaims claim;
        if (round == 0) {
            claim = internal::rexp_claim(0, precomputation, {});
        } else {
            claim.E = pairing_product(
                std::span(folded_g).subspan(half, half),
                std::span(crs.H).first(half));
            claim.F = pairing_product(
                std::span(folded_g).first(half),
                std::span(crs.H).subspan(half, half));
            claim.TL = pairing_product(
                std::span(folded_g).first(half),
                std::span(crs.H).first(half));
            claim.TR = pairing_product(
                std::span(folded_g).subspan(half, half),
                std::span(crs.H).first(half));
            result.dynamic_claims.push_back(claim);
        }
        internal::absorb_rexp_claim(
            transcript, round, dimension, claim);
        const Fr rho = transcript.challenge_nonzero(
            "vpipbf/challenge/rexp-r/v1", round);
        const Fr rho_inverse = inverse_nonzero(rho);
        result.rho[round] = rho;

        std::vector<G1> next_g;
        std::vector<G2> theta;
        next_g.reserve(half);
        theta.reserve(half);
        for (std::size_t i = 0; i < half; ++i) {
            next_g.push_back(g1_add(
                folded_g[i], g1_pow(folded_g[half + i], rho)));
            theta.push_back(g2_add(
                crs.H[i], g2_pow(crs.H[half + i], rho_inverse)));
        }
        auto& fresh = result.fresh[level];
        fresh.D0 = gt_mul(gt_mul(outer, gt_pow(claim.E, rho)),
                          gt_pow(claim.F, rho_inverse));
        fresh.D1 = gt_mul(claim.TL, gt_pow(claim.TR, rho));
        fresh.D2 = gt_mul(precomputation.pairing_x[level],
                          gt_pow(precomputation.delta2R[round], rho_inverse));
        fresh.Phi = next_g;
        fresh.Theta = std::move(theta);
        folded_g = std::move(next_g);
        outer = fresh.D1;
    }
    result.R = folded_g.at(0);
    transcript.absorb("vpipbf/rexp-result-R/v1", serialize(result.R));
    result.transcript_after_R = transcript.digest();
    result.r.resize(crs.d);
    for (std::size_t round = 0; round < crs.d; ++round)
        result.r[crs.d - round - 1] = result.rho[round];
    return result;
}

std::vector<Fr> replay_rho(const VpipBfCRS& crs,
                           const VpipBfPrecomputation& precomputation,
                           const VpipBfStatement& statement,
                           std::span<const RexpClaims> dynamic_claims,
                           const G1& result_R, Digest* after_R) {
    if (dynamic_claims.size() != crs.d - 1)
        throw std::invalid_argument("dynamic claim count");
    Transcript transcript(statement.digest);
    std::vector<Fr> rho(crs.d);
    for (std::size_t round = 0; round < crs.d; ++round) {
        const auto claim = internal::rexp_claim(
            round, precomputation, dynamic_claims);
        internal::absorb_rexp_claim(
            transcript, round, crs.n >> round, claim);
        rho[round] = transcript.challenge_nonzero(
            "vpipbf/challenge/rexp-r/v1", round);
    }
    transcript.absorb("vpipbf/rexp-result-R/v1", serialize(result_R));
    if (after_R) *after_R = transcript.digest();
    return rho;
}

} // namespace vpip_bf

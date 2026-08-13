#include "vme_ibf/phase1.hpp"

#include "vme_ibf/group_utils.hpp"
#include "internal/crypto.hpp"
#include "internal/protocol.hpp"

#include <stdexcept>

namespace vme_ibf {
namespace {

void validate_phase1_inputs(const VmeIbfCRS& crs,
                            const VmeIbfPrecomputation& precomputation,
                            const VmeIbfStatementInput& input) {
    if (!validate_crs(crs)
        || !validate_precomputation_shape(crs, precomputation)
        || !validate_precomputation_elements(precomputation)
        || !audit_precomputation(crs, precomputation)
        || input.x.size() != crs.n
        || input.digest != compute_statement_input_digest(crs, input.x))
        throw std::invalid_argument("invalid phase-1 input");
}

RexpClaims dynamic_claim(const VmeIbfCRS& crs,
                         std::span<const G1> current,
                         std::size_t half) {
    RexpClaims claim;
    claim.E = pairing_product(current.subspan(half, half),
                              std::span(crs.H).first(half));
    claim.F = pairing_product(current.first(half),
                              std::span(crs.H).subspan(half, half));
    claim.TL = pairing_product(current.first(half),
                               std::span(crs.H).first(half));
    claim.TR = pairing_product(current.subspan(half, half),
                               std::span(crs.H).first(half));
    return claim;
}

std::vector<G1> fold_g1(std::span<const G1> current, const Fr& rho) {
    const std::size_t half = current.size() / 2;
    std::vector<G1> folded;
    folded.reserve(half);
    for (std::size_t i = 0; i < half; ++i)
        folded.push_back(g1_add(current[i], g1_pow(current[half + i], rho)));
    return folded;
}

std::vector<G2> fold_g2(std::span<const G2> current, const Fr& rho_inverse) {
    const std::size_t half = current.size() / 2;
    std::vector<G2> folded;
    folded.reserve(half);
    for (std::size_t i = 0; i < half; ++i)
        folded.push_back(g2_add(current[i],
                                g2_pow(current[half + i], rho_inverse)));
    return folded;
}

bool valid_claim(const RexpClaims& claim) {
    return internal::valid_gt(claim.E) && internal::valid_gt(claim.F)
        && internal::valid_gt(claim.TL) && internal::valid_gt(claim.TR);
}

} // namespace

std::vector<Fr> tensor_vector(std::span<const Fr> challenges) {
    std::vector<Fr> weights(1);
    weights[0] = 1;
    for (std::size_t j = challenges.size(); j > 0; --j) {
        const Fr& challenge = challenges[j - 1];
        std::vector<Fr> next;
        next.reserve(weights.size() * 2);
        for (const Fr& weight : weights) {
            next.push_back(weight);
            Fr scaled;
            Fr::mul(scaled, weight, challenge);
            next.push_back(scaled);
        }
        weights = std::move(next);
    }
    return weights;
}

Digest compute_statement_digest(const VmeIbfCRS& crs,
                                const VmeIbfStatementInput& input,
                                const G2& x) {
    Bytes bytes;
    append_frame(bytes, "VME.BF.G2/STATEMENT/V2");
    append_frame(bytes, crs.digest);
    append_frame(bytes, input.digest);
    append_frame(bytes, serialize(x));
    return sha256(bytes);
}

Phase1Result prove_phase1(const VmeIbfCRS& crs,
                          const VmeIbfPrecomputation& precomputation,
                          const VmeIbfStatementInput& input) {
    validate_phase1_inputs(crs, precomputation, input);
    Phase1Result result;
    result.statement.x = input.x;
    result.statement.X = g2_multiexp_protocol(crs.H, input.x);
    result.statement.digest = compute_statement_digest(
        crs, input, result.statement.X);
    Transcript transcript(result.statement.digest);
    result.rho.resize(crs.d);
    result.fresh.resize(crs.d + 1);

    std::vector<G1> current = crs.G;
    GT outer = precomputation.pairing_x[0];
    for (std::size_t j = 0; j < crs.d; ++j) {
        const std::size_t dimension = crs.n >> j;
        const std::size_t half = dimension / 2;
        const std::size_t level = j + 1;
        const RexpClaims claim = j == 0
            ? internal::initial_rexp_claim(precomputation)
            : dynamic_claim(crs, current, half);
        if (j > 0) result.dynamic_claims.push_back(claim);
        internal::absorb_rexp_claim(transcript, j, dimension, claim);
        const Fr rho = internal::derive_rho(transcript, j);
        const Fr rho_inverse = inverse_nonzero(rho);
        result.rho[j] = rho;

        auto folded_g1 = fold_g1(current, rho);
        auto folded_g2 = fold_g2(std::span(crs.H).first(dimension), rho_inverse);
        FreshDoryInstance& fresh = result.fresh[level];
        fresh.D0 = gt_mul(gt_mul(outer, gt_pow(claim.E, rho)),
                          gt_pow(claim.F, rho_inverse));
        fresh.D1 = gt_mul(claim.TL, gt_pow(claim.TR, rho));
        fresh.D2 = gt_mul(precomputation.pairing_x[level],
                          gt_pow(precomputation.delta2R[j], rho_inverse));
        fresh.Phi = folded_g1;
        fresh.Theta = std::move(folded_g2);
        current = std::move(folded_g1);
        outer = fresh.D1;
    }
    if (current.size() != 1) throw std::logic_error("REXP did not terminate");
    result.R = current[0];
    internal::absorb_r(transcript, result.R);
    result.transcript_after_R = transcript.digest();
    result.r.resize(crs.d);
    for (std::size_t j = 0; j < crs.d; ++j)
        result.r[crs.d - j - 1] = result.rho[j];
    return result;
}

std::vector<Fr> replay_rho(const VmeIbfCRS& crs,
                           const VmeIbfPrecomputation& precomputation,
                           const VmeIbfStatement& statement,
                           std::span<const RexpClaims> dynamic_claims,
                           const G1& r,
                           Digest* after_r) {
    if (!validate_crs(crs)
        || !validate_precomputation_shape(crs, precomputation)
        || !validate_precomputation_elements(precomputation)
        || !validate_statement_shape(crs, statement)
        || !validate_statement_elements(statement)
        || !validate_statement_digest(crs, statement)
        || dynamic_claims.size() != crs.d - 1 || !valid_g1(r))
        throw std::invalid_argument("invalid REXP replay input");
    for (const auto& claim : dynamic_claims)
        if (!valid_claim(claim))
            throw std::invalid_argument("invalid REXP claim");
    Transcript transcript(statement.digest);
    std::vector<Fr> rho(crs.d);
    for (std::size_t j = 0; j < crs.d; ++j) {
        const RexpClaims claim = j == 0
            ? internal::initial_rexp_claim(precomputation)
            : dynamic_claims[j - 1];
        internal::absorb_rexp_claim(transcript, j, crs.n >> j, claim);
        rho[j] = internal::derive_rho(transcript, j);
    }
    internal::absorb_r(transcript, r);
    if (after_r) *after_r = transcript.digest();
    return rho;
}

} // namespace vme_ibf

#include "vpip_bf/proof.hpp"
#include "vpip_bf/phase1.hpp"
#include "vpip_bf/setup.hpp"

#include "internal/crypto.hpp"
#include "vpip_bf/group_utils.hpp"

namespace vpip_bf {

bool validate_crs_shape(const VpipBfCRS& crs) {
    try {
        return crs.n == internal::dimension_size(crs.d)
            && crs.G.size() == crs.n && crs.H.size() == crs.n;
    } catch (...) {
        return false;
    }
}

bool validate_crs_elements(const VpipBfCRS& crs) {
    for (const auto& point : crs.G)
        if (!valid_g1(point, true)) return false;
    for (const auto& point : crs.H)
        if (!valid_g2(point, true)) return false;
    return valid_g2(crs.Lprime, true);
}

bool validate_crs_digest(const VpipBfCRS& crs) {
    try { return crs.digest == compute_crs_digest(crs); }
    catch (...) { return false; }
}

bool validate_crs(const VpipBfCRS& crs) {
    return validate_crs_shape(crs) && validate_crs_elements(crs)
        && validate_crs_digest(crs);
}

bool validate_precomputation_shape(
    const VpipBfCRS& crs, const VpipBfPrecomputation& precomputation) {
    return validate_crs_shape(crs)
        && precomputation.pairing_x.size() == crs.d + 1
        && precomputation.delta1R.size() == crs.d
        && precomputation.delta2R.size() == crs.d;
}

bool validate_precomputation_elements(const VpipBfPrecomputation& p) {
    for (const auto& value : p.pairing_x)
        if (!internal::valid_gt(value)) return false;
    for (const auto& value : p.delta1R)
        if (!internal::valid_gt(value)) return false;
    for (const auto& value : p.delta2R)
        if (!internal::valid_gt(value)) return false;
    return true;
}

bool audit_precomputation(
    const VpipBfCRS& crs, const VpipBfPrecomputation& p) {
    if (!validate_crs(crs) || !validate_precomputation_shape(crs, p)
        || !validate_precomputation_elements(p)
        || p.digest != compute_precomputation_digest(crs, p)) return false;
    for (std::size_t k = 0; k <= crs.d; ++k) {
        const std::size_t dimension = crs.n >> k;
        if (p.pairing_x[k] != pairing_product(
                std::span(crs.G).first(dimension),
                std::span(crs.H).first(dimension))) return false;
    }
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t dimension = crs.n >> k;
        const std::size_t half = dimension / 2;
        if (p.delta1R[k] != pairing_product(
                std::span(crs.G).subspan(half, half),
                std::span(crs.H).first(half))) return false;
        if (p.delta2R[k] != pairing_product(
                std::span(crs.G).first(half),
                std::span(crs.H).subspan(half, half))) return false;
    }
    return true;
}

bool validate_statement_input(
    const VpipBfCRS& crs, const VpipBfStatementInput& input) {
    if (!validate_crs(crs) || input.X.size() != crs.n) return false;
    for (const auto& point : input.X)
        if (!valid_g1(point)) return false;
    try {
        return input.digest == compute_statement_input_digest(crs, input.X);
    } catch (...) {
        return false;
    }
}

bool validate_statement_shape(
    const VpipBfCRS& crs, const VpipBfStatement& statement) {
    return validate_crs_shape(crs) && statement.X.size() == crs.n;
}

bool validate_statement_elements(const VpipBfStatement& statement) {
    for (const auto& point : statement.X)
        if (!valid_g1(point)) return false;
    return internal::valid_gt(statement.C);
}

bool validate_statement_digest(
    const VpipBfCRS& crs, const VpipBfStatement& statement) {
    try {
        VpipBfStatementInput input{statement.X, {}};
        input.digest = compute_statement_input_digest(crs, input.X);
        return statement.digest
            == compute_statement_digest(crs, VpipBfPrecomputation{},
                                        input, statement.C);
    } catch (...) {
        return false;
    }
}

bool validate_proof_shape(const VpipBfCRS& crs, const VpipBfProof& proof) {
    return validate_crs_shape(crs)
        && proof.rexp_claims.size() == crs.d - 1
        && proof.dory_folds.size() == crs.d
        && proof.batch_U.size() == crs.d;
}

bool validate_proof_elements(const VpipBfProof& proof) {
    if (!valid_g1(proof.R) || !valid_g1(proof.PhiFinal)
        || !valid_g2(proof.ThetaFinal)) return false;
    for (const auto& claim : proof.rexp_claims)
        if (!internal::valid_gt(claim.E) || !internal::valid_gt(claim.F)
            || !internal::valid_gt(claim.TL)
            || !internal::valid_gt(claim.TR)) return false;
    for (const auto& fold : proof.dory_folds)
        if (!internal::valid_gt(fold.D1L) || !internal::valid_gt(fold.D1R)
            || !internal::valid_gt(fold.D2L)
            || !internal::valid_gt(fold.D2R)
            || !internal::valid_gt(fold.W1)
            || !internal::valid_gt(fold.W2)) return false;
    for (const auto& value : proof.batch_U)
        if (!internal::valid_gt(value)) return false;
    return true;
}

} // namespace vpip_bf

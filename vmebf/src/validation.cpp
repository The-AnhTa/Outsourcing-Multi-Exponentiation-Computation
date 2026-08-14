#include "vme_ibf/setup.hpp"

#include "vme_ibf/group_utils.hpp"
#include "internal/crypto.hpp"

namespace vme_ibf {
namespace {

Digest statement_digest(const VmeIbfCRS& crs, const VmeIbfStatement& statement) {
    const Digest input_digest = compute_statement_input_digest(crs, statement.x);
    Bytes bytes;
    append_frame(bytes, "VME.BF.G2/STATEMENT/V2");
    append_frame(bytes, crs.digest);
    append_frame(bytes, input_digest);
    append_frame(bytes, serialize(statement.X));
    return sha256(bytes);
}

} 

bool validate_crs_shape(const VmeIbfCRS& crs) {
    try {
        return crs.n == internal::dimension_size(crs.d)
            && crs.G.size() == crs.n && crs.H.size() == crs.n;
    } catch (...) {
        return false;
    }
}

bool validate_crs_elements(const VmeIbfCRS& crs) {
    for (const auto& point : crs.G)
        if (!valid_g1(point, true)) return false;
    for (const auto& point : crs.H)
        if (!valid_g2(point, true)) return false;
    return valid_g1(crs.L, true) && valid_g2(crs.Lprime, true);
}

bool validate_crs_digest(const VmeIbfCRS& crs) {
    try { return crs.digest == compute_crs_digest(crs); }
    catch (...) { return false; }
}

bool validate_crs(const VmeIbfCRS& crs) {
    return validate_crs_shape(crs) && validate_crs_elements(crs)
        && validate_crs_digest(crs);
}

bool validate_precomputation_shape(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation) {
    return validate_crs_shape(crs)
        && precomputation.pairing_x.size() == crs.d + 1
        && precomputation.delta1R.size() == crs.d
        && precomputation.delta2R.size() == crs.d;
}

bool validate_precomputation_elements(const VmeIbfPrecomputation& p) {
    for (const auto& value : p.pairing_x)
        if (!internal::valid_gt(value)) return false;
    for (const auto& value : p.delta1R)
        if (!internal::valid_gt(value)) return false;
    for (const auto& value : p.delta2R)
        if (!internal::valid_gt(value)) return false;
    return internal::valid_gt(p.pairing_LLprime);
}

bool audit_precomputation(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& p) {
    if (!validate_crs(crs) || !validate_precomputation_shape(crs, p)
        || !validate_precomputation_elements(p)) return false;
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
    GT pairing;
    mcl::bn::pairing(pairing, crs.L, crs.Lprime);
    return pairing == p.pairing_LLprime;
}

bool validate_statement_shape(
    const VmeIbfCRS& crs, const VmeIbfStatement& statement) {
    return validate_crs_shape(crs) && statement.x.size() == crs.n;
}

bool validate_statement_elements(const VmeIbfStatement& statement) {
    return valid_g2(statement.X);
}

bool validate_statement_digest(
    const VmeIbfCRS& crs, const VmeIbfStatement& statement) {
    try { return statement.digest == statement_digest(crs, statement); }
    catch (...) { return false; }
}

} 

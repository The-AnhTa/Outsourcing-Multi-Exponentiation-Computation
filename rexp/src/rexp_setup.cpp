#include "rexp/setup.hpp"

#include "internal/crypto.hpp"
#include "internal/rexp_transcript.hpp"

#include <stdexcept>

namespace rexp {
namespace {

DoryCRS make_level_crs(
    const PreparedPublicParameters& params, std::size_t level) {
    if (level > params.d()) throw std::invalid_argument("invalid generator level");
    DoryCRS result;
    result.d = params.d() - level;
    result.n = params.n() >> level;
    result.Gamma.assign(
        params.Gamma().begin(), params.Gamma().begin() + result.n);
    result.Lambda.assign(
        params.Lambda().begin(), params.Lambda().begin() + result.n);
    result.digest = ComputeDoryCRSDigest(result);
    return result;
}

DoryPrecomputation make_level_precomputation(
    const PreparedPublicParameters& params,
    std::size_t level,
    const Digest& digest) {
    if (level > params.d()) throw std::invalid_argument("invalid precomp level");
    DoryPrecomputation result;
    result.X.assign(params.X().begin() + level, params.X().end());
    result.Delta1R.assign(
        params.Delta1R().begin() + level, params.Delta1R().end());
    result.Delta2R.assign(
        params.Delta2R().begin() + level, params.Delta2R().end());
    result.crs_digest = digest;
    result.pairing_product_terms = 4 * (params.n() >> level) - 3;
    return result;
}

} 

RawRexpCRS GenerateRawCRS(std::size_t d, std::string_view seed) {
    const DoryCRS crs = GenerateDoryCRS(d, seed);
    return {crs.d, crs.n, crs.Gamma, crs.Lambda};
}

PreparedPublicParameters PreparePublicParameters(const RawRexpCRS& raw_crs) {
    initialize();
    internal::check_dimension(raw_crs.d, raw_crs.n, "Rexp");
    if (raw_crs.Gamma.size() != raw_crs.n
        || raw_crs.Lambda.size() != raw_crs.n) {
        throw std::invalid_argument("raw CRS vector length mismatch");
    }
    for (const auto& point : raw_crs.Gamma) {
        internal::require_point(point, true, "CRS G1 point");
    }
    for (const auto& point : raw_crs.Lambda) {
        internal::require_point(point, true, "CRS G2 point");
    }
    PreparedPublicParameters result;
    result.d_ = raw_crs.d;
    result.n_ = raw_crs.n;
    result.Gamma_ = raw_crs.Gamma;
    result.Lambda_ = raw_crs.Lambda;
    result.digest_ = internal::rexp_crs_digest(raw_crs);
    result.X_.reserve(result.d_ + 1);
    for (std::size_t level = 0; level <= result.d_; ++level) {
        const std::size_t dimension = result.n_ >> level;
        result.X_.push_back(internal::pairing_product(
            result.Gamma_, 0, result.Lambda_, 0, dimension));
        result.pairingProductTerms_ += dimension;
    }
    result.Delta1R_.reserve(result.d_);
    result.Delta2R_.reserve(result.d_);
    for (std::size_t level = 0; level < result.d_; ++level) {
        const std::size_t dimension = result.n_ >> level;
        const std::size_t half = dimension / 2;
        result.Delta1R_.push_back(internal::pairing_product(
            result.Gamma_, half, result.Lambda_, 0, half));
        result.Delta2R_.push_back(internal::pairing_product(
            result.Gamma_, 0, result.Lambda_, half, half));
        result.pairingProductTerms_ += 2 * half;
    }
    result.levelCRS_.resize(result.d_ + 1);
    result.levelPrecomputation_.resize(result.d_ + 1);
    for (std::size_t level = 1; level <= result.d_; ++level) {
        result.levelCRS_[level] = make_level_crs(result, level);
        result.levelPrecomputation_[level] = make_level_precomputation(
            result, level, result.levelCRS_[level].digest);
    }
    return result;
}

RawRexpStatement GenerateRawStatement(
    const PreparedPublicParameters& params) {
    initialize();
    RawRexpStatement result;
    result.H.reserve(params.n());
    G1 base;
    static constexpr char domain[] = "REXP-G1-PUBLIC-H-BASE-V1";
    mcl::bn::hashAndMapToG1(base, domain, sizeof(domain) - 1);
    internal::require_point(base, true, "statement base");
    for (std::size_t i = 0; i < params.n(); ++i) {
        Fr scalar;
        do {
            scalar.setByCSPRNG();
        } while (scalar.isZero());
        result.H.push_back(internal::g1_mul(base, scalar));
    }
    return result;
}

PreparedStatement PrepareStatement(
    const PreparedPublicParameters& params,
    const RawRexpStatement& raw_statement) {
    if (raw_statement.H.size() != params.n()) {
        throw std::invalid_argument("raw statement length mismatch");
    }
    for (const auto& point : raw_statement.H) {
        internal::require_point(point, false, "statement G1 point");
    }
    PreparedStatement result;
    result.H_ = raw_statement.H;
    result.crsDigest_ = params.digest();
    result.D1Initial_ = internal::pairing_product(
        result.H_, 0, params.Lambda(), 0, params.n());
    if (params.d() > 0) {
        const std::size_t half = params.n() / 2;
        result.E0_ = internal::pairing_product(
            result.H_, half, params.Lambda(), 0, half);
        result.F0_ = internal::pairing_product(
            result.H_, 0, params.Lambda(), half, half);
        result.TL0_ = internal::pairing_product(
            result.H_, 0, params.Lambda(), 0, half);
        result.TR0_ = result.E0_;
    } else {
        result.E0_.setOne();
        result.F0_.setOne();
        result.TL0_ = result.D1Initial_;
        result.TR0_.setOne();
    }
    result.digest_ = internal::rexp_statement_digest(
        params, result.H_, result.D1Initial_, result.E0_, result.F0_,
        result.TL0_, result.TR0_);
    return result;
}

RexpSetupResult Setup(std::size_t d, std::string_view seed) {
    RexpSetupResult result;
    result.rawCRS = GenerateRawCRS(d, seed);
    result.params = PreparePublicParameters(result.rawCRS);
    result.rawStatement = GenerateRawStatement(result.params);
    result.statement = PrepareStatement(result.params, result.rawStatement);
    result.proverInput.H = result.rawStatement.H;
    return result;
}

} 

#pragma once

#include "rexp/dory_verify.hpp"

#include <cstddef>
#include <string_view>
#include <vector>

namespace rexp {

struct RawRexpCRS {
    std::size_t d = 0, n = 0;
    std::vector<G1> Gamma;
    std::vector<G2> Lambda;
};

struct RawRexpStatement {
    std::vector<G1> H;
};


struct LegacyRexpStatement {
    std::vector<G1> H;
    GT D1Initial, E0, F0, TL0, TR0;
    Digest digest{};
};



struct LegacyRexpPrecomputation {
    std::vector<GT> X, Delta1R, Delta2R;
    Digest crsDigest{};
};

class PreparedPublicParameters {
public:
    std::size_t d() const { return d_; }
    std::size_t n() const { return n_; }
    const std::vector<G1>& Gamma() const { return Gamma_; }
    const std::vector<G2>& Lambda() const { return Lambda_; }
    const Digest& digest() const { return digest_; }
    const std::vector<GT>& X() const { return X_; }
    const std::vector<GT>& Delta1R() const { return Delta1R_; }
    const std::vector<GT>& Delta2R() const { return Delta2R_; }
    std::size_t pairingProductTerms() const { return pairingProductTerms_; }
    const DoryCRS& levelCRS(std::size_t level) const {
        return levelCRS_.at(level);
    }
    const DoryPrecomputation& levelPrecomputation(std::size_t level) const {
        return levelPrecomputation_.at(level);
    }

private:
    std::size_t d_ = 0, n_ = 0, pairingProductTerms_ = 0;
    std::vector<G1> Gamma_;
    std::vector<G2> Lambda_;
    std::vector<GT> X_, Delta1R_, Delta2R_;
    std::vector<DoryCRS> levelCRS_;
    std::vector<DoryPrecomputation> levelPrecomputation_;
    Digest digest_{};
    friend PreparedPublicParameters PreparePublicParameters(const RawRexpCRS&);
};

class PreparedStatement {
public:
    const std::vector<G1>& H() const { return H_; }
    const GT& D1Initial() const { return D1Initial_; }
    const GT& E0() const { return E0_; }
    const GT& F0() const { return F0_; }
    const GT& TL0() const { return TL0_; }
    const GT& TR0() const { return TR0_; }
    const Digest& digest() const { return digest_; }
    const Digest& crsDigest() const { return crsDigest_; }

private:
    std::vector<G1> H_;
    GT D1Initial_, E0_, F0_, TL0_, TR0_;
    Digest digest_{}, crsDigest_{};
    friend PreparedStatement PrepareStatement(
        const PreparedPublicParameters&, const RawRexpStatement&);
};

struct RexpProverInput {
    std::vector<G1> H;
};

struct RexpSetupResult {
    RawRexpCRS rawCRS;
    RawRexpStatement rawStatement;
    PreparedPublicParameters params;
    PreparedStatement statement;
    RexpProverInput proverInput;
};

} 

#pragma once

#include "rexp/dory_verify.hpp"

#include <string_view>

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

struct RexpProverInput { std::vector<G1> H; };
struct RexpRoundMessage { GT E, F, TL, TR; };
struct RexpProof {
    std::vector<RexpRoundMessage> dynamicRoundMessages;
    std::vector<DoryProof> doryProofs;
    G1 R;
};

struct RexpProofValidationMetrics {
    std::size_t gt_elements_checked = 0;
    double canonical_decode_ms = 0;
    double gt_subgroup_validation_ms = 0;
};

class ValidatedRexpProof {
public:
    const RexpProof& proof() const { return proof_; }
    std::size_t d() const { return d_; }
private:
    RexpProof proof_;
    std::size_t d_ = 0;
    friend ValidatedRexpProof ValidateRexpProof(
        const RexpProof&, std::size_t, RexpProofValidationMetrics*);
    friend ValidatedRexpProof DeserializeValidatedRexpProofWire(
        const std::vector<std::uint8_t>&, std::size_t,
        RexpProofValidationMetrics*);
};

struct RexpSetupResult {
    RawRexpCRS rawCRS;
    RawRexpStatement rawStatement;
    PreparedPublicParameters params;
    PreparedStatement statement;
    RexpProverInput proverInput;
};

struct RexpVerifyMetrics {
    std::size_t dory_verifications = 0, gt_multiexponentiations = 0;
    std::size_t dory_terminal_pairings = 0, final_pairings = 0;
    std::size_t actual_gt_bases = 0;
    double transcript_ms = 0, gt_multiexp_ms = 0;
    double dory_pairing_ms = 0, final_pairing_ms = 0, total_ms = 0;
#ifdef REXP_ENABLE_PROFILING
    struct PerDoryProfile {
        std::size_t round = 0, dimension = 0;
        std::size_t raw_gt_terms = 0, normalized_gt_terms = 0;
        std::size_t duplicate_bases_merged = 0;
        std::size_t zero_terms_removed = 0, identity_terms_removed = 0;
        std::size_t gt_multiexp_calls = 0, terminal_pairings = 0;
    };
    std::vector<PerDoryProfile> per_dory;
#endif
};

RawRexpCRS GenerateRawCRS(std::size_t d, std::string_view crs_seed);
RawRexpStatement GenerateRawStatement(const PreparedPublicParameters&);
PreparedPublicParameters PreparePublicParameters(const RawRexpCRS&);
PreparedStatement PrepareStatement(
    const PreparedPublicParameters&, const RawRexpStatement&);
RexpSetupResult Setup(std::size_t d, std::string_view crs_seed);

RexpProof Prove(const PreparedPublicParameters&, const PreparedStatement&,
                const RexpProverInput&);
bool VerifyPrepared(const PreparedPublicParameters&, const PreparedStatement&,
                    const RexpProof&, RexpVerifyMetrics* metrics = nullptr);
bool VerifyValidatedProof(
    const PreparedPublicParameters&, const PreparedStatement&,
    const ValidatedRexpProof&, RexpVerifyMetrics* metrics = nullptr);
bool VerifyOptimized(
    const PreparedPublicParameters&, const PreparedStatement&,
    const ValidatedRexpProof&);
bool VerifyReference(const PreparedPublicParameters&, const PreparedStatement&,
                     const RexpProof&);
bool Verify(const RawRexpCRS&, const RawRexpStatement&, const RexpProof&,
            RexpVerifyMetrics* metrics = nullptr);

std::vector<std::uint8_t> SerializeRexpCRS(const RawRexpCRS&);
std::vector<std::uint8_t> SerializeRexpPrecomputation(
    const PreparedPublicParameters&);
std::vector<std::uint8_t> SerializeRexpStatement(const RawRexpStatement&);
std::vector<std::uint8_t> SerializePreparedStatement(const PreparedStatement&);
std::vector<std::uint8_t> SerializeRexpProof(const RexpProof&, std::size_t d);
std::vector<std::uint8_t> SerializeRexpCRSWire(const RawRexpCRS&);
std::vector<std::uint8_t> SerializeRexpStatementWire(
    const RawRexpStatement&, std::size_t n);
std::vector<std::uint8_t> SerializeRexpProofWire(
    const RexpProof&, std::size_t d);
RawRexpCRS DeserializeRexpCRSWire(const std::vector<std::uint8_t>&);
RawRexpStatement DeserializeRexpStatementWire(
    const std::vector<std::uint8_t>&, std::size_t expected_n);
RexpProof DeserializeRexpProofWire(
    const std::vector<std::uint8_t>&, std::size_t expected_d,
    RexpProofValidationMetrics* metrics = nullptr);
ValidatedRexpProof DeserializeValidatedRexpProofWire(
    const std::vector<std::uint8_t>&, std::size_t expected_d,
    RexpProofValidationMetrics* metrics = nullptr);
bool IsValidGTSubgroup(const GT&);
ValidatedRexpProof ValidateRexpProof(
    const RexpProof&, std::size_t d,
    RexpProofValidationMetrics* metrics = nullptr);

}

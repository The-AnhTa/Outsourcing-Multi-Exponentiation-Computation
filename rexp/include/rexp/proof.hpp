#pragma once

#include "rexp/types.hpp"

namespace rexp {

struct RexpProofValidationMetrics;

struct RexpRoundMessage {
    GT E, F, TL, TR;
};

struct RexpProof {
    std::vector<RexpRoundMessage> dynamicRoundMessages;
    std::vector<DoryProof> doryProofs;
    G1 R;
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
    friend ValidatedRexpProof ValidateRexpProof(
        RexpProof&&, std::size_t, RexpProofValidationMetrics*);
    friend ValidatedRexpProof DeserializeValidatedRexpProofWire(
        const std::vector<std::uint8_t>&, std::size_t,
        RexpProofValidationMetrics*);
};

} 

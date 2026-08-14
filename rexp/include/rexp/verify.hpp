#pragma once

#include "rexp/metrics.hpp"
#include "rexp/proof.hpp"

namespace rexp {

bool VerifyPrepared(
    const PreparedPublicParameters&, const PreparedStatement&,
    const RexpProof&, RexpVerifyMetrics* metrics = nullptr);
bool VerifyValidatedProof(
    const PreparedPublicParameters&, const PreparedStatement&,
    const ValidatedRexpProof&, RexpVerifyMetrics* metrics = nullptr);
bool VerifyOptimized(
    const PreparedPublicParameters&, const PreparedStatement&,
    const ValidatedRexpProof&);
bool VerifyReference(
    const PreparedPublicParameters&, const PreparedStatement&, const RexpProof&);
bool Verify(
    const RawRexpCRS&, const RawRexpStatement&, const RexpProof&,
    RexpVerifyMetrics* metrics = nullptr);

bool IsValidGTSubgroup(const GT&);
ValidatedRexpProof ValidateRexpProof(
    const RexpProof&, std::size_t d,
    RexpProofValidationMetrics* metrics = nullptr);
ValidatedRexpProof ValidateRexpProof(
    RexpProof&&, std::size_t d,
    RexpProofValidationMetrics* metrics = nullptr);

} 

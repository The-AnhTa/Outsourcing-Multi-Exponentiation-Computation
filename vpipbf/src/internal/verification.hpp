#pragma once
#include "vpip_bf/verify_reference.hpp"
namespace vpip_bf {
bool validate_verification_inputs(const VpipBfCRS&,
    const VpipBfPrecomputation&, const VpipBfStatement&,
    const VpipBfProof&);
ReferenceVerificationTrace verify_core_unchecked(const VpipBfCRS&,
    const VpipBfPrecomputation&, const VpipBfStatement&,
    const VpipBfProof&, OnlineTimingBreakdown* = nullptr);
ReferenceVerificationTrace verify_core_symbolic_unchecked(const VpipBfCRS&,
    const VpipBfPrecomputation&, const VpipBfStatement&,
    const VpipBfProof&, OnlineTimingBreakdown* = nullptr);
}

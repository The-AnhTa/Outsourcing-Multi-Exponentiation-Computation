#pragma once

#include "vme_ibf/verify_online.hpp"

namespace vme_ibf::internal {

CombinedVerificationTrace verify_with_transcript_unchecked(
    const VmeIbfCRS&, const VmeIbfPrecomputation&,
    const VmeIbfStatement&, const VmeIbfProof&,
    const Digest& transcript_state, bool audit_reference = false);

} 

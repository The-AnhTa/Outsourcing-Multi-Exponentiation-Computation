#pragma once
#include "vpip_bf/phase1.hpp"

namespace vpip_bf {
Phase2Result prove_phase2(const VpipBfCRS&, const VpipBfPrecomputation&, const Phase1Result&);
}

#pragma once
#include "vme_ibf/phase1.hpp"

namespace vme_ibf {
Phase2Result prove_phase2(const VmeIbfCRS&, const VmeIbfPrecomputation&, const Phase1Result&);
}

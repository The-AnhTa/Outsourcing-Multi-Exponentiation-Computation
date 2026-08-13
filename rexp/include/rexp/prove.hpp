#pragma once

#include "rexp/proof.hpp"

namespace rexp {

RexpProof Prove(
    const PreparedPublicParameters&, const PreparedStatement&,
    const RexpProverInput&);

} // namespace rexp

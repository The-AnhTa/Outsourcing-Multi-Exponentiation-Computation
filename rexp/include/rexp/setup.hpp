#pragma once

#include "rexp/types.hpp"

namespace rexp {

RawRexpCRS GenerateRawCRS(std::size_t d, std::string_view crs_seed);
RawRexpStatement GenerateRawStatement(const PreparedPublicParameters&);
PreparedPublicParameters PreparePublicParameters(const RawRexpCRS&);
PreparedStatement PrepareStatement(
    const PreparedPublicParameters&, const RawRexpStatement&);
RexpSetupResult Setup(std::size_t d, std::string_view crs_seed);

} // namespace rexp

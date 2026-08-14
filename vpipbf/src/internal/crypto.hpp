#pragma once

#include "vpip_bf/types.hpp"

namespace vpip_bf::internal {
std::size_t dimension_size(std::size_t d);
bool canonical_fr(const Fr& value);
bool valid_gt(const GT& value);
}

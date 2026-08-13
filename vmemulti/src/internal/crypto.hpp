#pragma once

#include "vme_ibf/types.hpp"

#include <chrono>
#include <cstddef>

namespace vme_ibf::internal {

std::size_t dimension_size(std::size_t d);
Fr fr_one();
Fr fr_minus_one();
Fr fr_multiply(const Fr& left, const Fr& right);
bool valid_gt(const GT& value);

using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start, Clock::time_point end);

} // namespace vme_ibf::internal

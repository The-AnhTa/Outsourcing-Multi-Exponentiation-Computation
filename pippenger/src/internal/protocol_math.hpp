#pragma once

#include "pinkas/pinkas.hpp"

#include <cstddef>
#include <vector>

namespace pinkas::internal {

Group add(const Group& left, const Group& right);
Scalar scalar_add(const Scalar& left, const Scalar& right);
Group horner_reconstruct(const std::vector<Group>& row);

Group pippenger_msm(
    const std::vector<Group>& points,
    const std::vector<Scalar>& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length);

std::size_t ceil_log2(std::size_t value);

}  

#pragma once

#include "internal/constants.hpp"
#include "pippenger/pippenger.hpp"

#include <cstddef>
#include <vector>

namespace pippenger::internal {

Group add(const Group& left, const Group& right);
Group double_point(const Group& point);

CanonicalScalar canonical_bytes(const Scalar& scalar);

bool fits_unsigned_bits(
    const CanonicalScalar& scalar,
    std::size_t exponent_bit_length);

std::size_t window_digit(
    const CanonicalScalar& scalar,
    std::size_t shift,
    std::size_t width,
    std::size_t exponent_bit_length);

Group evaluate_pippenger_row(
    const std::vector<Group>& bases,
    const std::vector<CanonicalScalar>& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length,
    std::size_t window_count,
    std::size_t bucket_count);

}  

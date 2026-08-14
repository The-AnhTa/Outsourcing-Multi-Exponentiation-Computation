#pragma once

#include <mcl/bn.hpp>

#include <cstddef>
#include <vector>

namespace pippenger {

using Scalar = mcl::bn::Fr;
using Group = mcl::bn::G2;
using ScalarMatrix = std::vector<std::vector<Scalar>>;

void initialize();

std::size_t scalar_bit_length();
std::size_t maximum_window_width() noexcept;
std::size_t number_of_windows(std::size_t window_width);

std::vector<Group> multi_pippenger(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    std::size_t window_width);

std::vector<Group> multi_pippenger_bounded(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length);

std::vector<Group> naive_multi_exponentiation(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars);

}  // namespace pippenger

#pragma once
#include "rexpbf/types.hpp"
#include <span>

namespace rexpbf {
void initialize_bn254();
GT pairing_product(std::span<const G1> left, std::span<const G2> right);
}

#pragma once
#include "rexpbf/types.hpp"
#include <span>

namespace rexpbf {
GT gt_multiexp_pippenger(std::span<const GT> bases, std::span<const Fr> scalars);
GT gt_multiexp_naive(std::span<const GT> bases, std::span<const Fr> scalars);
}

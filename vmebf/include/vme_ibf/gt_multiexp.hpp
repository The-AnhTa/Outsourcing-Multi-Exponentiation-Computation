#pragma once
#include "vme_ibf/types.hpp"
#include <span>
namespace vme_ibf {
GT gt_multiexp_reference(std::span<const GT>,std::span<const Fr>);
GT gt_multiexp_pippenger(std::span<const GT>,std::span<const Fr>);
void set_gt_multiexp_window_override_for_benchmark(std::size_t);
}

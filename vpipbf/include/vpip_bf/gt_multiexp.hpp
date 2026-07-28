#pragma once
#include "vpip_bf/types.hpp"
#include <span>
namespace vpip_bf {
GT gt_multiexp_reference(std::span<const GT>,std::span<const Fr>);
GT gt_multiexp_native(std::span<const GT>,std::span<const Fr>);
GT gt_multiexp_pippenger(std::span<const GT>,std::span<const Fr>);
}

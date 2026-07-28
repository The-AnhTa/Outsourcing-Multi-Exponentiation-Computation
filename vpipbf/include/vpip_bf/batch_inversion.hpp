#pragma once
#include "vpip_bf/types.hpp"
#include <span>
namespace vpip_bf { std::vector<Fr> batch_invert_nonzero(std::span<const Fr> values); }

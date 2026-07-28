#pragma once
#include "vme_ibf/types.hpp"
#include <span>
namespace vme_ibf { std::vector<Fr> batch_invert_nonzero(std::span<const Fr> values); }

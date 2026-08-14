#include "internal/validation.hpp"

#include <stdexcept>

namespace pippenger::internal {

void validate_msm_dimensions(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  if (bases.empty()) {
    throw std::invalid_argument("n must be at least 1");
  }
  if (scalars.empty()) {
    throw std::invalid_argument("k must be at least 1");
  }
  for (const auto& row : scalars) {
    if (row.size() != bases.size()) {
      throw std::invalid_argument("every scalar row must contain n elements");
    }
  }
}

}  // namespace pippenger::internal

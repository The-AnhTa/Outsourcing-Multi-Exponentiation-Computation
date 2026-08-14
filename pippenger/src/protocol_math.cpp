#include "internal/protocol_math.hpp"

#include <stdexcept>

namespace pinkas::internal {

Group add(const Group& left, const Group& right) {
  Group result;
  Group::add(result, left, right);
  return result;
}

Scalar scalar_add(const Scalar& left, const Scalar& right) {
  Scalar result;
  Scalar::add(result, left, right);
  return result;
}

Group horner_reconstruct(const std::vector<Group>& row) {
  if (row.empty()) throw std::invalid_argument("cannot reconstruct empty row");
  Group result = row.back();
  for (std::size_t bit = row.size() - 1; bit-- > 0;) {
    Group doubled;
    Group::dbl(doubled, result);
    result = add(doubled, row[bit]);
  }
  return result;
}

Group pippenger_msm(
    const std::vector<Group>& points,
    const std::vector<Scalar>& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length) {
  ScalarMatrix one_instance(1, scalars);
  return pippenger::multi_pippenger_bounded(
             points, one_instance, window_width, exponent_bit_length)
      .front();
}

std::size_t ceil_log2(std::size_t value) {
  if (value == 0) throw std::invalid_argument("log2 of zero");
  std::size_t result = 0;
  for (--value; value != 0; value >>= 1) ++result;
  return result;
}

}  // namespace pinkas::internal

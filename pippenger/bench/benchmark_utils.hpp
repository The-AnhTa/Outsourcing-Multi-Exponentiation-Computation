#pragma once

#include "pippenger/pippenger.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace benchmark_utils {

inline std::size_t number(const char* text, std::string_view name) {
  const std::string value(text);
  if (value.empty() || value.front() == '-') {
    throw std::invalid_argument("invalid " + std::string(name));
  }
  std::size_t parsed = 0;
  const unsigned long long result = std::stoull(value, &parsed);
  if (parsed != value.size() ||
      result > std::numeric_limits<std::size_t>::max()) {
    throw std::invalid_argument("invalid " + std::string(name));
  }
  return static_cast<std::size_t>(result);
}

inline pippenger::Scalar random_nonzero_scalar() {
  pippenger::Scalar scalar;
  do {
    scalar.setByCSPRNG();
  } while (scalar.isZero());
  return scalar;
}

inline std::vector<pippenger::Group> random_bases(
    std::size_t count,
    std::string_view domain) {
  pippenger::Group generator;
  mcl::bn::hashAndMapToG2(generator, domain.data(), domain.size());
  std::vector<pippenger::Group> bases(count);
  for (auto& base : bases) {
    pippenger::Group::mul(base, generator, random_nonzero_scalar());
  }
  return bases;
}

inline pippenger::ScalarMatrix random_scalars(
    std::size_t rows,
    std::size_t columns) {
  pippenger::ScalarMatrix scalars(
      rows, std::vector<pippenger::Scalar>(columns));
  for (auto& row : scalars) {
    for (auto& scalar : row) scalar.setByCSPRNG();
  }
  return scalars;
}

}  // namespace benchmark_utils

#include "pippenger/pippenger.hpp"

#include "internal/constants.hpp"
#include "internal/pippenger_core.hpp"
#include "internal/validation.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

namespace pippenger {
namespace {

using internal::CanonicalScalar;

}

std::size_t number_of_windows(std::size_t window_width) {
  const std::size_t bit_length = scalar_bit_length();
  if (window_width == 0 || window_width > bit_length ||
      window_width > internal::kMaximumWindowWidth) {
    throw std::invalid_argument("window width exceeds the supported range");
  }
  return (bit_length + window_width - 1) / window_width;
}

std::size_t maximum_window_width() noexcept {
  return internal::kMaximumWindowWidth;
}

std::vector<Group> multi_pippenger(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    std::size_t window_width) {
  return multi_pippenger_bounded(
      bases, scalars, window_width, scalar_bit_length());
}

std::vector<Group> multi_pippenger_bounded(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length) {
  initialize();
  internal::validate_msm_dimensions(bases, scalars);

  if (exponent_bit_length == 0 ||
      exponent_bit_length > scalar_bit_length()) {
    throw std::invalid_argument(
        "exponent bit length must fit the scalar field");
  }
  if (window_width == 0 || window_width > exponent_bit_length ||
      window_width > maximum_window_width()) {
    throw std::invalid_argument(
        "window width exceeds the supported exponent range");
  }
  const std::size_t windows =
      (exponent_bit_length + window_width - 1) / window_width;
  const std::size_t bucket_count =
      (std::size_t{1} << window_width) - 1;

  if (bases.size() > std::vector<CanonicalScalar>{}.max_size() ||
      scalars.size() > std::vector<std::vector<CanonicalScalar>>{}.max_size() ||
      bucket_count > std::vector<Group>{}.max_size()) {
    throw std::length_error("MSM dimensions exceed allocation limits");
  }

  std::vector<std::vector<CanonicalScalar>> canonical(
      scalars.size(), std::vector<CanonicalScalar>(bases.size()));
  for (std::size_t j = 0; j < scalars.size(); ++j) {
    for (std::size_t i = 0; i < bases.size(); ++i) {
      canonical[j][i] = internal::canonical_bytes(scalars[j][i]);
      if (!internal::fits_unsigned_bits(
              canonical[j][i], exponent_bit_length)) {
        throw std::invalid_argument(
            "scalar exceeds the declared exponent bit length");
      }
    }
  }
  std::vector<Group> results(scalars.size());
  for (std::size_t j = 0; j < scalars.size(); ++j) {
    results[j] = internal::evaluate_pippenger_row(
        bases, canonical[j], window_width, exponent_bit_length, windows,
        bucket_count);
  }

  return results;
}

std::vector<Group> naive_multi_exponentiation(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  initialize();
  internal::validate_msm_dimensions(bases, scalars);

  std::vector<Group> results(scalars.size());
  for (std::size_t j = 0; j < scalars.size(); ++j) {
    results[j].clear();
    for (std::size_t i = 0; i < bases.size(); ++i) {
      Group term;
      Group::mul(term, bases[i], scalars[j][i]);
      results[j] = internal::add(results[j], term);
    }
  }
  return results;
}

}

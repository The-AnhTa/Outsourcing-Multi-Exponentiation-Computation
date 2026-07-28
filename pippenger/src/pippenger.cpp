#include "pippenger.hpp"

#include <array>
#include <cstdint>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace pippenger {
namespace {

constexpr std::size_t kMaxScalarBytes = 32;
using CanonicalScalar = std::array<std::uint8_t, kMaxScalarBytes>;

Group add(const Group& a, const Group& b) {
  Group result;
  Group::add(result, a, b);
  return result;
}

Group double_point(const Group& point) {
  Group result;
  Group::dbl(result, point);
  return result;
}

CanonicalScalar canonical_bytes(const Scalar& scalar) {
  CanonicalScalar bytes{};
  if (scalar.getLittleEndian(bytes.data(), bytes.size()) == 0) {
    throw std::runtime_error("failed to obtain canonical scalar");
  }
  return bytes;
}

std::size_t digit_at(
    const CanonicalScalar& scalar,
    std::size_t shift,
    std::size_t width,
    std::size_t exponent_bit_length) {
  std::size_t digit = 0;
  for (std::size_t bit = 0; bit < width; ++bit) {
    const std::size_t source_bit = shift + bit;
    if (source_bit >= exponent_bit_length) {
      break;
    }
    const std::size_t byte_index = source_bit / 8;
    const std::size_t bit_index = source_bit % 8;
    digit |= static_cast<std::size_t>(
                 (scalar[byte_index] >> bit_index) & std::uint8_t{1})
             << bit;
  }
  return digit;
}

bool fits_unsigned_bits(
    const CanonicalScalar& scalar,
    std::size_t exponent_bit_length) {
  const std::size_t complete_bytes = exponent_bit_length / 8;
  const std::size_t remaining_bits = exponent_bit_length % 8;
  std::size_t first_zero_byte = complete_bytes;
  if (remaining_bits != 0) {
    const std::uint8_t allowed =
        static_cast<std::uint8_t>(
            (std::uint16_t{1} << remaining_bits) - 1);
    if ((scalar[complete_bytes] &
         static_cast<std::uint8_t>(~allowed)) != 0) {
      return false;
    }
    first_zero_byte = complete_bytes + 1;
  }
  for (std::size_t i = first_zero_byte; i < scalar.size(); ++i) {
    if (scalar[i] != 0) return false;
  }
  return true;
}

void validate_dimensions(
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

}

void initialize() {
  static std::once_flag once;
  std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}

std::size_t scalar_bit_length() {
  initialize();
  return Scalar::getBitSize();
}

std::size_t number_of_windows(std::size_t window_width) {
  const std::size_t bit_length = scalar_bit_length();
  if (window_width == 0 || window_width > bit_length) {
    throw std::invalid_argument("w must be in [1, ceil(log2(p))]");
  }
  return (bit_length + window_width - 1) / window_width;
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
  validate_dimensions(bases, scalars);

  if (exponent_bit_length == 0 ||
      exponent_bit_length > scalar_bit_length()) {
    throw std::invalid_argument(
        "exponent bit length must fit the scalar field");
  }
  if (window_width == 0 || window_width > exponent_bit_length) {
    throw std::invalid_argument(
        "w must not exceed the exponent bit length");
  }
  const std::size_t windows =
      (exponent_bit_length + window_width - 1) / window_width;
  if (window_width >= std::numeric_limits<std::size_t>::digits) {
    throw std::invalid_argument("2^w does not fit in the bucket-index type");
  }
  const std::size_t bucket_count =
      (std::size_t{1} << window_width) - 1;

  std::vector<std::vector<CanonicalScalar>> canonical(
      scalars.size(), std::vector<CanonicalScalar>(bases.size()));
  for (std::size_t j = 0; j < scalars.size(); ++j) {
    for (std::size_t i = 0; i < bases.size(); ++i) {
      canonical[j][i] = canonical_bytes(scalars[j][i]);
      if (!fits_unsigned_bits(canonical[j][i], exponent_bit_length)) {
        throw std::invalid_argument(
            "scalar exceeds the declared exponent bit length");
      }
    }
  }



  std::vector<Group> buckets(bucket_count);
  std::vector<Group> results(scalars.size());

  for (std::size_t j = 0; j < scalars.size(); ++j) {
    Group accumulator;
    accumulator.clear();

    for (std::size_t ell = windows; ell-- > 0;) {
      if (ell + 1 < windows) {
        for (std::size_t doubling = 0; doubling < window_width; ++doubling) {
          accumulator = double_point(accumulator);
        }
      }

      for (auto& bucket : buckets) {
        bucket.clear();
      }

      const std::size_t shift = window_width * ell;
      for (std::size_t i = 0; i < bases.size(); ++i) {
        const std::size_t digit =
            digit_at(
                canonical[j][i], shift, window_width,
                exponent_bit_length);
        if (digit != 0) {
          buckets[digit - 1] = add(buckets[digit - 1], bases[i]);
        }
      }

      Group running_sum;
      running_sum.clear();
      for (std::size_t t = bucket_count; t-- > 0;) {
        running_sum = add(running_sum, buckets[t]);
        accumulator = add(accumulator, running_sum);
      }
    }
    results[j] = accumulator;
  }

  return results;
}

std::vector<Group> naive_multi_exponentiation(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  initialize();
  validate_dimensions(bases, scalars);

  std::vector<Group> results(scalars.size());
  for (std::size_t j = 0; j < scalars.size(); ++j) {
    results[j].clear();
    for (std::size_t i = 0; i < bases.size(); ++i) {
      Group term;
      Group::mul(term, bases[i], scalars[j][i]);
      results[j] = add(results[j], term);
    }
  }
  return results;
}

}

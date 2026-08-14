#include "internal/pippenger_core.hpp"

#include <cstdint>
#include <stdexcept>

namespace pippenger::internal {

Group add(const Group& left, const Group& right) {
  Group result;
  Group::add(result, left, right);
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

bool fits_unsigned_bits(
    const CanonicalScalar& scalar,
    std::size_t exponent_bit_length) {
  const std::size_t complete_bytes = exponent_bit_length / 8;
  const std::size_t remaining_bits = exponent_bit_length % 8;
  std::size_t first_zero_byte = complete_bytes;
  if (remaining_bits != 0) {
    const std::uint8_t allowed = static_cast<std::uint8_t>(
        (std::uint16_t{1} << remaining_bits) - 1);
    if ((scalar[complete_bytes] & static_cast<std::uint8_t>(~allowed)) != 0) {
      return false;
    }
    first_zero_byte = complete_bytes + 1;
  }
  for (std::size_t i = first_zero_byte; i < scalar.size(); ++i) {
    if (scalar[i] != 0) return false;
  }
  return true;
}

std::size_t window_digit(
    const CanonicalScalar& scalar,
    std::size_t shift,
    std::size_t width,
    std::size_t exponent_bit_length) {
  std::size_t digit = 0;
  for (std::size_t bit = 0; bit < width; ++bit) {
    const std::size_t source_bit = shift + bit;
    if (source_bit >= exponent_bit_length) break;
    const std::size_t byte_index = source_bit / 8;
    const std::size_t bit_index = source_bit % 8;
    digit |= static_cast<std::size_t>(
                 (scalar[byte_index] >> bit_index) & std::uint8_t{1})
             << bit;
  }
  return digit;
}

Group evaluate_pippenger_row(
    const std::vector<Group>& bases,
    const std::vector<CanonicalScalar>& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length,
    std::size_t window_count,
    std::size_t bucket_count) {
  std::vector<Group> buckets(bucket_count);
  Group accumulator;
  accumulator.clear();

  for (std::size_t window = window_count; window-- > 0;) {
    if (window + 1 < window_count) {
      for (std::size_t doubling = 0; doubling < window_width; ++doubling) {
        accumulator = double_point(accumulator);
      }
    }

    for (Group& bucket : buckets) bucket.clear();
    const std::size_t shift = window_width * window;
    for (std::size_t i = 0; i < bases.size(); ++i) {
      const std::size_t digit = window_digit(
          scalars[i], shift, window_width, exponent_bit_length);
      if (digit != 0) {
        buckets[digit - 1] = add(buckets[digit - 1], bases[i]);
      }
    }

    Group running_sum;
    running_sum.clear();
    for (std::size_t bucket = bucket_count; bucket-- > 0;) {
      running_sum = add(running_sum, buckets[bucket]);
      accumulator = add(accumulator, running_sum);
    }
  }
  return accumulator;
}

}  // namespace pippenger::internal

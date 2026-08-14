#include "internal/protocol_utils.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

namespace bp::internal {

Bytes u64be(std::uint64_t value) {
  Bytes out(8);
  for (std::size_t i = 0; i < out.size(); ++i)
    out[i] = static_cast<std::uint8_t>(value >> (56 - 8 * i));
  return out;
}

void append_raw(Bytes& out, std::span<const std::uint8_t> field) {
  out.insert(out.end(), field.begin(), field.end());
}

void frame(Bytes& out, std::span<const std::uint8_t> field) {
  append_raw(out, u64be(field.size()));
  append_raw(out, field);
}

void frame(Bytes& out, std::string_view field) {
  frame(out, {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

bool read_u64(std::span<const std::uint8_t> bytes, std::size_t offset,
              std::uint64_t& value) noexcept {
  if (offset > bytes.size() || bytes.size() - offset < 8) return false;
  value = 0;
  for (std::size_t i = 0; i < 8; ++i) value = (value << 8) | bytes[offset + i];
  return true;
}

bool power_of_two(std::size_t value) noexcept {
  return value != 0 && (value & (value - 1)) == 0;
}

std::size_t exact_log2(std::size_t value) {
  if (!power_of_two(value))
    throw std::invalid_argument("dimension must be a power of two");
  std::size_t result = 0;
  while (value > 1) {
    value >>= 1;
    ++result;
  }
  return result;
}

bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept {
  if (right > std::numeric_limits<std::size_t>::max() - left) return false;
  result = left + right;
  return true;
}

bool checked_mul(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept {
  if (left != 0 && right > std::numeric_limits<std::size_t>::max() / left)
    return false;
  result = left * right;
  return true;
}

Scalar fadd(const Scalar& left, const Scalar& right) {
  Scalar out;
  Scalar::add(out, left, right);
  return out;
}

Scalar fmul(const Scalar& left, const Scalar& right) {
  Scalar out;
  Scalar::mul(out, left, right);
  return out;
}

Scalar finv(const Scalar& value) {
  if (value.isZero()) throw std::invalid_argument("zero scalar");
  Scalar out;
  Scalar::inv(out, value);
  return out;
}

Scalar fneg(const Scalar& value) {
  Scalar out;
  Scalar::neg(out, value);
  return out;
}

Group gadd(const Group& left, const Group& right) {
  Group out;
  Group::add(out, left, right);
  return out;
}

Group gmul(const Group& point, const Scalar& scalar) {
  Group out;
  Group::mul(out, point, scalar);
  return out;
}

Group msm(std::span<const Group> points, std::span<const Scalar> scalars) {
  if (points.size() != scalars.size())
    throw std::invalid_argument("MSM length mismatch");
  Group out;
  out.clear();
  if (!points.empty()) {
    std::vector<Group> work(points.begin(), points.end());
    Group::mulVec(out, work.data(), scalars.data(), work.size());
  }
  return out;
}

Scalar inner_product(std::span<const Scalar> left,
                     std::span<const Scalar> right) {
  if (left.size() != right.size())
    throw std::invalid_argument("inner-product length mismatch");
  Scalar out;
  out.clear();
  for (std::size_t i = 0; i < left.size(); ++i)
    out = fadd(out, fmul(left[i], right[i]));
  return out;
}

bool valid_group(const Group& point, bool require_nonidentity) noexcept {
  return point.isValid() && point.isValidOrder() &&
         (!require_nonidentity || !point.isZero());
}

bool canonical_scalar(std::span<const std::uint8_t> bytes, Scalar& out) noexcept {
  try {
    Scalar candidate;
    if (bytes.size() != scalar_bytes() ||
        candidate.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
        serialize(candidate) != Bytes(bytes.begin(), bytes.end()))
      return false;
    out = candidate;
    return true;
  } catch (...) {
    return false;
  }
}

bool canonical_group(std::span<const std::uint8_t> bytes, Group& out,
                     bool require_nonidentity) noexcept {
  try {
    Group candidate;
    if (bytes.size() != group_bytes() ||
        candidate.deserialize(bytes.data(), bytes.size()) != bytes.size() ||
        !valid_group(candidate, require_nonidentity) ||
        serialize(candidate) != Bytes(bytes.begin(), bytes.end()))
      return false;
    out = candidate;
    return true;
  } catch (...) {
    return false;
  }
}

}  

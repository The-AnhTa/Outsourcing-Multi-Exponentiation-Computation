#include "internal/crypto.hpp"

#include <limits>
#include <vector>

namespace blsagg::internal {
namespace {

template<class Group>
bool canonical(const Group& value) {
  try {
    const Bytes encoded = serialize(value);
    Group decoded;
    return decoded.deserialize(encoded.data(), encoded.size()) ==
               encoded.size() &&
           decoded == value;
  } catch (...) {
    return false;
  }
}

}  // namespace

void append_raw(Bytes& out, std::span<const std::uint8_t> value) {
  out.insert(out.end(), value.begin(), value.end());
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

Fr fm(const Fr& left, const Fr& right) {
  Fr out;
  Fr::mul(out, left, right);
  return out;
}

Fr inv(const Fr& value) {
  if (value.isZero()) throw std::invalid_argument("zero challenge");
  Fr out;
  Fr::inv(out, value);
  return out;
}

G1 add(const G1& left, const G1& right) {
  G1 out;
  G1::add(out, left, right);
  return out;
}

G2 add(const G2& left, const G2& right) {
  G2 out;
  G2::add(out, left, right);
  return out;
}

G1 mul(const G1& point, const Fr& scalar) {
  G1 out;
  G1::mul(out, point, scalar);
  return out;
}

G2 mul(const G2& point, const Fr& scalar) {
  G2 out;
  G2::mul(out, point, scalar);
  return out;
}

GT gmul(const GT& left, const GT& right) {
  GT out;
  GT::mul(out, left, right);
  return out;
}

GT gpow(const GT& value, const Fr& exponent) {
  GT out;
  GT::pow(out, value, exponent);
  return out;
}

GT product(std::initializer_list<GT> values) {
  GT out;
  out.setOne();
  for (const auto& value : values) out = gmul(out, value);
  return out;
}

GT pair(const G1& left, const G2& right) {
  GT out;
  mcl::bn::pairing(out, left, right);
  return out;
}

G1 msm(std::span<const G1> points, std::span<const Fr> scalars) {
  if (points.size() != scalars.size())
    throw std::invalid_argument("G1 MSM length mismatch");
  G1 out;
  out.clear();
  if (!points.empty()) {
    std::vector<G1> work(points.begin(), points.end());
    G1::mulVec(out, work.data(), scalars.data(), work.size());
  }
  return out;
}

G2 msm(std::span<const G2> points, std::span<const Fr> scalars) {
  if (points.size() != scalars.size())
    throw std::invalid_argument("G2 MSM length mismatch");
  G2 out;
  out.clear();
  if (!points.empty()) {
    std::vector<G2> work(points.begin(), points.end());
    G2::mulVec(out, work.data(), scalars.data(), work.size());
  }
  return out;
}

bool valid(const G1& point, bool require_nonidentity) noexcept {
  return canonical(point) && point.isValid() && point.isValidOrder() &&
         (!require_nonidentity || !point.isZero());
}

bool valid(const G2& point, bool require_nonidentity) noexcept {
  return canonical(point) && point.isValid() && point.isValidOrder() &&
         (!require_nonidentity || !point.isZero());
}

bool valid(const GT& value) noexcept {
  return canonical(value) && mcl::bn::isValidGT(value);
}

}  // namespace blsagg::internal

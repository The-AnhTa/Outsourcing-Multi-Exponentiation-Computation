#pragma once

#include "blsagg/transcript.hpp"

#include <initializer_list>
#include <span>
#include <stdexcept>

namespace blsagg::internal {

template<class T>
Bytes encode(const T& value) {
  Bytes out(2048);
  const auto written = value.serialize(out.data(), out.size());
  if (written == 0) throw std::runtime_error("MCL serialization failed");
  out.resize(written);
  return out;
}

void append_raw(Bytes& out, std::span<const std::uint8_t> value);
bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept;
bool checked_mul(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept;

Fr fm(const Fr& left, const Fr& right);
Fr inv(const Fr& value);
G1 add(const G1& left, const G1& right);
G2 add(const G2& left, const G2& right);
G1 mul(const G1& point, const Fr& scalar);
G2 mul(const G2& point, const Fr& scalar);
GT gmul(const GT& left, const GT& right);
GT gpow(const GT& value, const Fr& exponent);
GT product(std::initializer_list<GT> values);
GT pair(const G1& left, const G2& right);
G1 msm(std::span<const G1> points, std::span<const Fr> scalars);
G2 msm(std::span<const G2> points, std::span<const Fr> scalars);

bool valid(const G1& point, bool require_nonidentity = false) noexcept;
bool valid(const G2& point, bool require_nonidentity = false) noexcept;
bool valid(const GT& value) noexcept;

}  // namespace blsagg::internal

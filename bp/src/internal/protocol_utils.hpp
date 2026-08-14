#pragma once

#include "bp/bp.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string_view>

namespace bp::internal {

Bytes u64be(std::uint64_t value);
void append_raw(Bytes& out, std::span<const std::uint8_t> field);
void frame(Bytes& out, std::span<const std::uint8_t> field);
void frame(Bytes& out, std::string_view field);
bool read_u64(std::span<const std::uint8_t> bytes, std::size_t offset,
              std::uint64_t& value) noexcept;

bool power_of_two(std::size_t value) noexcept;
std::size_t exact_log2(std::size_t value);
bool checked_add(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept;
bool checked_mul(std::size_t left, std::size_t right,
                 std::size_t& result) noexcept;

Scalar fadd(const Scalar& left, const Scalar& right);
Scalar fmul(const Scalar& left, const Scalar& right);
Scalar finv(const Scalar& value);
Scalar fneg(const Scalar& value);
Group gadd(const Group& left, const Group& right);
Group gmul(const Group& point, const Scalar& scalar);
Group msm(std::span<const Group> points, std::span<const Scalar> scalars);
Scalar inner_product(std::span<const Scalar> left,
                     std::span<const Scalar> right);

bool valid_group(const Group& point, bool require_nonidentity = false) noexcept;
bool canonical_scalar(std::span<const std::uint8_t> bytes, Scalar& out) noexcept;
bool canonical_group(std::span<const std::uint8_t> bytes, Group& out,
                     bool require_nonidentity = false) noexcept;

}  // namespace bp::internal

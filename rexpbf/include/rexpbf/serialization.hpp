#pragma once
#include "rexpbf/types.hpp"
#include <span>
#include <vector>

namespace rexpbf {
std::array<std::uint8_t, 8> encode_u64_be(std::uint64_t value);
void append_frame(std::vector<std::uint8_t>& out, std::span<const std::uint8_t> bytes);
std::vector<std::uint8_t> serialize_g1(const G1& point);
std::vector<std::uint8_t> serialize_g2(const G2& point);
std::vector<std::uint8_t> serialize_gt(const GT& value);
std::vector<std::uint8_t> serialize_fr(const Fr& value);
std::vector<std::uint8_t> serialize_crs(const CRS& crs);
std::vector<std::uint8_t> serialize_crs_wire(const CRS& crs);
std::vector<std::uint8_t> serialize_statement_wire(const CRS& crs, const Statement& statement);
G1 deserialize_g1(std::span<const std::uint8_t> bytes);
G2 deserialize_g2(std::span<const std::uint8_t> bytes);
GT deserialize_gt(std::span<const std::uint8_t> bytes);
CRS deserialize_crs_wire(std::span<const std::uint8_t> bytes);
Statement deserialize_statement_wire(std::span<const std::uint8_t> bytes, const CRS& crs);
Digest32 sha256(std::span<const std::uint8_t> bytes);
Digest32 compute_crs_digest(const CRS& crs);
Digest32 compute_statement_digest(const CRS& crs, const Statement& statement);
}

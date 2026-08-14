#pragma once
#include "vpip_bf/proof.hpp"
#include "vpip_bf/setup.hpp"
#include <span>

namespace vpip_bf {
inline constexpr std::size_t kMaxSerializedDimension = 20;
struct CanonicalElementSizes {std::size_t fr_bytes{},g1_bytes{},g2_bytes{},gt_bytes{};};
enum class DecodeError {None,WrongMagic,UnsupportedVersion,WrongCurve,InvalidDimension,SizeOverflow,WrongLength,Truncated,TrailingBytes,InvalidFr,InvalidG1,InvalidG2,InvalidGT,NonCanonicalEncoding,InvalidSubgroup,IdentityNotAllowed,InvalidCrs,CrsDigestMismatch};

CanonicalElementSizes canonical_element_sizes();
Bytes serialize_fr(const Fr&); Bytes serialize_g1(const G1&); Bytes serialize_g2(const G2&); Bytes serialize_gt(const GT&);
bool deserialize_fr(std::span<const std::uint8_t>,Fr&); bool deserialize_g1(std::span<const std::uint8_t>,G1&); bool deserialize_g2(std::span<const std::uint8_t>,G2&); bool deserialize_gt(std::span<const std::uint8_t>,GT&);
void append_u16_be(Bytes&,std::uint16_t); void append_u32_be(Bytes&,std::uint32_t); void append_u64_be(Bytes&,std::uint64_t);

std::size_t proof_mathematical_payload_bytes(std::size_t,const CanonicalElementSizes&);
std::size_t proof_wire_bytes(std::size_t,const CanonicalElementSizes&);
std::size_t crs_wire_bytes(std::size_t,const CanonicalElementSizes&);
std::size_t precomputation_wire_bytes(std::size_t,const CanonicalElementSizes&);
std::size_t statement_wire_bytes(std::size_t,const CanonicalElementSizes&);

Bytes serialize_crs(const VpipBfCRS&);
bool deserialize_crs(std::span<const std::uint8_t>,VpipBfCRS&,DecodeError* =nullptr);
Bytes serialize_precomputation(const VpipBfCRS&,const VpipBfPrecomputation&);
bool deserialize_precomputation(std::span<const std::uint8_t>,const VpipBfCRS&,VpipBfPrecomputation&,DecodeError* =nullptr);
bool validate_precomputation(const VpipBfCRS&,const VpipBfPrecomputation&);
Bytes serialize_statement(const VpipBfCRS&,const VpipBfStatement&);
bool deserialize_statement(std::span<const std::uint8_t>,const VpipBfCRS&,VpipBfStatement&,DecodeError* =nullptr);
Bytes serialize_proof(std::size_t,const VpipBfProof&);
bool deserialize_proof(std::span<const std::uint8_t>,const VpipBfCRS&,VpipBfProof&,DecodeError* =nullptr);
}

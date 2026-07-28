#pragma once
#include "vme_ibf/types.hpp"
#include <span>
#include <string>

namespace vme_ibf {
void initialize();
Bytes encode_u64_be(std::uint64_t value);
void append_frame(Bytes& out, std::span<const std::uint8_t> field);
void append_frame(Bytes& out, std::string_view field);
Digest sha256(std::span<const std::uint8_t> input);
Bytes serialize(const Fr&); Bytes serialize(const G1&); Bytes serialize(const G2&); Bytes serialize(const GT&);
bool valid_g1(const G1&, bool nonidentity = false); bool valid_g2(const G2&, bool nonidentity = false);
GT pairing_product(std::span<const G1>, std::span<const G2>);
G1 g1_multiexp_reference(std::span<const G1>, std::span<const Fr>);
G1 g1_multiexp(std::span<const G1>, std::span<const Fr>);
G2 g2_multiexp_reference(std::span<const G2>, std::span<const Fr>);
G2 g2_multiexp(std::span<const G2>, std::span<const Fr>);
G2 g2_multiexp_protocol(std::span<const G2>, std::span<const Fr>);
std::vector<G1> component_mul(std::span<const G1>, std::span<const G1>);
std::vector<G2> component_mul(std::span<const G2>, std::span<const G2>);
std::vector<G1> component_pow(std::span<const G1>, std::span<const Fr>);
std::vector<G2> component_pow(std::span<const G2>, std::span<const Fr>);
Fr inner_product(std::span<const Fr>, std::span<const Fr>);
Fr inverse_nonzero(const Fr&);
G1 g1_add(const G1&, const G1&); G2 g2_add(const G2&, const G2&);
G1 g1_pow(const G1&, const Fr&); G2 g2_pow(const G2&, const Fr&);
GT gt_mul(const GT&, const GT&); GT gt_pow(const GT&, const Fr&);
std::string hex(std::span<const std::uint8_t>);
std::vector<Fr> tensor_vector(std::span<const Fr> r);
}

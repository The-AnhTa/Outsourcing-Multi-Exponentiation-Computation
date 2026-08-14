#pragma once

#include "internal/constants.hpp"
#include "pinkas/pinkas.hpp"

#include <cstdint>
#include <span>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace pinkas::internal {

void append_raw(Bytes& output, std::span<const std::uint8_t> input);
Bytes encode_u64(std::uint64_t value);
void append_frame(Bytes& output, std::span<const std::uint8_t> field);
void append_frame(Bytes& output, std::string_view field);

template <class Element>
Bytes serialize_element(const Element& element) {
  Bytes bytes(512);
  const std::size_t written = element.serialize(bytes.data(), bytes.size());
  if (written == 0) throw std::runtime_error("mcl serialization failed");
  bytes.resize(written);
  return bytes;
}

template <class Element>
void append_element(Bytes& output, const Element& element) {
  append_frame(output, serialize_element(element));
}

Digest sha256(std::span<const std::uint8_t> input);
CanonicalScalar canonical_scalar_bytes(const Scalar& scalar);
bool scalar_bit(const CanonicalScalar& scalar, std::size_t bit);

Digest build_transcript(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof);

std::vector<std::vector<Scalar>> derive_challenges(
    const PublicParameters& parameters,
    const Digest& transcript,
    std::size_t instances);

}  

#include "internal/transcript.hpp"

#include <mcl/fp.hpp>

#include <cstdint>
#include <limits>
#include <stdexcept>

namespace pinkas::internal {

void append_raw(Bytes& output, std::span<const std::uint8_t> input) {
  output.insert(output.end(), input.begin(), input.end());
}

Bytes encode_u64(std::uint64_t value) {
  Bytes output(8);
  for (int shift = 56, i = 0; shift >= 0; shift -= 8, ++i) {
    output[static_cast<std::size_t>(i)] =
        static_cast<std::uint8_t>(value >> shift);
  }
  return output;
}

void append_frame(Bytes& output, std::span<const std::uint8_t> field) {
  append_raw(output, encode_u64(field.size()));
  append_raw(output, field);
}

void append_frame(Bytes& output, std::string_view field) {
  append_frame(
      output,
      {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

Digest sha256(std::span<const std::uint8_t> input) {
  if (input.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::length_error("SHA-256 input is too large");
  }
  Digest digest{};
  const std::uint32_t written = mcl::fp::sha256(
      digest.data(), static_cast<std::uint32_t>(digest.size()), input.data(),
      static_cast<std::uint32_t>(input.size()));
  if (written != digest.size()) throw std::runtime_error("SHA-256 failed");
  return digest;
}

CanonicalScalar canonical_scalar_bytes(const Scalar& scalar) {
  CanonicalScalar bytes{};
  if (scalar.getLittleEndian(bytes.data(), bytes.size()) == 0) {
    throw std::runtime_error("failed to obtain canonical scalar");
  }
  return bytes;
}

bool scalar_bit(const CanonicalScalar& scalar, std::size_t bit) {
  return ((scalar[bit / 8] >> (bit % 8)) & std::uint8_t{1}) != 0;
}

Digest build_transcript(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  Bytes transcript;
  append_frame(transcript, parameters.domain);
  append_frame(transcript, parameters.group_identifier);
  append_frame(transcript, encode_u64(bases.size()));
  append_frame(transcript, encode_u64(scalars.size()));
  append_frame(transcript, encode_u64(parameters.scalar_bit_length));
  append_frame(transcript, encode_u64(parameters.challenge_bit_length));
  for (const Group& base : bases) append_element(transcript, base);
  for (const auto& row : scalars) {
    for (const Scalar& scalar : row) append_element(transcript, scalar);
  }
  for (const Group& output : outputs) append_element(transcript, output);
  for (const auto& row : proof.W) {
    for (const Group& point : row) append_element(transcript, point);
  }
  return sha256(transcript);
}

namespace {

Scalar hash_to_challenge(
    const PublicParameters& parameters,
    std::size_t instance,
    std::size_t bit,
    const Digest& transcript) {
  Bytes input;
  append_frame(input, parameters.domain);
  append_frame(input, "r");
  append_frame(input, encode_u64(instance));
  append_frame(input, encode_u64(bit));
  append_frame(input, transcript);
  const Digest candidate = sha256(input);
  Scalar result;
  result.setBigEndianMod(candidate.data(), kChallengeBytes);
  return result;
}

}  // namespace

std::vector<std::vector<Scalar>> derive_challenges(
    const PublicParameters& parameters,
    const Digest& transcript,
    std::size_t instances) {
  std::vector<std::vector<Scalar>> challenges(
      instances, std::vector<Scalar>(parameters.scalar_bit_length));
  for (std::size_t instance = 0; instance < instances; ++instance) {
    for (std::size_t bit = 0; bit < parameters.scalar_bit_length; ++bit) {
      challenges[instance][bit] =
          hash_to_challenge(parameters, instance, bit, transcript);
    }
  }
  return challenges;
}

}  // namespace pinkas::internal

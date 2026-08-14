#include "pinkas/pinkas.hpp"

#include "internal/constants.hpp"
#include "internal/transcript.hpp"
#include "internal/validation.hpp"

#include <cstdint>
#include <limits>
#include <new>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace pinkas {
namespace {

class Reader {
 public:
  explicit Reader(std::span<const std::uint8_t> input) : input_(input) {}

  bool read_u64(std::uint64_t& value) {
    if (input_.size() - position_ < 8) return false;
    value = 0;
    for (std::size_t i = 0; i < 8; ++i) {
      value = (value << 8) | input_[position_++];
    }
    return true;
  }

  bool read_frame(std::span<const std::uint8_t>& field) {
    std::uint64_t length = 0;
    if (!read_u64(length) || length > input_.size() - position_) return false;
    field = input_.subspan(position_, static_cast<std::size_t>(length));
    position_ += static_cast<std::size_t>(length);
    return true;
  }

  std::size_t remaining() const { return input_.size() - position_; }
  bool finished() const { return remaining() == 0; }

 private:
  std::span<const std::uint8_t> input_;
  std::size_t position_ = 0;
};

}  

Bytes serialize_proof(
    const PublicParameters& parameters,
    const PinkasProof& proof) {
  if (!internal::valid_proof(parameters, proof)) {
    throw std::invalid_argument("invalid Pinkas proof serialization inputs");
  }
  Bytes output;
  internal::append_frame(output, internal::kProtocolVersion);
  internal::append_raw(output, internal::encode_u64(proof.W.size()));
  internal::append_raw(
      output, internal::encode_u64(parameters.scalar_bit_length));
  for (const auto& row : proof.W) {
    for (const Group& point : row) internal::append_element(output, point);
  }
  return output;
}

DecodeError deserialize_proof_detailed(
    const PublicParameters& parameters,
    std::span<const std::uint8_t> encoded,
  PinkasProof& proof) noexcept {
  proof = {};
  try {
    if (!internal::valid_parameters(parameters)) {
      return DecodeError::InvalidParameters;
    }
    if (encoded.empty()) return DecodeError::EmptyInput;

    Reader reader(encoded);
    std::span<const std::uint8_t> version;
    if (!reader.read_frame(version)) return DecodeError::Truncated;
    if (std::string_view(
            reinterpret_cast<const char*>(version.data()), version.size()) !=
        internal::kProtocolVersion) {
      return DecodeError::InvalidVersion;
    }

    std::uint64_t encoded_instances = 0;
    std::uint64_t encoded_bits = 0;
    if (!reader.read_u64(encoded_instances) ||
        !reader.read_u64(encoded_bits)) {
      return DecodeError::Truncated;
    }
    if (encoded_instances == 0 ||
        encoded_bits != parameters.scalar_bit_length) {
      return DecodeError::InvalidDimensions;
    }
    if (encoded_instances > std::numeric_limits<std::size_t>::max()) {
      return DecodeError::AllocationLimit;
    }

    const std::size_t instances =
        static_cast<std::size_t>(encoded_instances);
    const std::size_t bits = parameters.scalar_bit_length;
    if (instances > std::vector<std::vector<Group>>{}.max_size() ||
        bits > std::vector<Group>{}.max_size() ||
        instances > std::numeric_limits<std::size_t>::max() / bits) {
      return DecodeError::AllocationLimit;
    }
    const std::size_t point_count = instances * bits;
    constexpr std::size_t kMinimumFramedPointBytes = 9;
    if (point_count > reader.remaining() / kMinimumFramedPointBytes) {
      return DecodeError::Truncated;
    }

    PinkasProof parsed;
    parsed.W.resize(instances, std::vector<Group>(bits));
    for (auto& row : parsed.W) {
      for (Group& point : row) {
        std::span<const std::uint8_t> field;
        if (!reader.read_frame(field)) return DecodeError::Truncated;
        if (field.empty() ||
            point.deserialize(field.data(), field.size()) != field.size() ||
            !internal::valid_point(point)) {
          return DecodeError::InvalidPoint;
        }
        if (internal::serialize_element(point) !=
            Bytes(field.begin(), field.end())) {
          return DecodeError::NonCanonicalPoint;
        }
      }
    }
    if (!reader.finished()) return DecodeError::TrailingData;
    proof = std::move(parsed);
    return DecodeError::None;
  } catch (const std::bad_alloc&) {
    proof = {};
    return DecodeError::AllocationLimit;
  } catch (...) {
    proof = {};
    return DecodeError::InternalError;
  }
}

bool deserialize_proof(
    const PublicParameters& parameters,
    std::span<const std::uint8_t> encoded,
    PinkasProof& proof) {
  return deserialize_proof_detailed(parameters, encoded, proof) ==
         DecodeError::None;
}

bool verify_serialized(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    std::span<const std::uint8_t> encoded_proof) {
  PinkasProof proof;
  if (deserialize_proof_detailed(parameters, encoded_proof, proof) !=
      DecodeError::None) {
    return false;
  }
  ValidatedInputs validated;
  return prepare_validated_inputs(
             parameters, bases, scalars, outputs, proof, validated) &&
         verify_online_prevalidated(validated);
}

}  

#include "pinkas.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace pinkas {
namespace {

constexpr std::string_view kProtocolVersion = "PinkasFS128BN-v2";
constexpr std::string_view kGroupIdentifier = "BN254/G2/mcl-v3.00";
constexpr std::size_t kScalarBytes = 32;
constexpr std::size_t kChallengeBitLength = 128;
constexpr std::size_t kChallengeBytes = kChallengeBitLength / 8;
using Digest = std::array<std::uint8_t, 32>;
using CanonicalScalar = std::array<std::uint8_t, kScalarBytes>;

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

template <class T>
Bytes serialize_element(const T& element) {
  Bytes bytes(512);
  const std::size_t written = element.serialize(bytes.data(), bytes.size());
  if (written == 0) {
    throw std::runtime_error("mcl serialization failed");
  }
  bytes.resize(written);
  return bytes;
}

template <class T>
void append_element(Bytes& output, const T& element) {
  append_frame(output, serialize_element(element));
}

Digest sha256(std::span<const std::uint8_t> input) {
  if (input.size() > std::numeric_limits<std::uint32_t>::max()) {
    throw std::length_error("SHA-256 input is too large");
  }
  Digest digest{};
  const std::uint32_t written = mcl::fp::sha256(
      digest.data(),
      static_cast<std::uint32_t>(digest.size()),
      input.data(),
      static_cast<std::uint32_t>(input.size()));
  if (written != digest.size()) {
    throw std::runtime_error("SHA-256 failed");
  }
  return digest;
}

CanonicalScalar canonical_bytes(const Scalar& scalar) {
  CanonicalScalar bytes{};
  if (scalar.getLittleEndian(bytes.data(), bytes.size()) == 0) {
    throw std::runtime_error("failed to obtain canonical scalar");
  }
  return bytes;
}

bool scalar_bit(const CanonicalScalar& scalar, std::size_t bit) {
  return ((scalar[bit / 8] >> (bit % 8)) & std::uint8_t{1}) != 0;
}

Group add(const Group& a, const Group& b) {
  Group result;
  Group::add(result, a, b);
  return result;
}

Group double_point(const Group& point) {
  Group result;
  Group::dbl(result, point);
  return result;
}

Scalar scalar_add(const Scalar& a, const Scalar& b) {
  Scalar result;
  Scalar::add(result, a, b);
  return result;
}

bool valid_point(const Group& point) {
  return point.isValid() && point.isValidOrder();
}

bool valid_dimensions(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  if (parameters.group_identifier != kGroupIdentifier ||
      parameters.domain != kProtocolVersion ||
      parameters.scalar_bit_length != Scalar::getBitSize() ||
      parameters.challenge_bit_length != kChallengeBitLength ||
      parameters.msm_window_width == 0 ||
      bases.empty() ||
      scalars.empty() ||
      outputs.size() != scalars.size() ||
      proof.W.size() != scalars.size()) {
    return false;
  }
  for (std::size_t ell = 0; ell < scalars.size(); ++ell) {
    if (scalars[ell].size() != bases.size() ||
        proof.W[ell].size() != parameters.scalar_bit_length) {
      return false;
    }
  }
  return true;
}

bool validate_points(
    const std::vector<Group>& bases,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  for (const Group& base : bases) {
    if (!valid_point(base)) {
      return false;
    }
  }
  for (std::size_t ell = 0; ell < outputs.size(); ++ell) {
    if (!valid_point(outputs[ell])) {
      return false;
    }
    for (const Group& point : proof.W[ell]) {
      if (!valid_point(point)) {
        return false;
      }
    }
  }
  return true;
}

bool validate_statement_points(
    const std::vector<Group>& bases,
    const std::vector<Group>& outputs) {
  for (const Group& base : bases) {
    if (!valid_point(base)) {
      return false;
    }
  }
  for (const Group& output : outputs) {
    if (!valid_point(output)) {
      return false;
    }
  }
  return true;
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

  for (const Group& base : bases) {
    append_element(transcript, base);
  }
  for (const auto& row : scalars) {
    for (const Scalar& scalar : row) {
      append_element(transcript, scalar);
    }
  }
  for (const Group& output : outputs) {
    append_element(transcript, output);
  }
  for (const auto& row : proof.W) {
    for (const Group& point : row) {
      append_element(transcript, point);
    }
  }
  return sha256(transcript);
}

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

std::vector<std::vector<Scalar>> derive_challenges(
    const PublicParameters& parameters,
    const Digest& transcript,
    std::size_t instances) {
  std::vector<std::vector<Scalar>> challenges(
      instances, std::vector<Scalar>(parameters.scalar_bit_length));
  for (std::size_t ell = 0; ell < instances; ++ell) {
    for (std::size_t b = 0; b < parameters.scalar_bit_length; ++b) {
      challenges[ell][b] =
          hash_to_challenge(parameters, ell, b, transcript);
    }
  }
  return challenges;
}

Group horner_reconstruct(const std::vector<Group>& row) {
  Group result = row.back();
  for (std::size_t b = row.size() - 1; b-- > 0;) {
    result = double_point(result);
    result = add(result, row[b]);
  }
  return result;
}

Group pippenger_msm(
    const std::vector<Group>& points,
    const std::vector<Scalar>& scalars,
    std::size_t window_width,
    std::size_t exponent_bit_length) {
  ScalarMatrix one_instance(1, scalars);
  return pippenger::multi_pippenger_bounded(
             points, one_instance, window_width, exponent_bit_length)
      .front();
}

std::size_t ceil_log2(std::size_t value) {
  if (value == 0) throw std::invalid_argument("log2 of zero");
  std::size_t result = 0;
  for (--value; value != 0; value >>= 1) ++result;
  return result;
}

bool verify_arithmetic(const ValidatedInputs& validated) {
  const auto& parameters = *validated.parameters;
  const auto& bases = *validated.bases;
  const auto& scalars = *validated.scalars;
  const auto& outputs = *validated.outputs;
  const auto& proof = *validated.proof;
  const std::size_t k = scalars.size();
  const std::size_t m = parameters.scalar_bit_length;

  const Digest transcript =
      build_transcript(parameters, bases, scalars, outputs, proof);
  const auto challenges =
      derive_challenges(parameters, transcript, k);
  const std::size_t aggregate_exponent_bits =
      parameters.challenge_bit_length + ceil_log2(m);

  for (std::size_t ell = 0; ell < k; ++ell) {
    const Group w_prime = pippenger_msm(
        proof.W[ell], challenges[ell], parameters.msm_window_width,
        parameters.challenge_bit_length);

    std::vector<Scalar> exponents(bases.size());
    for (Scalar& exponent : exponents) exponent.clear();
    for (std::size_t i = 0; i < bases.size(); ++i) {
      const CanonicalScalar bytes = canonical_bytes(scalars[ell][i]);
      Scalar instance_sum;
      instance_sum.clear();
      for (std::size_t b = 0; b < m; ++b) {
        if (scalar_bit(bytes, b)) {
          instance_sum =
              scalar_add(instance_sum, challenges[ell][b]);
        }
      }
      exponents[i] = instance_sum;
    }
    const Group w_double_prime = pippenger_msm(
        bases, exponents, parameters.msm_window_width,
        aggregate_exponent_bits);
    if (w_prime != w_double_prime) return false;
    if (outputs[ell] != horner_reconstruct(proof.W[ell])) {
      return false;
    }
  }
  return true;
}

class Reader {
 public:
  explicit Reader(std::span<const std::uint8_t> input) : input_(input) {}

  bool read_u64(std::uint64_t& value) {
    if (position_ + 8 > input_.size()) {
      return false;
    }
    value = 0;
    for (std::size_t i = 0; i < 8; ++i) {
      value = (value << 8) | input_[position_++];
    }
    return true;
  }

  bool read_frame(std::span<const std::uint8_t>& field) {
    std::uint64_t length = 0;
    if (!read_u64(length) ||
        length > input_.size() - position_) {
      return false;
    }
    field = input_.subspan(position_, static_cast<std::size_t>(length));
    position_ += static_cast<std::size_t>(length);
    return true;
  }

  bool finished() const { return position_ == input_.size(); }

 private:
  std::span<const std::uint8_t> input_;
  std::size_t position_ = 0;
};

}

PublicParameters setup(std::size_t msm_window_width) {
  pippenger::initialize();

  static_cast<void>(pippenger::number_of_windows(msm_window_width));
  return {
      std::string(kGroupIdentifier),
      std::string(kProtocolVersion),
      Scalar::getBitSize(),
      kChallengeBitLength,
      msm_window_width};
}

ProverResult prove(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  if (bases.empty() || scalars.empty() ||
      parameters.group_identifier != kGroupIdentifier ||
      parameters.domain != kProtocolVersion ||
      parameters.scalar_bit_length != Scalar::getBitSize() ||
      parameters.challenge_bit_length != kChallengeBitLength) {
    throw std::invalid_argument("invalid Pinkas prover inputs");
  }
  for (const Group& base : bases) {
    if (!valid_point(base)) {
      throw std::invalid_argument("invalid base point");
    }
  }
  for (const auto& row : scalars) {
    if (row.size() != bases.size()) {
      throw std::invalid_argument("every scalar row must contain n elements");
    }
  }

  const std::size_t k = scalars.size();
  const std::size_t m = parameters.scalar_bit_length;
  PinkasProof proof;
  proof.W.resize(k, std::vector<Group>(m));
  for (auto& row : proof.W) {
    for (Group& point : row) {
      point.clear();
    }
  }
  std::vector<Group> outputs(k);

  for (std::size_t ell = 0; ell < k; ++ell) {
    for (std::size_t i = 0; i < bases.size(); ++i) {
      const CanonicalScalar bytes = canonical_bytes(scalars[ell][i]);
      for (std::size_t b = 0; b < m; ++b) {
        if (scalar_bit(bytes, b)) {
          proof.W[ell][b] = add(proof.W[ell][b], bases[i]);
        }
      }
    }
    outputs[ell] = horner_reconstruct(proof.W[ell]);
  }
  return {std::move(outputs), std::move(proof)};
}

bool prepare_validated_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof,
    ValidatedInputs& validated) {
  validated = {};
  try {
    if (!valid_dimensions(parameters, bases, scalars, outputs, proof) ||
        !validate_points(bases, outputs, proof)) {
      return false;
    }
    validated = {&parameters, &bases, &scalars, &outputs, &proof};
    return true;
  } catch (...) {
    return false;
  }
}

bool verify_online_prevalidated(const ValidatedInputs& validated) {
  if (!validated.parameters || !validated.bases || !validated.scalars ||
      !validated.outputs || !validated.proof) {
    return false;
  }
  try {
    return verify_arithmetic(validated);
  } catch (...) {
    return false;
  }
}

bool verify(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  ValidatedInputs validated;
  return prepare_validated_inputs(
             parameters, bases, scalars, outputs, proof, validated) &&
         verify_online_prevalidated(validated);
}

Bytes serialize_proof(
    const PublicParameters& parameters,
    const PinkasProof& proof) {
  Bytes output;
  append_frame(output, kProtocolVersion);
  append_raw(output, encode_u64(proof.W.size()));
  append_raw(output, encode_u64(parameters.scalar_bit_length));
  for (const auto& row : proof.W) {
    for (const Group& point : row) {
      append_element(output, point);
    }
  }
  return output;
}

bool deserialize_proof(
    const PublicParameters& parameters,
    std::span<const std::uint8_t> encoded,
    PinkasProof& proof) {
  proof = {};
  try {
    Reader reader(encoded);
    std::span<const std::uint8_t> version;
    std::uint64_t k = 0;
    std::uint64_t m = 0;
    if (!reader.read_frame(version) ||
        std::string_view(
            reinterpret_cast<const char*>(version.data()), version.size()) !=
            kProtocolVersion ||
        !reader.read_u64(k) ||
        !reader.read_u64(m) ||
        k == 0 ||
        m != parameters.scalar_bit_length ||
        m == 0 ||
        k > encoded.size() / static_cast<std::size_t>(m) / 9) {
      return false;
    }

    PinkasProof parsed;
    parsed.W.resize(
        static_cast<std::size_t>(k),
        std::vector<Group>(parameters.scalar_bit_length));
    for (auto& row : parsed.W) {
      for (Group& point : row) {
        std::span<const std::uint8_t> field;
        if (!reader.read_frame(field) ||
            point.deserialize(field.data(), field.size()) != field.size() ||
            serialize_element(point) != Bytes(field.begin(), field.end()) ||
            !valid_point(point)) {
          return false;
        }
      }
    }
    if (!reader.finished()) {
      return false;
    }
    proof = std::move(parsed);
    return true;
  } catch (...) {
    proof = {};
    return false;
  }
}

bool verify_serialized(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    std::span<const std::uint8_t> encoded_proof) {
  PinkasProof proof;
  if (!deserialize_proof(parameters, encoded_proof, proof) ||
      !valid_dimensions(parameters, bases, scalars, outputs, proof) ||
      !validate_statement_points(bases, outputs)) {
    return false;
  }
  const ValidatedInputs validated{
      &parameters, &bases, &scalars, &outputs, &proof};
  return verify_online_prevalidated(validated);
}

}

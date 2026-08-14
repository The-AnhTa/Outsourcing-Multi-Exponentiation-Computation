#include "pinkas/pinkas.hpp"
#include "pinkas.hpp"

#include "internal/transcript.hpp"

#include <mcl/fp.hpp>

#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <sstream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using pippenger::Group;
using pippenger::Scalar;
using pippenger::ScalarMatrix;

void check(bool condition, std::string_view message) {
  if (!condition) throw std::runtime_error(std::string(message));
}

template <class Function>
void check_invalid_argument(Function&& function, std::string_view message) {
  try {
    function();
  } catch (const std::invalid_argument&) {
    return;
  }
  throw std::runtime_error(std::string(message));
}

Scalar scalar(std::uint64_t value) {
  Scalar result;
  result.setStr(std::to_string(value));
  return result;
}

Group generator() {
  Group result;
  constexpr char domain[] = "pippenger-characterization-generator";
  mcl::bn::hashAndMapToG2(result, domain, sizeof(domain) - 1);
  return result;
}

std::vector<Group> deterministic_bases(std::size_t count) {
  const Group g = generator();
  std::vector<Group> bases(count);
  for (std::size_t i = 0; i < count; ++i) {
    Group::mul(bases[i], g, scalar(static_cast<std::uint64_t>(i + 1)));
  }
  return bases;
}

ScalarMatrix deterministic_scalars(std::size_t rows, std::size_t columns) {
  ScalarMatrix values(rows, std::vector<Scalar>(columns));
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      values[row][column] = scalar(static_cast<std::uint64_t>(
          (row + 1) * 17 + column * column + 3 * column));
    }
  }
  return values;
}

void test_pippenger_equivalence() {
  const auto bases = deterministic_bases(7);
  auto scalars = deterministic_scalars(3, bases.size());
  scalars[0][0].clear();
  const auto expected =
      pippenger::naive_multi_exponentiation(bases, scalars);
  constexpr std::array<std::size_t, 4> widths{1, 2, 4, 8};
  for (const std::size_t width : widths) {
    check(pippenger::multi_pippenger(bases, scalars, width) == expected,
          "Pippenger result differs from the naive result");
  }

  ScalarMatrix bounded{{scalar(0), scalar(1), scalar(31), scalar(255),
                        scalar(256), scalar(511), scalar(17)}};
  check(pippenger::multi_pippenger_bounded(bases, bounded, 3, 9) ==
            pippenger::naive_multi_exponentiation(bases, bounded),
        "bounded Pippenger result differs from the naive result");
  bounded[0][0] = scalar(512);
  check_invalid_argument(
      [&] { pippenger::multi_pippenger_bounded(bases, bounded, 3, 9); },
      "an out-of-range bounded scalar was accepted");
}

void test_pippenger_rejections() {
  const auto bases = deterministic_bases(2);
  const ScalarMatrix scalars{{scalar(1), scalar(2)}};
  check_invalid_argument(
      [&] { pippenger::number_of_windows(0); },
      "zero-width windows were accepted");
  check(pippenger::maximum_window_width() >= 8,
        "the supported window limit excludes the benchmark width");
  check_invalid_argument(
      [&] {
        pippenger::number_of_windows(
            pippenger::maximum_window_width() + 1);
      },
      "an unsafe window width was accepted");
  check_invalid_argument(
      [&] { pippenger::multi_pippenger({}, scalars, 1); },
      "empty bases were accepted");
  check_invalid_argument(
      [&] { pippenger::multi_pippenger(bases, {}, 1); },
      "empty scalar rows were accepted");
  check_invalid_argument(
      [&] {
        pippenger::multi_pippenger(
            bases, ScalarMatrix{{scalar(1)}}, 1);
      },
      "a ragged scalar matrix was accepted");
  check_invalid_argument(
      [&] { pippenger::multi_pippenger_bounded(bases, scalars, 1, 0); },
      "a zero exponent bound was accepted");
  check_invalid_argument(
      [&] {
        pippenger::multi_pippenger_bounded(
            bases, scalars, pippenger::maximum_window_width() + 1,
            pippenger::scalar_bit_length());
      },
      "bounded MSM accepted an unsafe window width");
}

void test_pinkas_round_trip() {
  const auto parameters = pinkas::setup(4);
  const auto bases = deterministic_bases(5);
  const auto scalars = deterministic_scalars(2, bases.size());
  const auto result = pinkas::prove(parameters, bases, scalars);

  const auto transcript = pinkas::internal::build_transcript(
      parameters, bases, scalars, result.Y, result.proof);
  std::ostringstream transcript_hex;
  transcript_hex << std::hex;
  for (const std::uint8_t byte : transcript) {
    transcript_hex.width(2);
    transcript_hex.fill('0');
    transcript_hex << static_cast<unsigned int>(byte);
  }
  constexpr std::string_view expected_transcript =
      "a7d8e98a580b36fc1efa09a58f758ff658e7988561c71aeb29d64e07d47e8949";
  check(transcript_hex.str() == expected_transcript,
        std::string("unexpected transcript digest: ") +
            transcript_hex.str());

  check(result.Y == pippenger::naive_multi_exponentiation(bases, scalars),
        "Pinkas outputs differ from direct MSM outputs");
  check(pinkas::verify(
            parameters, bases, scalars, result.Y, result.proof),
        "a valid Pinkas proof was rejected");

  pinkas::ValidatedInputs validated;
  check(pinkas::prepare_validated_inputs(
            parameters, bases, scalars, result.Y, result.proof, validated),
        "valid Pinkas inputs did not validate");
  check(pinkas::verify_online_prevalidated(validated),
        "prevalidated Pinkas verification failed");
  pinkas::ValidatedInputs empty_validated;
  check(!pinkas::verify_online_prevalidated(empty_validated),
        "a default validated object was accepted");

  auto changed_outputs = result.Y;
  auto changed_proof = result.proof;
  auto changed_scalars = scalars;
  pinkas::ValidatedInputs owning_validated;
  check(pinkas::prepare_validated_inputs(
            parameters, bases, changed_scalars, changed_outputs,
            changed_proof, owning_validated),
        "mutable Pinkas inputs did not validate");
  Group::add(changed_outputs[0], changed_outputs[0], generator());
  Group::add(changed_proof.W[0][0], changed_proof.W[0][0], generator());
  changed_scalars[0][0] = scalar(999);
  check(pinkas::verify_online_prevalidated(owning_validated),
        "validated state retained caller-owned references");

  const auto encoded = pinkas::serialize_proof(parameters, result.proof);
  check(encoded.size() == 36616,
        std::string("the deterministic Pinkas proof size changed: ") +
            std::to_string(encoded.size()));
  const auto proof_digest = pinkas::internal::sha256(encoded);
  std::ostringstream proof_hex;
  proof_hex << std::hex;
  for (const std::uint8_t byte : proof_digest) {
    proof_hex.width(2);
    proof_hex.fill('0');
    proof_hex << static_cast<unsigned int>(byte);
  }
  constexpr std::string_view expected_proof_digest =
      "4b7ae3e7375dfcb931cf406b3bdd4df6a92757a59df0ed66f1159c2f02b6f29a";
  check(proof_hex.str() == expected_proof_digest,
        std::string("unexpected serialized proof digest: ") +
            proof_hex.str());
  pinkas::PinkasProof decoded;
  check(pinkas::deserialize_proof(parameters, encoded, decoded),
        "a valid serialized Pinkas proof did not decode");
  check(decoded.W == result.proof.W,
        "Pinkas proof changed during serialization");
  check(pinkas::verify_serialized(
            parameters, bases, scalars, result.Y, encoded),
        "serialized Pinkas verification failed");

  auto tampered_outputs = result.Y;
  Group::add(tampered_outputs[0], tampered_outputs[0], generator());
  check(!pinkas::verify(
             parameters, bases, scalars, tampered_outputs, result.proof),
        "a tampered Pinkas output was accepted");

  auto tampered_scalars = scalars;
  tampered_scalars[0][0] = scalar(12345);
  check(!pinkas::verify(
             parameters, bases, tampered_scalars, result.Y, result.proof),
        "tampered Pinkas scalars were accepted");

  auto tampered_bases = bases;
  Group::add(tampered_bases[0], tampered_bases[0], generator());
  check(!pinkas::verify(
             parameters, tampered_bases, scalars, result.Y, result.proof),
        "tampered Pinkas bases were accepted");

  auto tampered_proof = result.proof;
  Group::add(tampered_proof.W[0][0], tampered_proof.W[0][0], generator());
  check(!pinkas::verify(
             parameters, bases, scalars, result.Y, tampered_proof),
        "a tampered Pinkas proof was accepted");

  auto short_proof = result.proof;
  short_proof.W[0].pop_back();
  check(!pinkas::verify(
             parameters, bases, scalars, result.Y, short_proof),
        "a dimensionally invalid Pinkas proof was accepted");

  auto truncated = encoded;
  truncated.pop_back();
  pinkas::PinkasProof rejected;
  check(pinkas::deserialize_proof_detailed(
            parameters, truncated, rejected) == pinkas::DecodeError::Truncated,
        "a truncated Pinkas proof was accepted");
  check(rejected.W.empty(), "failed decoding left a partial proof");

  auto trailing = encoded;
  trailing.push_back(0);
  check(pinkas::deserialize_proof_detailed(
            parameters, trailing, rejected) ==
            pinkas::DecodeError::TrailingData,
        "trailing Pinkas proof bytes were accepted");

  auto wrong_version = encoded;
  wrong_version[8] ^= 1;
  check(pinkas::deserialize_proof_detailed(
            parameters, wrong_version, rejected) ==
            pinkas::DecodeError::InvalidVersion,
        "an invalid proof version was accepted");

  auto invalid_dimensions = encoded;
  // Version frame occupies 8 length bytes plus 16 protocol bytes. The first
  // dimension follows and is an eight-byte big-endian instance count.
  invalid_dimensions[31] = 0;
  check(pinkas::deserialize_proof_detailed(
            parameters, invalid_dimensions, rejected) ==
            pinkas::DecodeError::InvalidDimensions,
        "zero proof instances were accepted");

  auto invalid_point = encoded;
  // The first framed G2 encoding begins after the version and dimensions.
  invalid_point[48] ^= 0xff;
  const auto invalid_point_error = pinkas::deserialize_proof_detailed(
      parameters, invalid_point, rejected);
  check(invalid_point_error == pinkas::DecodeError::InvalidPoint ||
            invalid_point_error == pinkas::DecodeError::NonCanonicalPoint,
        "an invalid point encoding was accepted");

  auto invalid_parameters = parameters;
  invalid_parameters.domain = "wrong-domain";
  check(pinkas::deserialize_proof_detailed(
            invalid_parameters, encoded, rejected) ==
            pinkas::DecodeError::InvalidParameters,
        "decoding accepted invalid public parameters");
  check_invalid_argument(
      [&] { pinkas::serialize_proof(invalid_parameters, result.proof); },
      "serialization accepted invalid public parameters");
  check_invalid_argument(
      [&] { pinkas::serialize_proof(parameters, short_proof); },
      "serialization accepted an invalid proof shape");
}

}  // namespace

int main() {
  try {
    pippenger::initialize();
    check(pippenger::scalar_bit_length() == Scalar::getBitSize(),
          "reported scalar bit length is inconsistent");
    test_pippenger_equivalence();
    test_pippenger_rejections();
    test_pinkas_round_trip();
    std::cout << "All pippenger/Pinkas tests passed.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Test failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

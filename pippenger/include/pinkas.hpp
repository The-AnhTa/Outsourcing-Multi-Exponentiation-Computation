#pragma once

#include "pippenger.hpp"

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace pinkas {

using Scalar = pippenger::Scalar;
using Group = pippenger::Group;
using ScalarMatrix = pippenger::ScalarMatrix;
using Bytes = std::vector<std::uint8_t>;

struct PublicParameters {
  std::string group_identifier;
  std::string domain;
  std::size_t scalar_bit_length = 0;
  std::size_t challenge_bit_length = 0;
  std::size_t msm_window_width = 0;
};

struct PinkasProof {
  std::vector<std::vector<Group>> W;
};

struct ProverResult {
  std::vector<Group> Y;
  PinkasProof proof;
};

struct ValidatedInputs {
  const PublicParameters* parameters = nullptr;
  const std::vector<Group>* bases = nullptr;
  const ScalarMatrix* scalars = nullptr;
  const std::vector<Group>* outputs = nullptr;
  const PinkasProof* proof = nullptr;
};

PublicParameters setup(std::size_t msm_window_width);

ProverResult prove(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars);

bool prepare_validated_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof,
    ValidatedInputs& validated);

bool verify_online_prevalidated(const ValidatedInputs& validated);

bool verify(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof);

Bytes serialize_proof(
    const PublicParameters& parameters,
    const PinkasProof& proof);

bool deserialize_proof(
    const PublicParameters& parameters,
    std::span<const std::uint8_t> encoded,
    PinkasProof& proof);



bool verify_serialized(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    std::span<const std::uint8_t> encoded_proof);

}

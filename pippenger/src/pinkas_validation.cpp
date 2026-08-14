#include "internal/validation.hpp"

#include "internal/constants.hpp"

namespace pinkas::internal {
namespace {

bool valid_points(const std::vector<Group>& points) {
  for (const Group& point : points) {
    if (!valid_point(point)) return false;
  }
  return true;
}

bool valid_scalar_dimensions(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  if (bases.empty() || scalars.empty()) return false;
  for (const auto& row : scalars) {
    if (row.size() != bases.size()) return false;
  }
  return true;
}

}  // namespace

bool valid_point(const Group& point) {
  return point.isValid() && point.isValidOrder();
}

bool valid_parameters(const PublicParameters& parameters) {
  try {
    pippenger::initialize();
    if (parameters.group_identifier != kGroupIdentifier ||
        parameters.domain != kProtocolVersion ||
        parameters.scalar_bit_length != Scalar::getBitSize() ||
        parameters.challenge_bit_length != kChallengeBitLength) {
      return false;
    }
    static_cast<void>(
        pippenger::number_of_windows(parameters.msm_window_width));
    return true;
  } catch (...) {
    return false;
  }
}

bool valid_prover_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  return valid_parameters(parameters) &&
         valid_scalar_dimensions(bases, scalars) && valid_points(bases);
}

bool valid_proof(
    const PublicParameters& parameters,
    const PinkasProof& proof) {
  if (!valid_parameters(parameters) || proof.W.empty()) return false;
  for (const auto& row : proof.W) {
    if (row.size() != parameters.scalar_bit_length || !valid_points(row)) {
      return false;
    }
  }
  return true;
}

bool valid_verification_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  if (!valid_prover_inputs(parameters, bases, scalars) ||
      outputs.size() != scalars.size() ||
      proof.W.size() != scalars.size() || !valid_points(outputs)) {
    return false;
  }
  for (const auto& row : proof.W) {
    if (row.size() != parameters.scalar_bit_length || !valid_points(row)) {
      return false;
    }
  }
  return true;
}

}  // namespace pinkas::internal

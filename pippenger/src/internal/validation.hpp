#pragma once

#include "pinkas/pinkas.hpp"

#include <vector>

namespace pippenger::internal {

void validate_msm_dimensions(
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars);

}  // namespace pippenger::internal

namespace pinkas::internal {

bool valid_point(const Group& point);
bool valid_parameters(const PublicParameters& parameters);

bool valid_prover_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars);

bool valid_verification_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof);

bool valid_proof(
    const PublicParameters& parameters,
    const PinkasProof& proof);

}  // namespace pinkas::internal

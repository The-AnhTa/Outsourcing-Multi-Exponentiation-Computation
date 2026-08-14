#pragma once

#include "blsagg/protocol.hpp"

namespace blsagg::internal {

Digest parameter_digest(const PublicParameters& parameters);
Digest precomputation_digest(const PublicParameters& parameters,
                             const Precomputation& precomputation);
Digest context_binding(const PublicParameters& parameters,
                       const Precomputation& precomputation,
                       const Statement& statement,
                       std::span<const G1> message_points);
Digest validated_proof_binding(const PublicParameters& parameters,
                               const Proof& proof);

bool valid_public_parameters(const PublicParameters& parameters);
bool valid_statement(const PublicParameters& parameters,
                     const Statement& statement);
bool valid_precomputation(const PublicParameters& parameters,
                          const Precomputation& precomputation);
bool valid_proof(const PublicParameters& parameters, const Proof& proof);

}  // namespace blsagg::internal

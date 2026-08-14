#include "internal/validation.hpp"

#include "blsagg/serialization.hpp"
#include "blsagg/transcript.hpp"
#include "internal/crypto.hpp"

#include <limits>
#include <unordered_set>

namespace blsagg::internal {
namespace {

struct BytesHash {
  std::size_t operator()(const Bytes& bytes) const noexcept {
    // FNV-1a is sufficient here: this hash table is only an optimization for
    // exact byte-vector equality, not a cryptographic commitment.
    std::size_t hash = sizeof(std::size_t) == 8
                           ? static_cast<std::size_t>(14695981039346656037ull)
                           : static_cast<std::size_t>(2166136261u);
    const std::size_t prime = sizeof(std::size_t) == 8
                                  ? static_cast<std::size_t>(1099511628211ull)
                                  : static_cast<std::size_t>(16777619u);
    for (const auto byte : bytes) {
      hash ^= byte;
      hash *= prime;
    }
    return hash;
  }
};

bool distinct_messages(const std::vector<Bytes>& messages) {
  std::unordered_set<Bytes, BytesHash> seen;
  seen.reserve(messages.size());
  for (const auto& message : messages) {
    if (!seen.insert(message).second) return false;
  }
  return true;
}

bool claims_equal(const RexpClaims& lhs, const RexpClaims& rhs) {
  return lhs.E == rhs.E && lhs.F == rhs.F && lhs.TL == rhs.TL &&
         lhs.TR == rhs.TR;
}

bool precomputations_equal(const Precomputation& lhs,
                           const Precomputation& rhs) {
  return lhs.gamma_chain == rhs.gamma_chain &&
         lhs.lambda_chain == rhs.lambda_chain && lhs.X == rhs.X &&
         lhs.delta1L == rhs.delta1L && lhs.delta1R == rhs.delta1R &&
         lhs.delta2L == rhs.delta2L && lhs.delta2R == rhs.delta2R &&
         claims_equal(lhs.g1_round0, rhs.g1_round0) &&
         claims_equal(lhs.g2_round0, rhs.g2_round0) &&
         lhs.digest == rhs.digest;
}

}  // namespace

Digest parameter_digest(const PublicParameters& parameters) {
  Bytes bytes;
  append_frame(bytes, "bls-agg-bf/pp/v1");
  append_frame(bytes, "BN254/mcl-v3.00");
  append_frame(bytes, parameters.mode == AggregationMode::BasicDistinct
                          ? "basic-distinct"
                          : "augmented");
  append_frame(bytes, encode_u64(parameters.k));
  append_frame(bytes, encode_u64(parameters.d));
  append_frame(bytes, serialize(parameters.H));
  for (const auto& point : parameters.Gamma)
    append_frame(bytes, serialize(point));
  for (const auto& point : parameters.Lambda)
    append_frame(bytes, serialize(point));
  append_frame(bytes, serialize(parameters.L));
  append_frame(bytes, serialize(parameters.Lprime));
  return sha256(bytes);
}

Digest precomputation_digest(const PublicParameters& parameters,
                             const Precomputation& precomputation) {
  Bytes bytes;
  append_frame(bytes, "bls-agg-bf/aux/v1");
  append_frame(bytes, parameters.digest);
  for (const auto& level : precomputation.gamma_chain)
    for (const auto& point : level) append_frame(bytes, serialize(point));
  for (const auto& level : precomputation.lambda_chain)
    for (const auto& point : level) append_frame(bytes, serialize(point));
  for (const auto* values : {&precomputation.X, &precomputation.delta1L,
                             &precomputation.delta1R,
                             &precomputation.delta2L,
                             &precomputation.delta2R})
    for (const auto& value : *values) append_frame(bytes, serialize(value));
  for (const auto* claims : {&precomputation.g1_round0,
                             &precomputation.g2_round0}) {
    append_frame(bytes, serialize(claims->E));
    append_frame(bytes, serialize(claims->F));
    append_frame(bytes, serialize(claims->TL));
    append_frame(bytes, serialize(claims->TR));
  }
  return sha256(bytes);
}

Digest context_binding(const PublicParameters& parameters,
                       const Precomputation& precomputation,
                       const Statement& statement,
                       std::span<const G1> message_points) {
  Bytes bytes;
  append_frame(bytes, "bls-agg-bf/validated-context/v1");
  append_frame(bytes, parameters.digest);
  append_frame(bytes, precomputation.digest);
  append_frame(bytes, serialize(statement.sigma_agg));
  for (const auto& message : statement.messages) append_frame(bytes, message);
  for (const auto& key : statement.public_keys)
    append_frame(bytes, serialize(key));
  for (const auto& point : message_points)
    append_frame(bytes, serialize(point));
  return sha256(bytes);
}

Digest validated_proof_binding(const PublicParameters& parameters,
                               const Proof& proof) {
  Bytes bytes;
  append_frame(bytes, "bls-agg-bf/validated-proof/v1");
  append_frame(bytes, parameters.digest);
  append_frame(bytes, serialize_proof(parameters, proof));
  return sha256(bytes);
}

bool valid_public_parameters(const PublicParameters& parameters) {
  if (parameters.d < 1 ||
      parameters.d >= std::numeric_limits<std::size_t>::digits ||
      parameters.k != (std::size_t{1} << parameters.d) ||
      parameters.Gamma.size() != parameters.k ||
      parameters.Lambda.size() != parameters.k ||
      (parameters.mode != AggregationMode::BasicDistinct &&
       parameters.mode != AggregationMode::Augmented))
    return false;
  if (!valid(parameters.H, true) || !valid(parameters.L, true) ||
      !valid(parameters.Lprime, true))
    return false;
  for (const auto& point : parameters.Gamma)
    if (!valid(point, true)) return false;
  for (const auto& point : parameters.Lambda)
    if (!valid(point, true)) return false;
  return parameter_digest(parameters) == parameters.digest;
}

bool valid_statement(const PublicParameters& parameters,
                     const Statement& statement) {
  if (statement.messages.size() != parameters.k ||
      statement.public_keys.size() != parameters.k ||
      !valid(statement.sigma_agg, true))
    return false;
  if (parameters.mode == AggregationMode::BasicDistinct &&
      !distinct_messages(statement.messages))
    return false;
  for (const auto& key : statement.public_keys)
    if (!valid(key, true)) return false;
  return true;
}

bool valid_precomputation(const PublicParameters& parameters,
                          const Precomputation& precomputation) {
  if (!valid_public_parameters(parameters) ||
      precomputation.digest !=
          precomputation_digest(parameters, precomputation))
    return false;
  try {
    return precomputations_equal(precomputation, precompute(parameters));
  } catch (...) {
    return false;
  }
}

bool valid_proof(const PublicParameters& parameters, const Proof& proof) {
  if (parameters.d < 1 ||
      proof.g1_rexp_claims.size() != parameters.d - 1 ||
      proof.g2_rexp_claims.size() != parameters.d - 1 ||
      proof.dory_steps.size() != parameters.d ||
      proof.insert_g1_u.size() != parameters.d ||
      proof.insert_g2_u.size() != parameters.d)
    return false;
  for (const auto* value :
       {&proof.cm_M, &proof.cm_pk, &proof.T, &proof.U1, &proof.U2})
    if (!valid(*value)) return false;
  if (!valid(proof.R_Gamma) || !valid(proof.R_Lambda) ||
      !valid(proof.Phi_final) || !valid(proof.Theta_final))
    return false;
  for (const auto* claims : {&proof.g1_rexp_claims,
                             &proof.g2_rexp_claims})
    for (const auto& claim : *claims)
      if (!valid(claim.E) || !valid(claim.F) || !valid(claim.TL) ||
          !valid(claim.TR))
        return false;
  for (const auto& step : proof.dory_steps)
    for (const auto* value :
         {&step.A1L, &step.A1R, &step.A2L, &step.A2R, &step.W1,
          &step.W2})
      if (!valid(*value)) return false;
  for (const auto& value : proof.insert_g1_u)
    if (!valid(value)) return false;
  for (const auto& value : proof.insert_g2_u)
    if (!valid(value)) return false;
  return true;
}

}  // namespace blsagg::internal

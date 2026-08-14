#include "pinkas/pinkas.hpp"

#include "internal/protocol_math.hpp"
#include "internal/transcript.hpp"
#include "internal/validation.hpp"

#include <utility>
#include <vector>

namespace pinkas {
namespace {

bool verify_arithmetic(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof) {
  const std::size_t instances = scalars.size();
  const std::size_t bits = parameters.scalar_bit_length;
  const internal::Digest transcript =
      internal::build_transcript(parameters, bases, scalars, outputs, proof);
  const auto challenges =
      internal::derive_challenges(parameters, transcript, instances);
  const std::size_t aggregate_exponent_bits =
      parameters.challenge_bit_length + internal::ceil_log2(bits);

  for (std::size_t instance = 0; instance < instances; ++instance) {
    const Group proof_sum = internal::pippenger_msm(
        proof.W[instance], challenges[instance],
        parameters.msm_window_width, parameters.challenge_bit_length);

    std::vector<Scalar> exponents(bases.size());
    for (Scalar& exponent : exponents) exponent.clear();
    for (std::size_t base = 0; base < bases.size(); ++base) {
      const auto bytes =
          internal::canonical_scalar_bytes(scalars[instance][base]);
      Scalar sum;
      sum.clear();
      for (std::size_t bit = 0; bit < bits; ++bit) {
        if (internal::scalar_bit(bytes, bit)) {
          sum = internal::scalar_add(sum, challenges[instance][bit]);
        }
      }
      exponents[base] = sum;
    }
    const Group statement_sum = internal::pippenger_msm(
        bases, exponents, parameters.msm_window_width,
        aggregate_exponent_bits);
    if (proof_sum != statement_sum ||
        outputs[instance] !=
            internal::horner_reconstruct(proof.W[instance])) {
      return false;
    }
  }
  return true;
}

}  // namespace

bool prepare_validated_inputs(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars,
    const std::vector<Group>& outputs,
    const PinkasProof& proof,
    ValidatedInputs& validated) {
  validated = {};
  try {
    if (!internal::valid_verification_inputs(
            parameters, bases, scalars, outputs, proof)) {
      return false;
    }
    ValidatedInputs result;
    result.parameters_ = parameters;
    result.bases_ = bases;
    result.scalars_ = scalars;
    result.outputs_ = outputs;
    result.proof_ = proof;
    result.ready_ = true;
    validated = std::move(result);
    return true;
  } catch (...) {
    validated = {};
    return false;
  }
}

bool verify_online_prevalidated(const ValidatedInputs& validated) {
  if (!validated.ready_) return false;
  try {
    return verify_arithmetic(
        validated.parameters_, validated.bases_, validated.scalars_,
        validated.outputs_, validated.proof_);
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

}  // namespace pinkas

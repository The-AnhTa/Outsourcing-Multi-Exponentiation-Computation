#include "pinkas/pinkas.hpp"

#include "internal/protocol_math.hpp"
#include "internal/transcript.hpp"
#include "internal/validation.hpp"

#include <stdexcept>
#include <utility>

namespace pinkas {

ProverResult prove(
    const PublicParameters& parameters,
    const std::vector<Group>& bases,
    const ScalarMatrix& scalars) {
  if (!internal::valid_prover_inputs(parameters, bases, scalars)) {
    throw std::invalid_argument("invalid Pinkas prover inputs");
  }

  const std::size_t instances = scalars.size();
  const std::size_t bits = parameters.scalar_bit_length;
  PinkasProof proof;
  proof.W.resize(instances, std::vector<Group>(bits));
  for (auto& row : proof.W) {
    for (Group& point : row) point.clear();
  }
  std::vector<Group> outputs(instances);

  for (std::size_t instance = 0; instance < instances; ++instance) {
    for (std::size_t base = 0; base < bases.size(); ++base) {
      const auto bytes =
          internal::canonical_scalar_bytes(scalars[instance][base]);
      for (std::size_t bit = 0; bit < bits; ++bit) {
        if (internal::scalar_bit(bytes, bit)) {
          proof.W[instance][bit] =
              internal::add(proof.W[instance][bit], bases[base]);
        }
      }
    }
    outputs[instance] = internal::horner_reconstruct(proof.W[instance]);
  }
  return {std::move(outputs), std::move(proof)};
}

}  // namespace pinkas

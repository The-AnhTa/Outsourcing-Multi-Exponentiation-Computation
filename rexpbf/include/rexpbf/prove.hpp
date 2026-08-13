#pragma once
#include "rexpbf/proof.hpp"
#include "rexpbf/setup.hpp"
#include "rexpbf/transcript.hpp"

namespace rexpbf {
ProveResult prove(const CRS& crs, const Precomputation& precomputation,
                  const Statement& statement, const ProverInput& prover_input);
ChallengeTrace replay_challenges(const CRS& crs, const Statement& statement,
                                 const Proof& proof,
                                 TranscriptMetrics* metrics = nullptr);
std::vector<std::uint8_t> serialize_proof_payload(const Proof& proof);
std::vector<std::uint8_t> serialize_proof_wire(const Proof& proof, std::size_t d, std::size_t n);
Proof deserialize_proof_wire(std::span<const std::uint8_t> bytes, std::size_t expected_d,
                             std::size_t expected_n);
}

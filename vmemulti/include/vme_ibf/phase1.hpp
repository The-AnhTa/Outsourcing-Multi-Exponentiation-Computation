#pragma once
#include "vme_ibf/setup.hpp"
#include "vme_ibf/transcript.hpp"
namespace vme_ibf {
std::vector<Fr> tensor_vector(std::span<const Fr> r);
Digest compute_statement_digest(const VmeIbfCRS&, const VmeIbfStatementInput&, const G2& X);
Phase1Result prove_phase1(const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeIbfStatementInput&);
Phase1Result prove_phase1_core(const VmeIbfCRS&, const VmeIbfPrecomputation&,
                              std::span<const Fr> aggregate_exponents,
                              const G2& aggregate_output, Transcript transcript);
std::vector<Fr> replay_rho(const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeIbfStatement&, std::span<const RexpClaims>, const G1& R, Digest* after_R = nullptr);
}

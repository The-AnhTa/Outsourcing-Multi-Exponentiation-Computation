#pragma once

#include "vme_ibf/setup.hpp"
#include "vme_ibf/transcript.hpp"

namespace vme_ibf {

Digest compute_statement_digest(const VmeIbfCRS&,
                                const VmeIbfStatementInput&,
                                const G2& X);
Phase1Result prove_phase1(const VmeIbfCRS&,
                          const VmeIbfPrecomputation&,
                          const VmeIbfStatementInput&);
std::vector<Fr> replay_rho(const VmeIbfCRS&,
                           const VmeIbfPrecomputation&,
                           const VmeIbfStatement&,
                           std::span<const RexpClaims>,
                           const G1& R,
                           Digest* after_R = nullptr);

}

#pragma once
#include "vpip_bf/setup.hpp"
#include "vpip_bf/transcript.hpp"
namespace vpip_bf {
std::vector<Fr> tensor_vector(std::span<const Fr> r);
Digest compute_statement_digest(const VpipBfCRS&, const VpipBfPrecomputation&,
                                const VpipBfStatementInput&, const GT& C);
Phase1Result prove_phase1(const VpipBfCRS&, const VpipBfPrecomputation&, const VpipBfStatementInput&);
std::vector<Fr> replay_rho(const VpipBfCRS&, const VpipBfPrecomputation&, const VpipBfStatement&, std::span<const RexpClaims>, const G1& R, Digest* after_R = nullptr);
}

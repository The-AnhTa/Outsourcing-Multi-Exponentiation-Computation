#pragma once

#include "protocol.hpp"
#include "vme_ibf/symbolic_gt.hpp"

namespace vme_ibf::internal {

struct EquationTimings {
    double transcript_ms{};
    double batch_inversion_ms{};
    double recurrence_ms{};
    double terminal_assembly_ms{};
};

struct VerificationEquations {
    ProtocolChallenges challenges;
    SymbolicGtExpression dory;
    SymbolicGtExpression rexp;
    G1 terminal_dory_g1;
    G2 terminal_dory_g2;
    EquationTimings timings;
};

bool validate_verification_objects(const VmeIbfCRS&,
                                   const VmeIbfPrecomputation&,
                                   const VmeIbfStatement&,
                                   const VmeIbfProof&,
                                   bool audit_precomputation_values);

bool build_verification_equations(const VmeIbfCRS&,
                                  const VmeIbfPrecomputation&,
                                  const VmeIbfStatement&,
                                  const VmeIbfProof&,
                                  VerificationEquations&);

const GT& resolve_gt_atom(GtAtomId, const VmeIbfPrecomputation&,
                          const VmeIbfProof&);
PairingInputs resolve_pairing_atom(PairingAtomId, const VmeIbfCRS&,
                                   const VmeIbfStatement&,
                                   const VmeIbfProof&,
                                   const VerificationEquations&);

} // namespace vme_ibf::internal

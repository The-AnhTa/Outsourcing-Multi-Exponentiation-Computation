#pragma once
#include "vme_ibf/types.hpp"

namespace vme_ibf {
struct VmeIbfProof {
  std::vector<RexpClaims> rexp_claims;
  G1 R;
  std::vector<DoryFoldProof> dory_folds;
  std::vector<GT> batch_U;
  G1 PhiFinal;
  G2 ThetaFinal;
};

VmeIbfProof assemble_public_proof(const Phase1Result&, const Phase2Result&);
}

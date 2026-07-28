#pragma once
#include "vpip_bf/types.hpp"

namespace vpip_bf {
struct VpipBfProof {
  std::vector<RexpClaims> rexp_claims;
  G1 R;
  std::vector<DoryFoldProof> dory_folds;
  std::vector<GT> batch_U;
  G1 PhiFinal;
  G2 ThetaFinal;
};
using Proof = VpipBfProof;

VpipBfProof assemble_public_proof(const Phase1Result&, const Phase2Result&);
}

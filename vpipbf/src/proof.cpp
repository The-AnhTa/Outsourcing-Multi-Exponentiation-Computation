#include "vpip_bf/proof.hpp"

namespace vpip_bf {
VpipBfProof assemble_public_proof(const Phase1Result&p1,const Phase2Result&p2){
 VpipBfProof p;p.rexp_claims=p1.dynamic_claims;p.R=p1.R;p.dory_folds=p2.proof.dory_folds;p.batch_U=p2.proof.batch_U;p.PhiFinal=p2.proof.PhiFinal;p.ThetaFinal=p2.proof.ThetaFinal;return p;
}
}

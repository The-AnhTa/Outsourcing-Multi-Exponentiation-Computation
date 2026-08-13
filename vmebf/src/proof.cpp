#include "vme_ibf/proof.hpp"

#include <stdexcept>

namespace vme_ibf {
VmeIbfProof assemble_public_proof(const Phase1Result&p1,const Phase2Result&p2){
 if(p2.proof.dory_folds.size()!=p1.dynamic_claims.size()+1||p2.proof.dory_folds.size()!=p2.proof.batch_U.size())throw std::invalid_argument("inconsistent proof phases");
 VmeIbfProof p;p.rexp_claims=p1.dynamic_claims;p.R=p1.R;p.dory_folds=p2.proof.dory_folds;p.batch_U=p2.proof.batch_U;p.PhiFinal=p2.proof.PhiFinal;p.ThetaFinal=p2.proof.ThetaFinal;return p;
}
}

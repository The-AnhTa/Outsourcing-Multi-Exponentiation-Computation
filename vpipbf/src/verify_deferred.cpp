#include "vpip_bf/verify_deferred.hpp"
#include "vpip_bf/verify_combined.hpp"
#include "vpip_bf/verify_online.hpp"
#include "vpip_bf/verify_reference.hpp"

namespace vpip_bf {
ReferenceVerificationTrace verify_core_symbolic_unchecked(const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&,OnlineTimingBreakdown* = nullptr);
std::optional<ValidatedVerificationInputs> PrevalidateVerificationInputs(
    const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&v){
  if(!validate_verification_inputs(c,p,s,v))return std::nullopt;
  return ValidatedVerificationInputs(c,p,s,v);
}

bool verify_online(const ValidatedVerificationInputs&i){
  return verify_core_symbolic_unchecked(i.crs(),i.precomputation(),i.statement(),i.proof()).accepted;
}
bool verify_online_with_timing(const ValidatedVerificationInputs&i,OnlineTimingBreakdown&t){
  t={};return verify_core_symbolic_unchecked(i.crs(),i.precomputation(),i.statement(),i.proof(),&t).accepted;
}
DeferredVerificationTrace verify_deferred_with_trace(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&v){auto inputs=PrevalidateVerificationInputs(c,p,s,v);DeferredVerificationTrace z;if(!inputs)return z;auto r=verify_core_unchecked(c,p,s,v);z.accepted=r.accepted;z.dory_accepted=r.dory_accepted;z.rexp_accepted=r.rexp_accepted;z.rho=r.rho;z.beta=r.beta;z.alpha=r.alpha;z.gamma=r.gamma;z.epsilon=r.epsilon;z.dory_residual_reference=r.dory_residual;z.dory_residual_pippenger=r.dory_residual;z.rexp_residual_reference=r.rexp_residual;z.rexp_residual_pippenger=r.rexp_residual;return z;}
bool verify_deferred(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&v){auto i=PrevalidateVerificationInputs(c,p,s,v);return i&&verify_online(*i);}
CombinedVerificationTrace verify_deferred_combined_with_trace(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&v){CombinedVerificationTrace z;auto i=PrevalidateVerificationInputs(c,p,s,v);if(!i)return z;auto r=verify_core_unchecked(c,p,s,v);z.accepted=r.accepted;z.rho=r.rho;z.beta=r.beta;z.alpha=r.alpha;z.gamma=r.gamma;z.epsilon=r.epsilon;z.evaluated_dory_residual=r.dory_residual;z.evaluated_rexp_residual=r.rexp_residual;return z;}
bool verify_deferred_combined(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&v){return verify_deferred(c,p,s,v);}
CombinedVerificationTrace verify_online_with_trace(const ValidatedVerificationInputs&i){CombinedVerificationTrace z;auto r=verify_core_symbolic_unchecked(i.crs(),i.precomputation(),i.statement(),i.proof());z.accepted=r.accepted;z.rho=r.rho;z.beta=r.beta;z.alpha=r.alpha;z.gamma=r.gamma;z.epsilon=r.epsilon;z.evaluated_dory_residual=r.dory_residual;z.evaluated_rexp_residual=r.rexp_residual;return z;}
}

#include "vpip_bf/verify_reference.hpp"
#include "internal/verification.hpp"
#include "internal/protocol.hpp"
#include "vpip_bf/serialization.hpp"
#include "vpip_bf/setup.hpp"
#include "vpip_bf/transcript.hpp"
#include "vpip_bf/gt_multiexp.hpp"
#include <array>
#include <atomic>
#include <chrono>

namespace vpip_bf {
namespace {
using Clock=std::chrono::steady_clock;
double elapsed(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
std::atomic<std::size_t> core_calls{};
Fr fm(const Fr&a,const Fr&b){Fr z;Fr::mul(z,a,b);return z;}
GT prod(std::initializer_list<GT>x){GT z;z.setOne();for(const auto&q:x)z=gt_mul(z,q);return z;}
}

void reset_verification_core_call_count_for_testing(){core_calls=0;}
std::size_t verification_core_call_count_for_testing(){return core_calls.load();}

bool validate_verification_inputs(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&proof){
 try{
  return validate_crs(c)&&audit_precomputation(c,p)
      &&validate_statement_shape(c,s)&&validate_statement_elements(s)
      &&validate_statement_digest(c,s)&&validate_proof_shape(c,proof)
      &&validate_proof_elements(proof);
 }catch(...){return false;}
}

ReferenceVerificationTrace verify_core_unchecked(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&proof,OnlineTimingBreakdown*timing){
 ReferenceVerificationTrace tr;core_calls.fetch_add(1);auto total_start=Clock::now();
 try{
  auto q0=Clock::now();Transcript T(s.digest);auto q1=Clock::now();if(timing)timing->transcript_init_ms+=elapsed(q0,q1);
  tr.rho.resize(c.d);std::vector<DoryTargetState>fresh(c.d+1);GT outer=p.pairing_x[0];
  for(size_t j=0;j<c.d;++j){RexpClaims msg=internal::rexp_claim(j,p,proof.rexp_claims);auto a=Clock::now();internal::absorb_rexp_claim(T,j,c.n>>j,msg);auto b=Clock::now();Fr rho=T.challenge_nonzero("vpipbf/challenge/rexp-r/v1",j);auto e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ri=inverse_nonzero(rho);tr.rho[j]=rho;size_t t=j+1;a=Clock::now();fresh[t].D0=prod({outer,gt_pow(msg.E,rho),gt_pow(msg.F,ri)});fresh[t].D1=gt_mul(msg.TL,gt_pow(msg.TR,rho));fresh[t].D2=gt_mul(p.pairing_x[t],gt_pow(p.delta2R[j],ri));outer=fresh[t].D1;b=Clock::now();if(timing)timing->gt_multiexp_ms+=elapsed(a,b);}
  auto a=Clock::now();T.absorb("vpipbf/rexp-result-R/v1",serialize(proof.R));auto b=Clock::now();if(timing)timing->transcript_replay_ms+=elapsed(a,b);
  std::vector<Fr>r(c.d);for(size_t j=0;j<c.d;++j)r[c.d-j-1]=tr.rho[j];a=Clock::now();auto weights=tensor_vector(r);b=Clock::now();if(timing)timing->tensor_reconstruction_ms+=elapsed(a,b);if(weights.size()!=c.n)return {};
  a=Clock::now();G1 Y=g1_multiexp(s.X,weights);b=Clock::now();if(timing)timing->g1_msm_y_ms+=elapsed(a,b);
  std::array<Bytes,3>ini{encode_u64_be(c.d),encode_u64_be(c.n),serialize(Y)};a=Clock::now();T.absorb("vpipbf/initial-targets/v1",ini);b=Clock::now();if(timing)timing->transcript_replay_ms+=elapsed(a,b);
  DoryTargetState agg;mcl::bn::pairing(agg.D0,Y,c.Lprime);agg.D1=s.C;mcl::bn::pairing(agg.D2,proof.R,c.Lprime);tr.beta.resize(c.d);tr.alpha.resize(c.d);tr.gamma.resize(c.d+1);tr.gamma[0].clear();
  for(size_t t=1;t<=c.d;++t){size_t k=t-1,m=c.n>>k,h=m/2;const auto&f=proof.dory_folds[k];std::array<Bytes,6>bm{encode_u64_be(k),encode_u64_be(m),serialize(f.D1L),serialize(f.D1R),serialize(f.D2L),serialize(f.D2R)};a=Clock::now();T.absorb("VPIP.BF/DORY-BETA-MESSAGE/V2",bm);b=Clock::now();Fr beta=T.challenge_nonzero("VPIP.BF/BETA/V2",k);auto e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr bi=inverse_nonzero(beta);tr.beta[k]=beta;std::array<Bytes,4>am{encode_u64_be(k),encode_u64_be(m),serialize(f.W1),serialize(f.W2)};a=Clock::now();T.absorb("VPIP.BF/DORY-ALPHA-MESSAGE/V2",am);b=Clock::now();Fr alpha=T.challenge_nonzero("VPIP.BF/ALPHA/V2",k);e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ai=inverse_nonzero(alpha);tr.alpha[k]=alpha;
    a=Clock::now();DoryTargetState bar;std::array<GT,4>d0b{agg.D1,agg.D2,f.W1,f.W2};std::array<Fr,4>d0e{bi,beta,alpha,ai};bar.D0=prod({agg.D0,p.pairing_x[k],gt_multiexp_native(d0b,d0e)});std::array<GT,3>d1b{f.D1L,p.pairing_x[k+1],p.delta1R[k]};std::array<Fr,3>d1e{alpha,fm(alpha,beta),beta};bar.D1=gt_mul(f.D1R,gt_multiexp_native(d1b,d1e));std::array<GT,3>d2b{f.D2L,p.pairing_x[k+1],p.delta2R[k]};std::array<Fr,3>d2e{ai,fm(ai,bi),bi};bar.D2=gt_mul(f.D2R,gt_multiexp_native(d2b,d2e));b=Clock::now();if(timing)timing->dory_recurrence_ms+=elapsed(a,b);
    const GT&u=proof.batch_U[k];std::array<Bytes,3>um{encode_u64_be(t),encode_u64_be(h),serialize(u)};a=Clock::now();T.absorb("VPIP.BF/BATCH-U/V2",um);b=Clock::now();Fr gamma=T.challenge_nonzero("VPIP.BF/GAMMA/V2",t);e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}tr.gamma[t]=gamma;Fr gs=fm(gamma,gamma);a=Clock::now();std::array<GT,2>batch0b{bar.D0,u};std::array<Fr,2>batch0e{gs,gamma};agg.D0=gt_mul(gt_multiexp_native(batch0b,batch0e),fresh[t].D0);agg.D1=gt_mul(gt_pow(bar.D1,gamma),fresh[t].D1);agg.D2=gt_mul(gt_pow(bar.D2,gamma),fresh[t].D2);b=Clock::now();if(timing)timing->dory_recurrence_ms+=elapsed(a,b);
  }
  std::array<Bytes,2>fin{serialize(proof.PhiFinal),serialize(proof.ThetaFinal)};a=Clock::now();T.absorb("VPIP.BF/DORY-FINAL/V2",fin);b=Clock::now();Fr epsilon=T.challenge_nonzero("VPIP.BF/EPSILON/V2",c.d);auto e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ei=inverse_nonzero(epsilon);tr.epsilon=epsilon;
  a=Clock::now();GT lhs=prod({agg.D0,gt_pow(agg.D1,ei),gt_pow(agg.D2,epsilon),p.pairing_x[c.d]});G1 g1=g1_add(proof.PhiFinal,g1_pow(c.G[0],epsilon));G2 g2=g2_add(proof.ThetaFinal,g2_pow(c.H[0],ei));GT rhs;mcl::bn::pairing(rhs,g1,g2);GT rhs_inv;GT::inv(rhs_inv,rhs);tr.dory_residual=gt_mul(lhs,rhs_inv);tr.dory_accepted=lhs==rhs;b=Clock::now();if(timing)timing->terminal_dory_ms+=elapsed(a,b);
  a=Clock::now();GT rexp;mcl::bn::pairing(rexp,proof.R,c.H[0]);GT outer_inv;GT::inv(outer_inv,outer);tr.rexp_residual=gt_mul(rexp,outer_inv);tr.rexp_accepted=rexp==outer;b=Clock::now();if(timing)timing->terminal_rexp_ms+=elapsed(a,b);
  tr.final_aggregate=agg;tr.final_rexp_d1=outer;tr.accepted=tr.dory_accepted&&tr.rexp_accepted;
  if(timing){timing->online_verify_ms=elapsed(total_start,Clock::now());double known=timing->transcript_init_ms+timing->transcript_replay_ms+timing->challenge_derivation_ms+timing->tensor_reconstruction_ms+timing->g1_msm_y_ms+timing->dory_recurrence_ms+timing->gt_multiexp_ms+timing->terminal_dory_ms+timing->terminal_rexp_ms;timing->other_online_ms=timing->online_verify_ms-known;}
  return tr;
 }catch(...){if(timing)timing->online_verify_ms=elapsed(total_start,Clock::now());return {};}
}




ReferenceVerificationTrace verify_core_symbolic_unchecked(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&proof,OnlineTimingBreakdown*timing){
 ReferenceVerificationTrace tr;core_calls.fetch_add(1);auto total_start=Clock::now();
 try{
  auto q0=Clock::now();Transcript T(s.digest);auto q1=Clock::now();if(timing)timing->transcript_init_ms+=elapsed(q0,q1);
  tr.rho.resize(c.d);std::vector<DoryTargetState>fresh(c.d+1);GT outer=p.pairing_x[0];
  for(size_t j=0;j<c.d;++j){
   RexpClaims msg=internal::rexp_claim(j,p,proof.rexp_claims);
   auto a=Clock::now();internal::absorb_rexp_claim(T,j,c.n>>j,msg);auto b=Clock::now();Fr rho=T.challenge_nonzero("vpipbf/challenge/rexp-r/v1",j);auto e=Clock::now();
   if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ri=inverse_nonzero(rho);tr.rho[j]=rho;size_t t=j+1;a=Clock::now();
   fresh[t].D0=prod({outer,gt_pow(msg.E,rho),gt_pow(msg.F,ri)});fresh[t].D1=gt_mul(msg.TL,gt_pow(msg.TR,rho));fresh[t].D2=gt_mul(p.pairing_x[t],gt_pow(p.delta2R[j],ri));outer=fresh[t].D1;
   b=Clock::now();if(timing)timing->gt_multiexp_ms+=elapsed(a,b);
  }
  auto a=Clock::now();T.absorb("vpipbf/rexp-result-R/v1",serialize(proof.R));auto b=Clock::now();if(timing)timing->transcript_replay_ms+=elapsed(a,b);
  std::vector<Fr>r(c.d);for(size_t j=0;j<c.d;++j)r[c.d-j-1]=tr.rho[j];a=Clock::now();auto weights=tensor_vector(r);b=Clock::now();if(timing)timing->tensor_reconstruction_ms+=elapsed(a,b);if(weights.size()!=c.n)return {};
  a=Clock::now();G1 Y=g1_multiexp(s.X,weights);b=Clock::now();if(timing)timing->g1_msm_y_ms+=elapsed(a,b);
  std::array<Bytes,3>ini{encode_u64_be(c.d),encode_u64_be(c.n),serialize(Y)};a=Clock::now();T.absorb("vpipbf/initial-targets/v1",ini);b=Clock::now();if(timing)timing->transcript_replay_ms+=elapsed(a,b);

  DoryTargetState initial;mcl::bn::pairing(initial.D0,Y,c.Lprime);initial.D1=s.C;mcl::bn::pairing(initial.D2,proof.R,c.Lprime);
  const size_t initial0=0,px0=3,delta10=px0+c.d+1,delta20=delta10+c.d,fold0=delta20+c.d,u0=fold0+6*c.d,fresh0=u0+c.d,atom_count=fresh0+3*c.d;
  std::vector<GT>atoms(atom_count);atoms[0]=initial.D0;atoms[1]=initial.D1;atoms[2]=initial.D2;
  for(size_t k=0;k<=c.d;++k)atoms[px0+k]=p.pairing_x[k];
  for(size_t k=0;k<c.d;++k){atoms[delta10+k]=p.delta1R[k];atoms[delta20+k]=p.delta2R[k];const auto&f=proof.dory_folds[k];atoms[fold0+6*k]=f.D1L;atoms[fold0+6*k+1]=f.D1R;atoms[fold0+6*k+2]=f.D2L;atoms[fold0+6*k+3]=f.D2R;atoms[fold0+6*k+4]=f.W1;atoms[fold0+6*k+5]=f.W2;atoms[u0+k]=proof.batch_U[k];atoms[fresh0+3*k]=fresh[k+1].D0;atoms[fresh0+3*k+1]=fresh[k+1].D1;atoms[fresh0+3*k+2]=fresh[k+1].D2;}
  using Expr=std::vector<Fr>;Fr one;one=1;auto zero=[&](){Expr x(atom_count);for(auto&v:x)v.clear();return x;};auto atom=[&](Expr&x,size_t i,const Fr&v){Fr::add(x[i],x[i],v);};auto scaled=[&](Expr&dst,const Expr&src,const Fr&v){for(size_t i=0;i<atom_count;++i){Fr z;Fr::mul(z,src[i],v);Fr::add(dst[i],dst[i],z);}};
  std::array<Expr,3>agg{zero(),zero(),zero()};atom(agg[0],initial0,one);atom(agg[1],initial0+1,one);atom(agg[2],initial0+2,one);
  tr.beta.resize(c.d);tr.alpha.resize(c.d);tr.gamma.resize(c.d+1);tr.gamma[0].clear();
  for(size_t t=1;t<=c.d;++t){
   size_t k=t-1,m=c.n>>k,h=m/2;const auto&f=proof.dory_folds[k];std::array<Bytes,6>bm{encode_u64_be(k),encode_u64_be(m),serialize(f.D1L),serialize(f.D1R),serialize(f.D2L),serialize(f.D2R)};
   a=Clock::now();T.absorb("VPIP.BF/DORY-BETA-MESSAGE/V2",bm);b=Clock::now();Fr beta=T.challenge_nonzero("VPIP.BF/BETA/V2",k);auto e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr bi=inverse_nonzero(beta);tr.beta[k]=beta;
   std::array<Bytes,4>am{encode_u64_be(k),encode_u64_be(m),serialize(f.W1),serialize(f.W2)};a=Clock::now();T.absorb("VPIP.BF/DORY-ALPHA-MESSAGE/V2",am);b=Clock::now();Fr alpha=T.challenge_nonzero("VPIP.BF/ALPHA/V2",k);e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ai=inverse_nonzero(alpha);tr.alpha[k]=alpha;
   a=Clock::now();std::array<Expr,3>bar{zero(),zero(),zero()};scaled(bar[0],agg[0],one);atom(bar[0],px0+k,one);scaled(bar[0],agg[1],bi);scaled(bar[0],agg[2],beta);atom(bar[0],fold0+6*k+4,alpha);atom(bar[0],fold0+6*k+5,ai);
   atom(bar[1],fold0+6*k,alpha);atom(bar[1],fold0+6*k+1,one);atom(bar[1],px0+k+1,fm(alpha,beta));atom(bar[1],delta10+k,beta);
   atom(bar[2],fold0+6*k+2,ai);atom(bar[2],fold0+6*k+3,one);atom(bar[2],px0+k+1,fm(ai,bi));atom(bar[2],delta20+k,bi);
   b=Clock::now();if(timing)timing->dory_recurrence_ms+=elapsed(a,b);
   const GT&u=proof.batch_U[k];std::array<Bytes,3>um{encode_u64_be(t),encode_u64_be(h),serialize(u)};a=Clock::now();T.absorb("VPIP.BF/BATCH-U/V2",um);b=Clock::now();Fr gamma=T.challenge_nonzero("VPIP.BF/GAMMA/V2",t);e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}tr.gamma[t]=gamma;Fr gs=fm(gamma,gamma);
   a=Clock::now();std::array<Expr,3>next{zero(),zero(),zero()};scaled(next[0],bar[0],gs);atom(next[0],u0+k,gamma);atom(next[0],fresh0+3*k,one);scaled(next[1],bar[1],gamma);atom(next[1],fresh0+3*k+1,one);scaled(next[2],bar[2],gamma);atom(next[2],fresh0+3*k+2,one);agg=std::move(next);b=Clock::now();if(timing)timing->dory_recurrence_ms+=elapsed(a,b);
  }
  std::array<Bytes,2>fin{serialize(proof.PhiFinal),serialize(proof.ThetaFinal)};a=Clock::now();T.absorb("VPIP.BF/DORY-FINAL/V2",fin);b=Clock::now();Fr epsilon=T.challenge_nonzero("VPIP.BF/EPSILON/V2",c.d);auto e=Clock::now();if(timing){timing->transcript_replay_ms+=elapsed(a,b);timing->challenge_derivation_ms+=elapsed(b,e);}Fr ei=inverse_nonzero(epsilon);tr.epsilon=epsilon;
  a=Clock::now();Expr lhs_expr=zero();scaled(lhs_expr,agg[0],one);scaled(lhs_expr,agg[1],ei);scaled(lhs_expr,agg[2],epsilon);atom(lhs_expr,px0+c.d,one);GT lhs=gt_multiexp_pippenger(atoms,lhs_expr);b=Clock::now();if(timing)timing->dory_recurrence_ms+=elapsed(a,b);
  a=Clock::now();G1 g1=g1_add(proof.PhiFinal,g1_pow(c.G[0],epsilon));G2 g2=g2_add(proof.ThetaFinal,g2_pow(c.H[0],ei));GT rhs;mcl::bn::pairing(rhs,g1,g2);GT rhs_inv;GT::inv(rhs_inv,rhs);tr.dory_residual=gt_mul(lhs,rhs_inv);tr.dory_accepted=lhs==rhs;b=Clock::now();if(timing)timing->terminal_dory_ms+=elapsed(a,b);
  a=Clock::now();GT rexp;mcl::bn::pairing(rexp,proof.R,c.H[0]);GT outer_inv;GT::inv(outer_inv,outer);tr.rexp_residual=gt_mul(rexp,outer_inv);tr.rexp_accepted=rexp==outer;b=Clock::now();if(timing)timing->terminal_rexp_ms+=elapsed(a,b);
  tr.final_rexp_d1=outer;tr.accepted=tr.dory_accepted&&tr.rexp_accepted;
  if(timing){timing->online_verify_ms=elapsed(total_start,Clock::now());double known=timing->transcript_init_ms+timing->transcript_replay_ms+timing->challenge_derivation_ms+timing->tensor_reconstruction_ms+timing->g1_msm_y_ms+timing->dory_recurrence_ms+timing->gt_multiexp_ms+timing->terminal_dory_ms+timing->terminal_rexp_ms;timing->other_online_ms=timing->online_verify_ms-known;}
  return tr;
 }catch(...){if(timing)timing->online_verify_ms=elapsed(total_start,Clock::now());return {};}
}

ReferenceVerificationTrace verify_reference_diagnostic(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&proof){if(!validate_verification_inputs(c,p,s,proof))return {};return verify_core_unchecked(c,p,s,proof);}
bool verify_reference(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,const VpipBfProof&proof){return verify_reference_diagnostic(c,p,s,proof).accepted;}
}

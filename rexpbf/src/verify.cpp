#include "rexpbf/verify.hpp"
#include "rexpbf/gt_multiexp.hpp"
#include "rexpbf/pairing.hpp"
#include "internal/crypto.hpp"
#include "internal/proof.hpp"
#include "internal/symbolic.hpp"
#include <chrono>
#include <stdexcept>

namespace rexpbf {
namespace {
using Clock=std::chrono::steady_clock;
double millis(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
Fr one(){return internal::fr_one();}
Fr inv(const Fr&x){return internal::fr_inverse_nonzero(x);}
G1 g1mul(const G1&p,const Fr&x){return internal::g1_multiply(p,x);} G2 g2mul(const G2&p,const Fr&x){return internal::g2_multiply(p,x);}
G1 g1add(const G1&a,const G1&b){return internal::g1_add(a,b);} G2 g2add(const G2&a,const G2&b){return internal::g2_add(a,b);}
}
std::vector<Fr> batch_invert_nonzero(std::span<const Fr> v){
    initialize_bn254();if(v.empty())return{};std::vector<Fr> prefix(v.size()),out(v.size());Fr acc=one();
    for(std::size_t i=0;i<v.size();++i){if(v[i].isZero())throw std::invalid_argument("cannot invert zero");Fr::mul(acc,acc,v[i]);prefix[i]=acc;}
    Fr back=inv(acc);for(std::size_t i=v.size();i-->0;){Fr before=i?prefix[i-1]:one();Fr::mul(out[i],back,before);Fr::mul(back,back,v[i]);}return out;
}
template<bool Profile,bool Diagnostics>
static bool core_impl(const ValidatedVerificationInputs& token,VerifyDiagnostics* diagnostics,VerifyCoreBreakdown* profile){
 const auto&c=token.crs();const auto&p=token.precomputation();const auto&s=token.statement();const auto&q=token.proof();
 Clock::time_point total,phase;
 if constexpr(Profile){total=Clock::now();phase=total;}
 try{
  ChallengeTrace ch;
  if constexpr(Profile){TranscriptMetrics tm;ch=internal::replay_challenges_prevalidated(c,s,q,&tm);auto now=Clock::now();profile->transcript_replay_ms=millis(phase,now);phase=now;profile->transcript_serialization_ms=tm.serialization_ms;profile->sha256_ms=tm.sha256_ms;profile->challenge_to_field_ms=tm.challenge_to_field_ms;profile->transcript_entries=tm.transcript_entries;profile->transcript_bytes_absorbed=tm.bytes_absorbed;profile->sha256_calls=tm.sha256_calls;profile->challenge_derivations=tm.challenge_derivations;profile->rejection_sampling_retries=tm.rejection_sampling_retries;}
  else ch=internal::replay_challenges_prevalidated(c,s,q);
  std::vector<Fr> all;all.reserve(ch.rho.size()+ch.beta.size()+ch.alpha.size()+1);all.insert(all.end(),ch.rho.begin(),ch.rho.end());all.insert(all.end(),ch.beta.begin(),ch.beta.end());all.insert(all.end(),ch.alpha.begin(),ch.alpha.end());all.push_back(ch.epsilon);
  auto inverses=batch_invert_nonzero(all);
  if constexpr(Profile){auto now=Clock::now();profile->batch_inversion_ms=millis(phase,now);phase=now;}
  std::size_t pos=0;const auto inverse_span=std::span<const Fr>(inverses);auto ri=inverse_span.subspan(pos,c.d);pos+=c.d;
  auto bi=inverse_span.subspan(pos,c.d-1);pos+=c.d-1;auto ai=inverse_span.subspan(pos,c.d-1);Fr ei=inverses.back();
  internal::SymbolicMetrics symbolic_metrics;
  auto ex=internal::build_symbolic_expressions(
      p,s,q,ch,ri,bi,ai,Profile ? &symbolic_metrics : nullptr);
  internal::GTExpression expr=ex.accumulated;auto z=ex.left;z.scale(ei);expr.multiply(z);z=ex.right;z.scale(ch.epsilon);expr.multiply(z);expr.multiply_atom(p.x[c.d],one());
  if constexpr(Profile){profile->symbolic_initialization_ms=symbolic_metrics.initialization_ms;profile->symbolic_dory_fold_ms=symbolic_metrics.dory_fold_ms;profile->symbolic_rexp_fresh_ms=symbolic_metrics.fresh_rexp_ms;profile->symbolic_batch_ms=symbolic_metrics.batch_ms;profile->normalization_ms=symbolic_metrics.normalization_ms;phase=Clock::now();profile->fresh_rexp_instances=c.d;profile->dory_fold_steps=c.d-1;profile->batching_u_messages=q.steps.size();profile->gamma_challenges=ch.gamma.size();profile->rexp_challenges=ch.rho.size();profile->dory_fold_challenges=ch.beta.size()+ch.alpha.size();profile->eta_challenges=0;profile->final_dory_checks=1;profile->final_rexp_checks=1;profile->combined_terminal_checks=0;}
  const std::size_t expr_terms_before=expr.terms().size();
  const std::size_t outer_terms_before=ex.outer.terms().size();
  const std::size_t expr_zeros=expr.normalize();
  const std::size_t outer_zeros=ex.outer.normalize();
  const std::size_t expr_terms_after=expr.terms().size();
  const std::size_t outer_terms_after=ex.outer.terms().size();
  const std::size_t terms_before=expr_terms_before+outer_terms_before;
  const std::size_t zeros=expr_zeros+outer_zeros;
  const std::size_t terms_after=expr_terms_after+outer_terms_after;
  if constexpr(Profile){auto now=Clock::now();profile->normalization_ms+=millis(phase,now);phase=now;profile->normalization_calls=symbolic_metrics.normalization_calls+2;profile->total_terms_before_normalization=symbolic_metrics.terms_before_normalization+terms_before;profile->total_terms_after_normalization=symbolic_metrics.terms_after_normalization+terms_after;profile->zero_terms_removed=symbolic_metrics.zero_terms_removed+zeros;profile->duplicate_terms_coalesced=profile->total_terms_before_normalization-profile->total_terms_after_normalization-profile->zero_terms_removed;}
  std::vector<GT>b,outer_b;std::vector<Fr>sc,outer_sc;b.reserve(expr.terms().size());sc.reserve(expr.terms().size());outer_b.reserve(ex.outer.terms().size());outer_sc.reserve(ex.outer.terms().size());
  for(auto&t:expr.terms()){b.push_back(*t.base);sc.push_back(t.scalar);}for(auto&t:ex.outer.terms()){outer_b.push_back(*t.base);outer_sc.push_back(t.scalar);}
  GT lhs=gt_multiexp_pippenger(b,sc);GT fresh_d1=gt_multiexp_pippenger(outer_b,outer_sc);
  if constexpr(Profile){auto now=Clock::now();profile->gt_multiexp_ms=millis(phase,now);phase=now;}
  G1 A=g1add(q.phi_final,g1mul(c.gamma[0],ch.epsilon));G2 B=g2add(q.theta_final,g2mul(c.lambda[0],ei));
  if constexpr(Profile){auto now=Clock::now();profile->final_dory_g1_g2_ms=millis(phase,now);profile->terminal_g1_g2_ms=profile->final_dory_g1_g2_ms;phase=now;}
  GT dory_pairing;mcl::bn::pairing(dory_pairing,A,B);
  if constexpr(Profile){auto now=Clock::now();profile->final_dory_pairing_ms=millis(phase,now);phase=now;}
  const G1& rexp_g1=q.r_final;const G2& rexp_g2=c.lambda[0];
  if constexpr(Profile){auto now=Clock::now();profile->final_rexp_preparation_ms=millis(phase,now);phase=now;}
  GT rexp_pairing;mcl::bn::pairing(rexp_pairing,rexp_g1,rexp_g2);
  if constexpr(Profile){auto now=Clock::now();profile->final_rexp_pairing_ms=millis(phase,now);profile->terminal_pairing_ms=profile->final_dory_pairing_ms+profile->final_rexp_pairing_ms;phase=now;}
  const bool final_dory_ok=lhs==dory_pairing;const bool final_rexp_ok=rexp_pairing==fresh_d1;
  if constexpr(Diagnostics){auto&d=*diagnostics;d={};d.gt_terms_before_coalescing=expr_terms_before;d.gt_zero_scalar_terms=expr_zeros;d.gt_terms_after_coalescing=expr_terms_after;d.gt_nonzero_scalar_terms=expr_terms_after;d.gt_multiexp_calls=2;d.terminal_pairing_terms=2;d.gt_bases=b;d.gt_scalars=sc;d.terminal_g1={A,q.r_final};d.terminal_g2={B,c.lambda[0]};d.gt_result=lhs;d.terminal_pairing_result=dory_pairing;}
  if constexpr(Profile){auto now=Clock::now();profile->identity_check_ms=millis(phase,now);profile->gt_multiexp_calls=2;profile->terminal_pairing_terms=2;profile->subgroup_validation_calls=0;profile->total_ms=millis(total,now);const double accounted=profile->transcript_replay_ms+profile->batch_inversion_ms+profile->symbolic_initialization_ms+profile->symbolic_dory_fold_ms+profile->symbolic_rexp_fresh_ms+profile->symbolic_batch_ms+profile->normalization_ms+profile->gt_multiexp_ms+profile->final_dory_g1_g2_ms+profile->final_dory_pairing_ms+profile->final_rexp_preparation_ms+profile->final_rexp_pairing_ms+profile->identity_check_ms;profile->other_profiling_overhead_ms=profile->total_ms-accounted;}
  return final_dory_ok && final_rexp_ok;
 }catch(...){return false;}
}
bool verify_prevalidated(const ValidatedVerificationInputs&t,VerifyCoreBreakdown*b){if(b){*b={};return core_impl<true,false>(t,nullptr,b);}return core_impl<false,false>(t,nullptr,nullptr);}
bool verify_online(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
 auto token=validate_verification_inputs(c,p,s,q);return token&&verify_prevalidated(*token);
}
bool verify_online_with_breakdown(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,VerifyCoreBreakdown&b){
 b={};auto token=validate_verification_inputs(c,p,s,q);return token&&verify_prevalidated(*token,&b);
}
bool verify_with_diagnostics(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,VerifyDiagnostics&d){
 auto token=validate_verification_inputs(c,p,s,q);return token&&core_impl<false,true>(*token,&d,nullptr);
}
bool verify(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
 return verify_online(c,p,s,q);
}
}

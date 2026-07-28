#include "rexpbf/verify.hpp"
#include "rexpbf/gt_multiexp.hpp"
#include "rexpbf/pairing.hpp"
#include "rexpbf/serialization.hpp"
#include <algorithm>
#include <chrono>
#include <map>
#include <limits>
#include <stdexcept>

namespace rexpbf {
namespace {
using Clock=std::chrono::steady_clock;
double millis(Clock::time_point a,Clock::time_point b){return std::chrono::duration<double,std::milli>(b-a).count();}
Fr one(){Fr x;x=1;return x;}
Fr mul(const Fr&a,const Fr&b){Fr z;Fr::mul(z,a,b);return z;} Fr inv(const Fr&x){Fr z;Fr::inv(z,x);return z;}
GT pw(const GT&a,const Fr&x){GT z;GT::pow(z,a,x);return z;} GT mm(const GT&a,const GT&b){GT z;GT::mul(z,a,b);return z;}
GT prod(std::initializer_list<GT>x){GT z;z.setOne();for(auto&a:x)z=mm(z,a);return z;}
G1 g1mul(const G1&p,const Fr&x){G1 z;G1::mul(z,p,x);return z;} G2 g2mul(const G2&p,const Fr&x){G2 z;G2::mul(z,p,x);return z;}
G1 g1add(const G1&a,const G1&b){G1 z;G1::add(z,a,b);return z;} G2 g2add(const G2&a,const G2&b){G2 z;G2::add(z,a,b);return z;}
bool valid_proof(const CRS&c,const Proof&p){
    if(p.steps.size()!=c.d-1||!p.phi_final.isValid()||!p.phi_final.isValidOrder()
       ||!p.theta_final.isValid()||!p.theta_final.isValidOrder()||!p.r_final.isValid()||!p.r_final.isValidOrder()
       ||p.r_final.isZero())return false;
    auto gt=[](const GT&x){return mcl::bn::isValidGT(x);};
    for(auto&s:p.steps)if(!gt(s.dory_fold.d1_left)||!gt(s.dory_fold.d1_right)||!gt(s.dory_fold.d2_left)
      ||!gt(s.dory_fold.d2_right)||!gt(s.dory_fold.w1)||!gt(s.dory_fold.w2)||!gt(s.rexp_round.e)
      ||!gt(s.rexp_round.f)||!gt(s.rexp_round.t_left)||!gt(s.rexp_round.t_right)||!gt(s.u))return false;
    return true;
}
bool inputs(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
    return validate_crs(c)&&validate_precomputation_shape(c,p)&&validate_statement_shape(c,s)
      &&p.crs_digest==c.digest&&s.crs_digest==c.digest&&valid_proof(c,q);
}
struct Expressions {GTExpression a0,a1,a2,outer;};
template<bool Profile>
Expressions symbolic(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,const ChallengeTrace&ch,
                     const std::vector<Fr>&ri,const std::vector<Fr>&bi,const std::vector<Fr>&ai,double* normalization_ms){
    (void)c;
    auto a0=GTExpression::atom(s.d1_initial);a0.multiply_atom(s.e0,ch.rho[0]);a0.multiply_atom(s.f0,ri[0]);
    auto a1=GTExpression::atom(s.t_left0);a1.multiply_atom(s.t_right0,ch.rho[0]);auto outer=a1;
    auto a2=GTExpression::atom(p.x[1]);a2.multiply_atom(p.delta2_right[0],ri[0]);
    for(std::size_t i=0;i<q.steps.size();++i){std::size_t t=i+2,ell=i+1;auto&z=q.steps[i];
      auto f0=a0;f0.multiply_atom(p.x[ell],one());auto e=a1;e.scale(bi[i]);f0.multiply(e);e=a2;e.scale(ch.beta[i]);f0.multiply(e);
      f0.multiply_atom(z.dory_fold.w1,ch.alpha[i]);f0.multiply_atom(z.dory_fold.w2,ai[i]);
      GTExpression f1;f1.multiply_atom(z.dory_fold.d1_left,ch.alpha[i]);f1.multiply_atom(z.dory_fold.d1_right,one());
      f1.multiply_atom(p.x[t],mul(ch.alpha[i],ch.beta[i]));f1.multiply_atom(p.delta1_right[ell],ch.beta[i]);
      GTExpression f2;f2.multiply_atom(z.dory_fold.d2_left,ai[i]);f2.multiply_atom(z.dory_fold.d2_right,one());
      f2.multiply_atom(p.x[t],mul(ai[i],bi[i]));f2.multiply_atom(p.delta2_right[ell],bi[i]);
      auto g0=outer;g0.multiply_atom(z.rexp_round.e,ch.rho[t-1]);g0.multiply_atom(z.rexp_round.f,ri[t-1]);
      GTExpression g1;g1.multiply_atom(z.rexp_round.t_left,one());g1.multiply_atom(z.rexp_round.t_right,ch.rho[t-1]);
      GTExpression g2;g2.multiply_atom(p.x[t],one());g2.multiply_atom(p.delta2_right[ell],ri[t-1]);
      Fr gamma2=mul(ch.gamma[i],ch.gamma[i]);f0.scale(gamma2);f0.multiply_atom(z.u,ch.gamma[i]);f0.multiply(g0);
      f1.scale(ch.gamma[i]);f1.multiply(g1);f2.scale(ch.gamma[i]);f2.multiply(g2);
      Clock::time_point normalization_start;
      if constexpr(Profile)normalization_start=Clock::now();
      f0.normalize();f1.normalize();f2.normalize();g1.normalize();
      if constexpr(Profile)*normalization_ms+=millis(normalization_start,Clock::now());
      a0=std::move(f0);a1=std::move(f1);a2=std::move(f2);outer=std::move(g1);
    }return{a0,a1,a2,outer};
}
}
std::optional<ValidatedVerificationInputs> validate_verification_inputs(
 const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,ValidationBreakdown* bp){
 ValidationBreakdown local;auto&b=bp?*bp:local;b={};auto total=Clock::now();auto phase=total;
 try{
  bool shape=c.d>=1&&c.d<std::numeric_limits<std::size_t>::digits&&c.n==(std::size_t{1}<<c.d)&&c.gamma.size()==c.n&&c.lambda.size()==c.n;
  b.crs_shape_ms=millis(phase,Clock::now());if(!shape)return std::nullopt;
  phase=Clock::now();bool digest=c.digest==compute_crs_digest(c);b.crs_digest_ms=millis(phase,Clock::now());if(!digest)return std::nullopt;
  phase=Clock::now();for(auto&x:c.gamma){++b.g1_elements_checked;++b.crs_g1_checked;if(!x.isValid()||!x.isValidOrder()||x.isZero())return std::nullopt;}
  for(auto&x:c.lambda){++b.g2_elements_checked;++b.crs_g2_checked;if(!x.isValid()||!x.isValidOrder()||x.isZero())return std::nullopt;}b.crs_group_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();shape=p.x.size()==c.d+1&&p.delta1_right.size()==c.d&&p.delta2_right.size()==c.d&&p.pairing_terms==4*std::uint64_t(c.n)-3;b.precomputation_shape_ms=millis(phase,Clock::now());if(!shape)return std::nullopt;
  phase=Clock::now();bool bind=p.crs_digest==c.digest;b.precomputation_digest_binding_ms=millis(phase,Clock::now());if(!bind)return std::nullopt;
  auto check_gt=[&](const GT&x){auto a=Clock::now();bool ok=mcl::bn::isValidGT(x);b.gt_subgroup_validation_ms+=millis(a,Clock::now());++b.gt_elements_checked;++b.gt_subgroup_checks;return ok;};
  phase=Clock::now();for(auto&v:p.x){++b.precomputation_gt_checked;if(!check_gt(v))return std::nullopt;}for(auto&v:p.delta1_right){++b.precomputation_gt_checked;if(!check_gt(v))return std::nullopt;}for(auto&v:p.delta2_right){++b.precomputation_gt_checked;if(!check_gt(v))return std::nullopt;}b.precomputation_gt_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();shape=s.crs_digest==c.digest&&s.h.size()==c.n&&s.pairing_terms==3*std::uint64_t(c.n);b.statement_shape_ms=millis(phase,Clock::now());if(!shape)return std::nullopt;
  phase=Clock::now();digest=s.digest==compute_statement_digest(c,s);b.statement_digest_ms=millis(phase,Clock::now());if(!digest)return std::nullopt;
  phase=Clock::now();for(auto&x:s.h){++b.g1_elements_checked;++b.statement_g1_checked;if(!x.isValid()||!x.isValidOrder()||x.isZero())return std::nullopt;}b.statement_g1_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();for(auto*x:{&s.d1_initial,&s.e0,&s.f0,&s.t_left0,&s.t_right0}){++b.statement_gt_checked;if(!check_gt(*x))return std::nullopt;}b.statement_gt_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();shape=q.steps.size()==c.d-1;b.proof_shape_ms=millis(phase,Clock::now());if(!shape)return std::nullopt;
  phase=Clock::now();for(auto*x:{&q.phi_final,&q.r_final}){++b.g1_elements_checked;++b.proof_g1_checked;if(!x->isValid()||!x->isValidOrder()||x->isZero())return std::nullopt;}b.proof_g1_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();++b.g2_elements_checked;++b.proof_g2_checked;if(!q.theta_final.isValid()||!q.theta_final.isValidOrder()||q.theta_final.isZero())return std::nullopt;b.proof_g2_validation_ms=millis(phase,Clock::now());
  phase=Clock::now();for(auto&z:q.steps)for(auto*x:{&z.dory_fold.d1_left,&z.dory_fold.d1_right,&z.dory_fold.d2_left,&z.dory_fold.d2_right,&z.dory_fold.w1,&z.dory_fold.w2,&z.rexp_round.e,&z.rexp_round.f,&z.rexp_round.t_left,&z.rexp_round.t_right,&z.u}){++b.proof_gt_checked;if(!check_gt(*x))return std::nullopt;}b.proof_gt_validation_ms=millis(phase,Clock::now());
  b.total_ms=millis(total,Clock::now());return ValidatedVerificationInputs(c,p,s,q);
 }catch(...){b.total_ms=millis(total,Clock::now());return std::nullopt;}
}
GTExpression GTExpression::atom(const GT&b){GTExpression e;e.multiply_atom(b,one());return e;}
void GTExpression::multiply(const GTExpression&o){terms_.insert(terms_.end(),o.terms_.begin(),o.terms_.end());}
void GTExpression::multiply_atom(const GT&b,const Fr&s){terms_.push_back({&b,s});}
void GTExpression::scale(const Fr&s){for(auto&t:terms_)Fr::mul(t.scalar,t.scalar,s);}
std::size_t GTExpression::normalize(){std::map<const GT*,Fr> m;for(auto&t:terms_){auto it=m.find(t.base);if(it==m.end())m.emplace(t.base,t.scalar);else Fr::add(it->second,it->second,t.scalar);}
  std::size_t zeros=0;terms_.clear();for(auto&[b,s]:m){if(s.isZero())++zeros;else terms_.push_back({b,s});}return zeros;}
std::vector<Fr> batch_invert_nonzero(std::span<const Fr> v){
    initialize_bn254();if(v.empty())return{};std::vector<Fr> prefix(v.size()),out(v.size());Fr acc=one();
    for(std::size_t i=0;i<v.size();++i){if(v[i].isZero())throw std::invalid_argument("cannot invert zero");Fr::mul(acc,acc,v[i]);prefix[i]=acc;}
    Fr back=inv(acc);for(std::size_t i=v.size();i-->0;){Fr before=i?prefix[i-1]:one();Fr::mul(out[i],back,before);Fr::mul(back,back,v[i]);}return out;
}
bool verify_reference(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
  try{if(!inputs(c,p,s,q)||c.lambda[0].isZero())return false;auto ch=replay_challenges(c,s,q);Fr ri=inv(ch.rho[0]);
   GT a0=prod({s.d1_initial,pw(s.e0,ch.rho[0]),pw(s.f0,ri)});
   GT a1=mm(s.t_left0,pw(s.t_right0,ch.rho[0]));GT outer=a1;GT a2=mm(p.x[1],pw(p.delta2_right[0],ri));
  for(std::size_t i=0;i<q.steps.size();++i){auto&z=q.steps[i];std::size_t t=i+2,ell=i+1;Fr bi=inv(ch.beta[i]),ai=inv(ch.alpha[i]);ri=inv(ch.rho[t-1]);
   GT f0=prod({a0,p.x[ell],pw(a1,bi),pw(a2,ch.beta[i]),pw(z.dory_fold.w1,ch.alpha[i]),pw(z.dory_fold.w2,ai)});
   GT f1=prod({pw(z.dory_fold.d1_left,ch.alpha[i]),z.dory_fold.d1_right,pw(p.x[t],mul(ch.alpha[i],ch.beta[i])),pw(p.delta1_right[ell],ch.beta[i])});
   GT f2=prod({pw(z.dory_fold.d2_left,ai),z.dory_fold.d2_right,pw(p.x[t],mul(ai,bi)),pw(p.delta2_right[ell],bi)});
   GT g0=prod({outer,pw(z.rexp_round.e,ch.rho[t-1]),pw(z.rexp_round.f,ri)});
    GT g1=mm(z.rexp_round.t_left,pw(z.rexp_round.t_right,ch.rho[t-1]));GT g2=mm(p.x[t],pw(p.delta2_right[ell],ri));
    Fr g2s=mul(ch.gamma[i],ch.gamma[i]);a0=prod({pw(f0,g2s),pw(z.u,ch.gamma[i]),g0});a1=mm(pw(f1,ch.gamma[i]),g1);a2=mm(pw(f2,ch.gamma[i]),g2);outer=g1;}
  Fr ei=inv(ch.epsilon);GT lhs=prod({a0,pw(a1,ei),pw(a2,ch.epsilon),p.x[c.d]});
  G1 A=g1add(q.phi_final,g1mul(c.gamma[0],ch.epsilon));G2 B=g2add(q.theta_final,g2mul(c.lambda[0],ei));
   GT pd,pr;mcl::bn::pairing(pd,A,B);mcl::bn::pairing(pr,q.r_final,c.lambda[0]);return lhs==pd&&pr==outer;
 }catch(...){return false;}
}
template<bool Profile,bool Diagnostics>
static bool core_impl(const ValidatedVerificationInputs& token,VerifyDiagnostics* diagnostics,VerifyCoreBreakdown* profile){
 const auto&c=token.crs();const auto&p=token.precomputation();const auto&s=token.statement();const auto&q=token.proof();
 Clock::time_point total,phase;
 if constexpr(Profile){total=Clock::now();phase=total;}
 try{
  ChallengeTrace ch;
  if constexpr(Profile){TranscriptMetrics tm;ch=replay_challenges(c,s,q,&tm,true);auto now=Clock::now();profile->transcript_replay_ms=millis(phase,now);phase=now;profile->transcript_serialization_ms=tm.serialization_ms;profile->sha256_ms=tm.sha256_ms;profile->challenge_to_field_ms=tm.challenge_to_field_ms;profile->transcript_entries=tm.transcript_entries;profile->transcript_bytes_absorbed=tm.bytes_absorbed;profile->sha256_calls=tm.sha256_calls;profile->challenge_derivations=tm.challenge_derivations;profile->rejection_sampling_retries=tm.rejection_sampling_retries;}
  else ch=replay_challenges(c,s,q,nullptr,true);
  std::vector<Fr> all;all.insert(all.end(),ch.rho.begin(),ch.rho.end());all.insert(all.end(),ch.beta.begin(),ch.beta.end());all.insert(all.end(),ch.alpha.begin(),ch.alpha.end());all.push_back(ch.epsilon);
  auto inverses=batch_invert_nonzero(all);
  if constexpr(Profile){auto now=Clock::now();profile->batch_inversion_ms=millis(phase,now);phase=now;}
  std::size_t pos=0;std::vector<Fr> ri(inverses.begin(),inverses.begin()+c.d);pos+=c.d;
  std::vector<Fr> bi(inverses.begin()+pos,inverses.begin()+pos+c.d-1);pos+=c.d-1;std::vector<Fr> ai(inverses.begin()+pos,inverses.begin()+pos+c.d-1);Fr ei=inverses.back();
  double recurrence_normalization=0;auto ex=symbolic<Profile>(c,p,s,q,ch,ri,bi,ai,&recurrence_normalization);
  GTExpression expr=ex.a0;auto z=ex.a1;z.scale(ei);expr.multiply(z);z=ex.a2;z.scale(ch.epsilon);expr.multiply(z);expr.multiply_atom(p.x[c.d],one());
  if constexpr(Profile){auto now=Clock::now();profile->symbolic_dory_fold_ms=millis(phase,now)-recurrence_normalization;profile->normalization_ms=recurrence_normalization;phase=now;profile->fresh_rexp_instances=c.d;profile->dory_fold_steps=c.d-1;profile->batching_u_messages=q.steps.size();profile->gamma_challenges=ch.gamma.size();profile->rexp_challenges=ch.rho.size();profile->dory_fold_challenges=ch.beta.size()+ch.alpha.size();profile->eta_challenges=0;profile->final_dory_checks=1;profile->final_rexp_checks=1;profile->combined_terminal_checks=0;}
  const std::size_t terms_before=expr.terms().size()+ex.outer.terms().size();
  const std::size_t zeros=expr.normalize()+ex.outer.normalize();
  const std::size_t terms_after=expr.terms().size()+ex.outer.terms().size();
  if constexpr(Profile){auto now=Clock::now();profile->normalization_ms+=millis(phase,now);phase=now;profile->normalization_calls=4*(c.d-1)+2;profile->total_terms_before_normalization=terms_before;profile->total_terms_after_normalization=terms_after;profile->zero_terms_removed=zeros;profile->duplicate_terms_coalesced=terms_before-terms_after-zeros;}
  std::vector<GT>b,outer_b;std::vector<Fr>sc,outer_sc;
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
  if constexpr(Diagnostics){auto&d=*diagnostics;d={};d.gt_terms_before_coalescing=terms_before;d.gt_zero_scalar_terms=zeros;d.gt_terms_after_coalescing=terms_after+zeros;d.gt_nonzero_scalar_terms=terms_after;d.gt_multiexp_calls=2;d.terminal_pairing_terms=2;d.gt_bases=b;d.gt_scalars=sc;d.terminal_g1={A,q.r_final};d.terminal_g2={B,c.lambda[0]};d.gt_result=lhs;d.terminal_pairing_result=dory_pairing;}
  if constexpr(Profile){auto now=Clock::now();profile->identity_check_ms=millis(phase,now);profile->gt_multiexp_calls=2;profile->terminal_pairing_terms=2;profile->subgroup_validation_calls=0;profile->total_ms=millis(total,now);const double accounted=profile->transcript_replay_ms+profile->batch_inversion_ms+profile->symbolic_dory_fold_ms+profile->normalization_ms+profile->gt_multiexp_ms+profile->final_dory_g1_g2_ms+profile->final_dory_pairing_ms+profile->final_rexp_preparation_ms+profile->final_rexp_pairing_ms+profile->identity_check_ms;profile->other_profiling_overhead_ms=profile->total_ms-accounted;}
  return final_dory_ok && final_rexp_ok;
 }catch(...){return false;}
}
bool verify_prevalidated(const ValidatedVerificationInputs&t,VerifyCoreBreakdown*b){if(b){*b={};return core_impl<true,false>(t,nullptr,b);}return core_impl<false,false>(t,nullptr,nullptr);}
bool verify_online(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
 return core_impl<false,false>(ValidatedVerificationInputs(c,p,s,q),nullptr,nullptr);
}
bool verify_online_with_breakdown(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,VerifyCoreBreakdown&b){
 b={};return core_impl<true,false>(ValidatedVerificationInputs(c,p,s,q),nullptr,&b);
}
bool verify_with_diagnostics(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q,VerifyDiagnostics&d){
 auto token=validate_verification_inputs(c,p,s,q);return token&&core_impl<false,true>(*token,&d,nullptr);
}
bool verify(const CRS&c,const Precomputation&p,const Statement&s,const Proof&q){
 auto token=validate_verification_inputs(c,p,s,q);return token&&verify_prevalidated(*token);
}
}

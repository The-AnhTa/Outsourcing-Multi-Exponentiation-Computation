#include "vme_ibf/verify_reference.hpp"
#include "vme_ibf/serialization.hpp"
#include "vme_ibf/setup.hpp"
#include "vme_ibf/transcript.hpp"
#include <mcl/gmp_util.hpp>
#include <array>
#include <limits>
#include <stdexcept>

namespace vme_ibf {
namespace {
Fr fm(const Fr&a,const Fr&b){Fr z;Fr::mul(z,a,b);return z;}
GT prod(std::initializer_list<GT>x){GT z;z.setOne();for(const auto&q:x)z=gt_mul(z,q);return z;}
template<class T>bool canonical(const T&x){try{auto b=serialize(x);T y;return y.deserialize(b.data(),b.size())==b.size()&&y==x;}catch(...){return false;}}
bool valid_gt(const GT&x){


 if(!canonical(x))return false;mpz_class order;mcl::gmp::setStr(order,Fr::getModulo(),10);GT acc;GT::pow(acc,x,order);GT one;one.setOne();return acc==one;
}
Digest statement_input_digest(const VmeIbfCRS&c,std::span<const Fr>x){return compute_statement_input_digest(c,x);}
Digest statement_digest(const VmeIbfCRS&c,const Digest&id,const G2&X){Bytes b;append_frame(b,"VME.BF.G2/STATEMENT/V2");append_frame(b,c.digest);append_frame(b,id);append_frame(b,serialize(X));return sha256(b);}
bool valid_inputs(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,const VmeIbfStatement&s,const VmeIbfProof&proof){
 if(c.d<1||c.d>=std::numeric_limits<size_t>::digits||c.n!=(size_t{1}<<c.d)||c.G.size()!=c.n||c.H.size()!=c.n||s.x.size()!=c.n)return false;
 if(p.pairing_x.size()!=c.d+1||p.delta1R.size()!=c.d||p.delta2R.size()!=c.d||proof.rexp_claims.size()!=c.d-1||proof.dory_folds.size()!=c.d||proof.batch_U.size()!=c.d)return false;
 for(const auto&x:c.G)if(!valid_g1(x,true))return false;for(const auto&x:c.H)if(!valid_g2(x,true))return false;if(!valid_g1(c.L,true)||!valid_g2(c.Lprime,true)||!valid_g2(s.X)||!valid_g1(proof.R)||!valid_g1(proof.PhiFinal)||!valid_g2(proof.ThetaFinal))return false;
 for(const auto&x:p.pairing_x)if(!valid_gt(x))return false;for(const auto&x:p.delta1R)if(!valid_gt(x))return false;for(const auto&x:p.delta2R)if(!valid_gt(x))return false;if(!valid_gt(p.pairing_LLprime))return false;
 for(const auto&r:proof.rexp_claims)if(!valid_gt(r.E)||!valid_gt(r.F)||!valid_gt(r.TL)||!valid_gt(r.TR))return false;for(const auto&f:proof.dory_folds)if(!valid_gt(f.D1L)||!valid_gt(f.D1R)||!valid_gt(f.D2L)||!valid_gt(f.D2R)||!valid_gt(f.W1)||!valid_gt(f.W2))return false;for(const auto&u:proof.batch_U)if(!valid_gt(u))return false;
 if(compute_crs_digest(c)!=c.digest)return false;Digest id=statement_input_digest(c,s.x);if(statement_digest(c,id,s.X)!=s.digest)return false;return validate_precomputation(c,p);
}
void absorb_rexp(Transcript&T,size_t j,size_t m,const RexpClaims&r){std::array<Bytes,6>f{encode_u64_be(j),encode_u64_be(m),serialize(r.E),serialize(r.F),serialize(r.TL),serialize(r.TR)};T.absorb("VME.BF.G2/REXP-G1-CLAIMS/V2",f);}
}

ReferenceVerificationTrace verify_reference_diagnostic(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,const VmeIbfStatement&s,const VmeIbfProof&proof){
 ReferenceVerificationTrace tr;try{if(!valid_inputs(c,p,s,proof))return tr;Transcript T(s.digest);tr.rho.resize(c.d);std::vector<DoryTargetState>fresh(c.d+1);GT outer=p.pairing_x[0];
 for(size_t j=0;j<c.d;++j){RexpClaims r=j==0?RexpClaims{p.delta1R[0],p.delta2R[0],p.pairing_x[1],p.delta1R[0]}:proof.rexp_claims[j-1];absorb_rexp(T,j,c.n>>j,r);Fr rho=T.challenge_nonzero("VME.BF.G2/RHO/V2",j),ri=inverse_nonzero(rho);tr.rho[j]=rho;size_t t=j+1;fresh[t].D0=prod({outer,gt_pow(r.E,rho),gt_pow(r.F,ri)});fresh[t].D1=gt_mul(r.TL,gt_pow(r.TR,rho));fresh[t].D2=gt_mul(p.pairing_x[t],gt_pow(p.delta2R[j],ri));outer=fresh[t].D1;}
 T.absorb("VME.BF.G2/R-G1/V2",serialize(proof.R));std::vector<Fr>r(c.d);for(size_t j=0;j<c.d;++j)r[c.d-j-1]=tr.rho[j];auto weights=tensor_vector(r);if(weights.size()!=c.n)return ReferenceVerificationTrace{};Fr q=inner_product(weights,s.x);std::array<Bytes,3>ini{encode_u64_be(c.d),encode_u64_be(c.n),serialize(q)};T.absorb("VME.BF.G2/VME-INITIAL/V2",ini);DoryTargetState agg;agg.D0=gt_pow(p.pairing_LLprime,q);mcl::bn::pairing(agg.D1,c.L,s.X);mcl::bn::pairing(agg.D2,proof.R,c.Lprime);tr.beta.resize(c.d);tr.alpha.resize(c.d);tr.gamma.resize(c.d+1);tr.gamma[0].clear();
 for(size_t t=1;t<=c.d;++t){size_t k=t-1,m=c.n>>k,h=m/2;const auto&f=proof.dory_folds[k];std::array<Bytes,6>bm{encode_u64_be(k),encode_u64_be(m),serialize(f.D1L),serialize(f.D1R),serialize(f.D2L),serialize(f.D2R)};T.absorb("VME.BF.G2/DORY-BETA-MESSAGE/V2",bm);Fr beta=T.challenge_nonzero("VME.BF.G2/BETA/V2",k),bi=inverse_nonzero(beta);tr.beta[k]=beta;std::array<Bytes,4>am{encode_u64_be(k),encode_u64_be(m),serialize(f.W1),serialize(f.W2)};T.absorb("VME.BF.G2/DORY-ALPHA-MESSAGE/V2",am);Fr alpha=T.challenge_nonzero("VME.BF.G2/ALPHA/V2",k),ai=inverse_nonzero(alpha);tr.alpha[k]=alpha;DoryTargetState bar;bar.D0=prod({agg.D0,p.pairing_x[k],gt_pow(agg.D1,bi),gt_pow(agg.D2,beta),gt_pow(f.W1,alpha),gt_pow(f.W2,ai)});bar.D1=prod({gt_pow(f.D1L,alpha),f.D1R,gt_pow(p.pairing_x[k+1],fm(alpha,beta)),gt_pow(p.delta1R[k],beta)});bar.D2=prod({gt_pow(f.D2L,ai),f.D2R,gt_pow(p.pairing_x[k+1],fm(ai,bi)),gt_pow(p.delta2R[k],bi)});const GT&u=proof.batch_U[k];std::array<Bytes,3>um{encode_u64_be(t),encode_u64_be(h),serialize(u)};T.absorb("VME.BF.G2/BATCH-U/V2",um);Fr gamma=T.challenge_nonzero("VME.BF.G2/GAMMA/V2",t);tr.gamma[t]=gamma;Fr gs=fm(gamma,gamma);agg.D0=prod({gt_pow(bar.D0,gs),gt_pow(u,gamma),fresh[t].D0});agg.D1=gt_mul(gt_pow(bar.D1,gamma),fresh[t].D1);agg.D2=gt_mul(gt_pow(bar.D2,gamma),fresh[t].D2);}
 std::array<Bytes,2>fin{serialize(proof.PhiFinal),serialize(proof.ThetaFinal)};T.absorb("VME.BF.G2/DORY-FINAL/V2",fin);Fr epsilon=T.challenge_nonzero("VME.BF.G2/EPSILON/V2",c.d),ei=inverse_nonzero(epsilon);tr.epsilon=epsilon;GT lhs=prod({agg.D0,gt_pow(agg.D1,ei),gt_pow(agg.D2,epsilon),p.pairing_x[c.d]});G1 g1=g1_add(proof.PhiFinal,g1_pow(c.G[0],epsilon));G2 g2=g2_add(proof.ThetaFinal,g2_pow(c.H[0],ei));GT rhs;mcl::bn::pairing(rhs,g1,g2);GT rhs_inv;GT::inv(rhs_inv,rhs);tr.dory_residual=gt_mul(lhs,rhs_inv);tr.dory_accepted=lhs==rhs;GT rexp;mcl::bn::pairing(rexp,proof.R,c.H[0]);GT outer_inv;GT::inv(outer_inv,outer);tr.rexp_residual=gt_mul(rexp,outer_inv);tr.rexp_accepted=rexp==outer;tr.final_aggregate=agg;tr.final_rexp_d1=outer;tr.accepted=tr.dory_accepted&&tr.rexp_accepted;return tr;}catch(...){return ReferenceVerificationTrace{};}
}
bool verify_reference(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,const VmeIbfStatement&s,const VmeIbfProof&proof){return verify_reference_diagnostic(c,p,s,proof).accepted;}
}

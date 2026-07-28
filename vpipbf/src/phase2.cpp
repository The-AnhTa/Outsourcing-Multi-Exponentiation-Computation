#include "vpip_bf/phase2.hpp"
#include <array>
#include <stdexcept>

namespace vpip_bf {
namespace {
Fr fr_mul(const Fr&a,const Fr&b){Fr z;Fr::mul(z,a,b);return z;}
GT products(std::initializer_list<GT> xs){GT z;z.setOne();for(const auto&x:xs)z=gt_mul(z,x);return z;}
void absorb_initial(Transcript&T,size_t d,size_t n,const G1&Y){std::array<Bytes,3>f{encode_u64_be(d),encode_u64_be(n),serialize(Y)};T.absorb("vpipbf/initial-targets/v1",f);}
void absorb_beta(Transcript&T,size_t k,size_t m,const DoryFoldProof&f){std::array<Bytes,6>x{encode_u64_be(k),encode_u64_be(m),serialize(f.D1L),serialize(f.D1R),serialize(f.D2L),serialize(f.D2R)};T.absorb("VPIP.BF/DORY-BETA-MESSAGE/V2",x);}
void absorb_alpha(Transcript&T,size_t k,size_t m,const DoryFoldProof&f){std::array<Bytes,4>x{encode_u64_be(k),encode_u64_be(m),serialize(f.W1),serialize(f.W2)};T.absorb("VPIP.BF/DORY-ALPHA-MESSAGE/V2",x);}
void absorb_u(Transcript&T,size_t t,size_t h,const GT&u){std::array<Bytes,3>x{encode_u64_be(t),encode_u64_be(h),serialize(u)};T.absorb("VPIP.BF/BATCH-U/V2",x);}
}

Phase2Result prove_phase2(const VpipBfCRS&c,const VpipBfPrecomputation&p,const Phase1Result&p1){
 if(c.d<1||c.n!=(size_t{1}<<c.d)||c.G.size()!=c.n||c.H.size()!=c.n||p1.statement.X.size()!=c.n||p1.r.size()!=c.d||p1.fresh.size()!=c.d+1||p.pairing_x.size()!=c.d+1||p.delta1R.size()!=c.d||p.delta2R.size()!=c.d)throw std::invalid_argument("malformed Phase-II input");
 Phase2Result out;out.proof.dory_folds.reserve(c.d);out.proof.batch_U.reserve(c.d);out.challenges.beta.resize(c.d);out.challenges.alpha.resize(c.d);out.challenges.gamma.resize(c.d+1);out.challenges.gamma[0].clear();out.folded.resize(c.d+1);out.aggregate.resize(c.d+1);
 Transcript T=Transcript::resume(p1.transcript_after_R);auto s=tensor_vector(p1.r);if(s.size()!=c.n)throw std::logic_error("tensor length mismatch");out.Y=g1_multiexp(p1.statement.X,s);absorb_initial(T,c.d,c.n,out.Y);
 DoryInstanceState agg;agg.witness.Phi=p1.statement.X;agg.witness.Theta.reserve(c.n);for(size_t i=0;i<c.n;++i)agg.witness.Theta.push_back(g2_pow(c.Lprime,s[i]));mcl::bn::pairing(agg.target.D0,out.Y,c.Lprime);agg.target.D1=p1.statement.C;mcl::bn::pairing(agg.target.D2,p1.R,c.Lprime);out.initial_instance=agg;out.aggregate[0]=agg;
 for(size_t t=1;t<=c.d;++t){size_t k=t-1,m=c.n>>k,h=m/2;if(agg.witness.Phi.size()!=m||agg.witness.Theta.size()!=m)throw std::logic_error("aggregate dimension mismatch");DoryFoldProof fp;auto phi=std::span(agg.witness.Phi);auto theta=std::span(agg.witness.Theta);auto gn=std::span(c.G).first(h);auto ln=std::span(c.H).first(h);fp.D1L=pairing_product(phi.first(h),ln);fp.D1R=pairing_product(phi.subspan(h,h),ln);fp.D2L=pairing_product(gn,theta.first(h));fp.D2R=pairing_product(gn,theta.subspan(h,h));absorb_beta(T,k,m,fp);Fr beta=T.challenge_nonzero("VPIP.BF/BETA/V2",k),bi=inverse_nonzero(beta);out.challenges.beta[k]=beta;
  std::vector<G1> pc;std::vector<G2>tc;pc.reserve(m);tc.reserve(m);for(size_t i=0;i<m;++i){pc.push_back(g1_add(agg.witness.Phi[i],g1_pow(c.G[i],beta)));tc.push_back(g2_add(agg.witness.Theta[i],g2_pow(c.H[i],bi)));}fp.W1=pairing_product(std::span(pc).first(h),std::span(tc).subspan(h,h));fp.W2=pairing_product(std::span(pc).subspan(h,h),std::span(tc).first(h));absorb_alpha(T,k,m,fp);Fr alpha=T.challenge_nonzero("VPIP.BF/ALPHA/V2",k),ai=inverse_nonzero(alpha);out.challenges.alpha[k]=alpha;
  DoryInstanceState bar;bar.witness.Phi.reserve(h);bar.witness.Theta.reserve(h);for(size_t i=0;i<h;++i){bar.witness.Phi.push_back(g1_add(g1_pow(pc[i],alpha),pc[h+i]));bar.witness.Theta.push_back(g2_add(g2_pow(tc[i],ai),tc[h+i]));}bar.target.D0=products({agg.target.D0,p.pairing_x[k],gt_pow(agg.target.D1,bi),gt_pow(agg.target.D2,beta),gt_pow(fp.W1,alpha),gt_pow(fp.W2,ai)});bar.target.D1=products({gt_pow(fp.D1L,alpha),fp.D1R,gt_pow(p.pairing_x[k+1],fr_mul(alpha,beta)),gt_pow(p.delta1R[k],beta)});bar.target.D2=products({gt_pow(fp.D2L,ai),fp.D2R,gt_pow(p.pairing_x[k+1],fr_mul(ai,bi)),gt_pow(p.delta2R[k],bi)});out.folded[t]=bar;out.proof.dory_folds.push_back(fp);
  const auto&fresh=p1.fresh[t];if(fresh.Phi.size()!=h||fresh.Theta.size()!=h)throw std::logic_error("fresh dimension mismatch");GT u=gt_mul(pairing_product(bar.witness.Phi,fresh.Theta),pairing_product(fresh.Phi,bar.witness.Theta));out.proof.batch_U.push_back(u);absorb_u(T,t,h,u);Fr gamma=T.challenge_nonzero("VPIP.BF/GAMMA/V2",t);out.challenges.gamma[t]=gamma;Fr g2=fr_mul(gamma,gamma);
  DoryInstanceState next;next.target.D0=products({gt_pow(bar.target.D0,g2),gt_pow(u,gamma),fresh.D0});next.target.D1=gt_mul(gt_pow(bar.target.D1,gamma),fresh.D1);next.target.D2=gt_mul(gt_pow(bar.target.D2,gamma),fresh.D2);next.witness.Phi.reserve(h);next.witness.Theta.reserve(h);for(size_t i=0;i<h;++i){next.witness.Phi.push_back(g1_add(g1_pow(bar.witness.Phi[i],gamma),fresh.Phi[i]));next.witness.Theta.push_back(g2_add(g2_pow(bar.witness.Theta[i],gamma),fresh.Theta[i]));}agg=next;out.aggregate[t]=next;
 }
 if(agg.witness.Phi.size()!=1||agg.witness.Theta.size()!=1)throw std::logic_error("Dory did not terminate");out.proof.PhiFinal=agg.witness.Phi[0];out.proof.ThetaFinal=agg.witness.Theta[0];std::array<Bytes,2>final{serialize(out.proof.PhiFinal),serialize(out.proof.ThetaFinal)};T.absorb("VPIP.BF/DORY-FINAL/V2",final);out.challenges.epsilon=T.challenge_nonzero("VPIP.BF/EPSILON/V2",c.d);out.final_transcript_digest=T.digest();out.final_aggregate_target=agg.target;return out;
}
}

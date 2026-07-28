#include "vpip_bf/phase1.hpp"
#include <array>
#include <stdexcept>

namespace vpip_bf {
std::vector<Fr> tensor_vector(std::span<const Fr>r){std::vector<Fr>w(1);w[0]=1;for(size_t jj=r.size();jj>0;--jj){const Fr&x=r[jj-1];std::vector<Fr>n;n.reserve(w.size()*2);for(auto&a:w){n.push_back(a);Fr q;Fr::mul(q,a,x);n.push_back(q);}w=std::move(n);}return w;}
Digest compute_statement_digest(const VpipBfCRS&c,const VpipBfPrecomputation&,const VpipBfStatementInput&i,const GT&C){Bytes b;append_frame(b,"vpipbf/statement/v1");append_frame(b,c.digest);append_frame(b,i.digest);for(auto&x:i.X)append_frame(b,serialize(x));for(auto&x:c.H)append_frame(b,serialize(x));append_frame(b,serialize(c.Lprime));append_frame(b,serialize(C));return sha256(b);}
namespace {
void absorb_claim(Transcript&T,size_t j,size_t m,const RexpClaims&c){std::array<Bytes,6>f{encode_u64_be(j),encode_u64_be(m),serialize(c.E),serialize(c.F),serialize(c.TL),serialize(c.TR)};T.absorb("vpipbf/rexp-message/v1",f);}
RexpClaims claims_for(size_t j,const VpipBfPrecomputation&p,std::span<const RexpClaims>dyn){if(j==0)return {p.delta1R[0],p.delta2R[0],p.pairing_x[1],p.delta1R[0]};return dyn[j-1];}
}
Phase1Result prove_phase1(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatementInput&i){
 if(c.d<1||c.n!=(size_t{1}<<c.d)||i.X.size()!=c.n||p.pairing_x.size()!=c.d+1||p.delta1R.size()!=c.d||p.delta2R.size()!=c.d)throw std::invalid_argument("malformed Phase-I input");
 Phase1Result o;o.statement.X=i.X;o.statement.C=pairing_product(i.X,c.H);o.statement.digest=compute_statement_digest(c,p,i,o.statement.C);Transcript T(o.statement.digest);o.rho.resize(c.d);o.fresh.resize(c.d+1);std::vector<G1>ac=c.G;GT outer=p.pairing_x[0];
 for(size_t j=0;j<c.d;++j){size_t m=c.n>>j,h=m/2,t=j+1;RexpClaims cl;if(j==0)cl=claims_for(0,p,{});else{cl.E=pairing_product(std::span(ac).subspan(h,h),std::span(c.H).first(h));cl.F=pairing_product(std::span(ac).first(h),std::span(c.H).subspan(h,h));cl.TL=pairing_product(std::span(ac).first(h),std::span(c.H).first(h));cl.TR=pairing_product(std::span(ac).subspan(h,h),std::span(c.H).first(h));o.dynamic_claims.push_back(cl);}absorb_claim(T,j,m,cl);Fr rho=T.challenge_nonzero("vpipbf/challenge/rexp-r/v1",j),ri=inverse_nonzero(rho);o.rho[j]=rho;std::vector<G1>an;std::vector<G2>theta;an.reserve(h);theta.reserve(h);for(size_t z=0;z<h;++z){an.push_back(g1_add(ac[z],g1_pow(ac[h+z],rho)));theta.push_back(g2_add(c.H[z],g2_pow(c.H[h+z],ri)));}auto&f=o.fresh[t];f.D0=gt_mul(gt_mul(outer,gt_pow(cl.E,rho)),gt_pow(cl.F,ri));f.D1=gt_mul(cl.TL,gt_pow(cl.TR,rho));f.D2=gt_mul(p.pairing_x[t],gt_pow(p.delta2R[j],ri));f.Phi=an;f.Theta=std::move(theta);ac=std::move(an);outer=f.D1;}
 o.R=ac.at(0);T.absorb("vpipbf/rexp-result-R/v1",serialize(o.R));o.transcript_after_R=T.digest();o.r.resize(c.d);for(size_t j=0;j<c.d;++j)o.r[c.d-j-1]=o.rho[j];return o;
}
std::vector<Fr> replay_rho(const VpipBfCRS&c,const VpipBfPrecomputation&p,const VpipBfStatement&s,std::span<const RexpClaims>dyn,const G1&R,Digest*after){if(dyn.size()!=c.d-1)throw std::invalid_argument("dynamic claim count");Transcript T(s.digest);std::vector<Fr>rho(c.d);for(size_t j=0;j<c.d;++j){auto cl=claims_for(j,p,dyn);absorb_claim(T,j,c.n>>j,cl);rho[j]=T.challenge_nonzero("vpipbf/challenge/rexp-r/v1",j);}T.absorb("vpipbf/rexp-result-R/v1",serialize(R));if(after)*after=T.digest();return rho;}
}

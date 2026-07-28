#include "vme_ibf/phase1.hpp"
#include <stdexcept>

namespace vme_ibf {
std::vector<Fr> tensor_vector(std::span<const Fr>r){std::vector<Fr>w(1);w[0]=1;for(size_t jj=r.size();jj>0;--jj){const Fr&x=r[jj-1];std::vector<Fr>n;n.reserve(w.size()*2);for(auto&a:w){n.push_back(a);Fr q;Fr::mul(q,a,x);n.push_back(q);}w=std::move(n);}return w;}
Digest compute_statement_digest(const VmeIbfCRS&c,const VmeIbfStatementInput&i,const G2&X){Bytes b;append_frame(b,"VME.BF.G2/STATEMENT/V2");append_frame(b,c.digest);append_frame(b,i.digest);append_frame(b,serialize(X));return sha256(b);}
static void absorb_claim(Transcript&T,size_t j,size_t m,const RexpClaims&c){std::array<Bytes,6> f{encode_u64_be(j),encode_u64_be(m),serialize(c.E),serialize(c.F),serialize(c.TL),serialize(c.TR)};T.absorb("VME.BF.G2/REXP-G1-CLAIMS/V2",f);}
static RexpClaims claims_for(size_t j,const VmeIbfPrecomputation&p,std::span<const RexpClaims>dyn){if(j==0)return {p.delta1R[0],p.delta2R[0],p.pairing_x[1],p.delta1R[0]};return dyn[j-1];}
Phase1Result prove_phase1(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,const VmeIbfStatementInput&i){
 if(c.d<1||c.n!=(size_t{1}<<c.d)||i.x.size()!=c.n||p.pairing_x.size()!=c.d+1||p.delta1R.size()!=c.d||p.delta2R.size()!=c.d)throw std::invalid_argument("malformed Phase-I input");
 G2 X=g2_multiexp_protocol(c.H,i.x);Digest statement_digest=compute_statement_digest(c,i,X);
 auto out=prove_phase1_core(c,p,i.x,X,Transcript(statement_digest));
 out.statement.digest=statement_digest;
 return out;
}
Phase1Result prove_phase1_core(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,std::span<const Fr>x,const G2&X,Transcript T){
 if(c.d<1||c.n!=(size_t{1}<<c.d)||x.size()!=c.n||p.pairing_x.size()!=c.d+1||p.delta1R.size()!=c.d||p.delta2R.size()!=c.d||!valid_g2(X))throw std::invalid_argument("malformed aggregate Phase-I input");
 Phase1Result o;o.statement.x.assign(x.begin(),x.end());o.statement.X=X;o.statement.digest=T.digest();o.rho.resize(c.d);o.fresh.resize(c.d+1);std::vector<G1>ac=c.G;GT outer=p.pairing_x[0];
 for(size_t j=0;j<c.d;++j){size_t m=c.n>>j,h=m/2,t=j+1;RexpClaims cl;if(j==0)cl=claims_for(0,p,{});else{cl.E=pairing_product(std::span(ac).subspan(h,h),std::span(c.H).first(h));cl.F=pairing_product(std::span(ac).first(h),std::span(c.H).subspan(h,h));cl.TL=pairing_product(std::span(ac).first(h),std::span(c.H).first(h));cl.TR=pairing_product(std::span(ac).subspan(h,h),std::span(c.H).first(h));o.dynamic_claims.push_back(cl);}absorb_claim(T,j,m,cl);Fr rho=T.challenge_nonzero("VME.BF.G2/RHO/V2",j),ri=inverse_nonzero(rho);o.rho[j]=rho;std::vector<G1>an;std::vector<G2>theta;an.reserve(h);theta.reserve(h);for(size_t z=0;z<h;++z){an.push_back(g1_add(ac[z],g1_pow(ac[h+z],rho)));theta.push_back(g2_add(c.H[z],g2_pow(c.H[h+z],ri)));}auto&f=o.fresh[t];f.D0=gt_mul(gt_mul(outer,gt_pow(cl.E,rho)),gt_pow(cl.F,ri));f.D1=gt_mul(cl.TL,gt_pow(cl.TR,rho));f.D2=gt_mul(p.pairing_x[t],gt_pow(p.delta2R[j],ri));f.Phi=an;f.Theta=std::move(theta);ac=std::move(an);outer=f.D1;}
 if(ac.size()!=1)throw std::logic_error("Rexp did not terminate");o.R=ac[0];T.absorb("VME.BF.G2/R-G1/V2",serialize(o.R));o.transcript_after_R=T.digest();o.r.resize(c.d);for(size_t j=0;j<c.d;++j)o.r[c.d-j-1]=o.rho[j];return o;
}
std::vector<Fr> replay_rho(const VmeIbfCRS&c,const VmeIbfPrecomputation&p,const VmeIbfStatement&s,std::span<const RexpClaims>dyn,const G1&R,Digest*after){if(dyn.size()!=c.d-1)throw std::invalid_argument("dynamic claim count");Transcript T(s.digest);std::vector<Fr>rho(c.d);for(size_t j=0;j<c.d;++j){auto cl=claims_for(j,p,dyn);absorb_claim(T,j,c.n>>j,cl);rho[j]=T.challenge_nonzero("VME.BF.G2/RHO/V2",j);}T.absorb("VME.BF.G2/R-G1/V2",serialize(R));if(after)*after=T.digest();return rho;}
}

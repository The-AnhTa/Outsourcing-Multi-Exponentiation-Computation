#include "rexpbf/prove.hpp"
#include "rexpbf/pairing.hpp"
#include "rexpbf/serialization.hpp"
#include "rexpbf/transcript.hpp"
#include "internal/crypto.hpp"
#include "internal/protocol.hpp"
#include <stdexcept>

namespace rexpbf {
namespace {
Fr inv(const Fr& x) { return internal::fr_inverse_nonzero(x); }
Fr mul(const Fr& a,const Fr& b) { return internal::fr_multiply(a,b); }
G1 scale(const G1& p,const Fr& x) { return internal::g1_multiply(p,x); }
G2 scale(const G2& p,const Fr& x) { return internal::g2_multiply(p,x); }
G1 add(const G1&a,const G1&b) { return internal::g1_add(a,b); }
G2 add(const G2&a,const G2&b) { return internal::g2_add(a,b); }
GT power(const GT&a,const Fr&x) { return internal::gt_power(a,x); }
GT product(const GT&a,const GT&b) { return internal::gt_multiply(a,b); }
GT products(std::initializer_list<GT> xs) { return internal::gt_product(xs); }
std::vector<G1> fold_h(std::span<const G1> h,const Fr& rho) {
    const auto n=h.size()/2; std::vector<G1> out; out.reserve(n);
    for(std::size_t i=0;i<n;++i) out.push_back(add(h[i],scale(h[n+i],rho))); return out;
}
std::vector<G2> fresh_theta(std::span<const G2> l,const Fr& rho_inv) {
    const auto n=l.size()/2; std::vector<G2> out; out.reserve(n);
    for(std::size_t i=0;i<n;++i) out.push_back(add(l[i],scale(l[n+i],rho_inv))); return out;
}
void validate_inputs(const CRS& c,const Precomputation&p,const Statement&s,const ProverInput&i) {
    if(!validate_crs(c)||!validate_precomputation_shape(c,p)||!validate_precomputation_elements(p)
       ||!validate_statement_shape(c,s)||!validate_statement_elements(s)||!validate_statement_digest(c,s))
        throw std::invalid_argument("invalid Setup input");
    if(p.crs_digest!=c.digest||s.crs_digest!=c.digest) throw std::invalid_argument("CRS digest mismatch");
    if(i.h!=s.h) throw std::invalid_argument("prover input differs from public statement");
}
}

ProveResult prove(const CRS& c,const Precomputation&p,const Statement&s,const ProverInput&pi) {
    validate_inputs(c,p,s,pi);
    Transcript tr(internal::transcript_domain); internal::initialize_protocol_transcript(tr,c,s);
    ProveResult result; result.r.resize(c.d); result.proof.steps.reserve(c.d-1);
    std::vector<G1> h_current=s.h;
    GT outer_d1=s.d1_initial, ad0,ad1,ad2;
    std::vector<G1> aphi; std::vector<G2> atheta;
    Fr rho=internal::derive_initial_rho(tr,c); result.r[c.d-1]=rho; Fr rho_inv=inv(rho);
    auto h_next=fold_h(h_current,rho); auto theta=fresh_theta(lambda_level(c,0),rho_inv);
    ad0=products({outer_d1,power(s.e0,rho),power(s.f0,rho_inv)});
    ad1=product(s.t_left0,power(s.t_right0,rho));
    ad2=product(p.x[1],power(p.delta2_right[0],rho_inv));
    aphi=h_next; atheta=std::move(theta); h_current=std::move(h_next); outer_d1=ad1;

    for(std::size_t t=2;t<=c.d;++t) {
        const std::size_t ell=t-1,j=t-1,m=c.n>>ell,half=m/2;
        BatchFoldStep step;
        auto ph=std::span<const G1>(aphi);
        auto th=std::span<const G2>(atheta);
        auto gn=gamma_level(c,t);
        auto ln=lambda_level(c,t);
        internal::append_round_metadata(tr,"REXP-BF-G1-DORY-FIRST-V1",t,j,m,half);
        step.dory_fold.d1_left=pairing_product(ph.first(half),ln);
        step.dory_fold.d1_right=pairing_product(ph.subspan(half),ln);
        step.dory_fold.d2_left=pairing_product(gn,th.first(half));
        step.dory_fold.d2_right=pairing_product(gn,th.subspan(half));
        tr.append_gt("D1Left",step.dory_fold.d1_left); tr.append_gt("D1Right",step.dory_fold.d1_right);
        tr.append_gt("D2Left",step.dory_fold.d2_left); tr.append_gt("D2Right",step.dory_fold.d2_right);
        Fr beta=tr.challenge_nonzero_fr("REXP-BF-G1-DORY-BETA-V1",ell), beta_inv=inv(beta);
        std::vector<G1> pc; std::vector<G2> tc; pc.reserve(m);tc.reserve(m);
        auto gl=gamma_level(c,ell);
        auto ll=lambda_level(c,ell);
        for(std::size_t i=0;i<m;++i) { pc.push_back(add(aphi[i],scale(gl[i],beta))); tc.push_back(add(atheta[i],scale(ll[i],beta_inv))); }
        step.dory_fold.w1=pairing_product(std::span<const G1>(pc).first(half),std::span<const G2>(tc).subspan(half));
        step.dory_fold.w2=pairing_product(std::span<const G1>(pc).subspan(half),std::span<const G2>(tc).first(half));
        internal::append_round_metadata(tr,"REXP-BF-G1-DORY-SECOND-V1",t,j,m,half);
        tr.append_gt("W1",step.dory_fold.w1);tr.append_gt("W2",step.dory_fold.w2);
        Fr alpha=tr.challenge_nonzero_fr("REXP-BF-G1-DORY-ALPHA-V1",ell),alpha_inv=inv(alpha);
        std::vector<G1> folded_phi;std::vector<G2> folded_theta;folded_phi.reserve(half);folded_theta.reserve(half);
        for(std::size_t i=0;i<half;++i) {
            folded_phi.push_back(add(scale(pc[i],alpha),pc[half+i]));
            folded_theta.push_back(add(scale(tc[i],alpha_inv),tc[half+i]));
        }
        GT fd0=products({ad0,p.x[ell],power(ad1,beta_inv),power(ad2,beta),
                         power(step.dory_fold.w1,alpha),power(step.dory_fold.w2,alpha_inv)});
        GT fd1=products({power(step.dory_fold.d1_left,alpha),step.dory_fold.d1_right,
                         power(p.x[t],mul(alpha,beta)),power(p.delta1_right[ell],beta)});
        GT fd2=products({power(step.dory_fold.d2_left,alpha_inv),step.dory_fold.d2_right,
                         power(p.x[t],mul(alpha_inv,beta_inv)),power(p.delta2_right[ell],beta_inv)});

        auto hc=std::span<const G1>(h_current);
        auto lj=lambda_level(c,j);
        internal::append_round_metadata(tr,"REXP-BF-G1-REXP-ROUND-V1",t,j,m,half);
        step.rexp_round.e=pairing_product(hc.subspan(half),lj.first(half));
        step.rexp_round.f=pairing_product(hc.first(half),lj.subspan(half));
        step.rexp_round.t_left=pairing_product(hc.first(half),ln);
        step.rexp_round.t_right=pairing_product(hc.subspan(half),ln);
        tr.append_gt("E",step.rexp_round.e);tr.append_gt("F",step.rexp_round.f);
        tr.append_gt("TLeft",step.rexp_round.t_left);tr.append_gt("TRight",step.rexp_round.t_right);
        rho=tr.challenge_nonzero_fr("REXP-BF-G1-RHO-V1",j);result.r[c.d-j-1]=rho;rho_inv=inv(rho);
        auto fresh_phi=fold_h(h_current,rho);auto ft=fresh_theta(lj,rho_inv);
        GT fresh0=products({outer_d1,power(step.rexp_round.e,rho),power(step.rexp_round.f,rho_inv)});
        GT fresh1=product(step.rexp_round.t_left,power(step.rexp_round.t_right,rho));
        GT fresh2=product(p.x[t],power(p.delta2_right[j],rho_inv));
        step.u=product(pairing_product(folded_phi,ft),pairing_product(fresh_phi,folded_theta));
        internal::append_round_metadata(tr,"REXP-BF-G1-BATCH-V1",t,j,m,half);tr.append_gt("U",step.u);
        Fr gamma=tr.challenge_nonzero_fr("REXP-BF-G1-GAMMA-V1",t),gamma2=mul(gamma,gamma);
        ad0=products({power(fd0,gamma2),power(step.u,gamma),fresh0});
        ad1=product(power(fd1,gamma),fresh1);ad2=product(power(fd2,gamma),fresh2);
        aphi.clear();atheta.clear();aphi.reserve(half);atheta.reserve(half);
        for(std::size_t i=0;i<half;++i) {
            aphi.push_back(add(scale(folded_phi[i],gamma),fresh_phi[i]));
            atheta.push_back(add(scale(folded_theta[i],gamma),ft[i]));
        }
        h_current=std::move(fresh_phi);outer_d1=fresh1;result.proof.steps.push_back(std::move(step));
    }
    if(aphi.size()!=1||atheta.size()!=1||h_current.size()!=1) throw std::logic_error("terminal dimension invariant failed");
    result.proof.phi_final=aphi[0];result.proof.theta_final=atheta[0];
    internal::append_round_metadata(tr,"REXP-BF-G1-DORY-FINAL-V1",c.d,c.d-1,1,1);
    tr.append_g1("phi_final",result.proof.phi_final);tr.append_g2("theta_final",result.proof.theta_final);
    (void)tr.challenge_nonzero_fr("REXP-BF-G1-DORY-EPSILON-V1",c.d);
    result.proof.r_final=h_current[0];
    tr.append_g1("R",result.proof.r_final);
    return result;
}

}

#include "blsagg/protocol.hpp"
#include "blsagg/serialization.hpp"
#include "blsagg/transcript.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>
#include <type_traits>
#include <unordered_set>

namespace blsagg {
namespace {

Fr fm(const Fr& a, const Fr& b) { Fr z; Fr::mul(z, a, b); return z; }
Fr inv(const Fr& a) {
  if (a.isZero()) throw std::invalid_argument("zero challenge");
  Fr z; Fr::inv(z, a); return z;
}
G1 add(const G1& a, const G1& b) { G1 z; G1::add(z, a, b); return z; }
G2 add(const G2& a, const G2& b) { G2 z; G2::add(z, a, b); return z; }
G1 mul(const G1& a, const Fr& x) { G1 z; G1::mul(z, a, x); return z; }
G2 mul(const G2& a, const Fr& x) { G2 z; G2::mul(z, a, x); return z; }
GT gmul(const GT& a, const GT& b) { GT z; GT::mul(z, a, b); return z; }
GT gpow(const GT& a, const Fr& x) { GT z; GT::pow(z, a, x); return z; }
GT product(std::initializer_list<GT> xs) {
  GT z; z.setOne();
  for (const auto& x : xs) z = gmul(z, x);
  return z;
}
GT pair(const G1& a, const G2& b) { GT z; mcl::bn::pairing(z, a, b); return z; }

struct SymbolicGTArena {
  std::vector<GT> bases;
};
struct SymbolicGT {
  std::vector<Fr> coefficients;
};
struct SymbolicTarget {
  SymbolicGT D0, D1, D2;
};
Fr fr_one() { Fr x; x = 1; return x; }
SymbolicGT symbolic_atom(SymbolicGTArena& arena, const GT& base) {
  arena.bases.push_back(base);
  SymbolicGT out; out.coefficients.resize(arena.bases.size());
  for (auto& x : out.coefficients) x.clear();
  out.coefficients.back() = fr_one();
  return out;
}
SymbolicGT symbolic_scale(const SymbolicGT& in, const Fr& scalar) {
  SymbolicGT out; out.coefficients.resize(in.coefficients.size());
  for (std::size_t i=0;i<in.coefficients.size();++i)
    Fr::mul(out.coefficients[i],in.coefficients[i],scalar);
  return out;
}
SymbolicGT symbolic_sum(std::initializer_list<SymbolicGT> inputs) {
  std::size_t n=0;for(const auto& x:inputs)n=std::max(n,x.coefficients.size());
  SymbolicGT out;out.coefficients.resize(n);for(auto& x:out.coefficients)x.clear();
  for(const auto& input:inputs)for(std::size_t i=0;i<input.coefficients.size();++i)
    Fr::add(out.coefficients[i],out.coefficients[i],input.coefficients[i]);
  return out;
}
GT symbolic_evaluate(const SymbolicGTArena& arena, const SymbolicGT& expression) {
  std::vector<GT> bases;std::vector<Fr> coefficients;
  bases.reserve(expression.coefficients.size());coefficients.reserve(expression.coefficients.size());
  for(std::size_t i=0;i<expression.coefficients.size();++i)if(!expression.coefficients[i].isZero()){
    bases.push_back(arena.bases.at(i));coefficients.push_back(expression.coefficients[i]);}
  GT out;out.setOne();
  if(!bases.empty())GT::powVec(out,bases.data(),coefficients.data(),bases.size());
  return out;
}
DoryTarget symbolic_evaluate(const SymbolicGTArena& arena, const SymbolicTarget& target) {
  return {symbolic_evaluate(arena,target.D0),symbolic_evaluate(arena,target.D1),
          symbolic_evaluate(arena,target.D2)};
}
SymbolicTarget symbolic_batch(SymbolicGTArena& arena,const SymbolicTarget& old,
                              const SymbolicTarget& fresh,const GT& u,const Fr& gamma){
  return {symbolic_sum({symbolic_scale(old.D0,fm(gamma,gamma)),
                        symbolic_scale(symbolic_atom(arena,u),gamma),fresh.D0}),
          symbolic_sum({symbolic_scale(old.D1,gamma),fresh.D1}),
          symbolic_sum({symbolic_scale(old.D2,gamma),fresh.D2})};
}
SymbolicTarget symbolic_fold(SymbolicGTArena& arena,const Precomputation&a,
                             const SymbolicTarget&in,const DoryStep&s,std::size_t level,
                             const Fr&beta,const Fr&alpha){
  const auto bi=inv(beta),ai=inv(alpha);
  return {symbolic_sum({in.D0,symbolic_atom(arena,a.X[level]),symbolic_scale(in.D1,bi),
            symbolic_scale(in.D2,beta),symbolic_scale(symbolic_atom(arena,s.W1),alpha),
            symbolic_scale(symbolic_atom(arena,s.W2),ai)}),
          symbolic_sum({symbolic_scale(symbolic_atom(arena,s.A1L),alpha),
            symbolic_atom(arena,s.A1R),
            symbolic_scale(symbolic_atom(arena,a.delta1L[level]),fm(alpha,beta)),
            symbolic_scale(symbolic_atom(arena,a.delta1R[level]),beta)}),
          symbolic_sum({symbolic_scale(symbolic_atom(arena,s.A2L),ai),
            symbolic_atom(arena,s.A2R),
            symbolic_scale(symbolic_atom(arena,a.delta2L[level]),fm(ai,bi)),
            symbolic_scale(symbolic_atom(arena,a.delta2R[level]),bi)})};
}

template<class G> bool canonical(const G& x) {
  try {
    const auto b = serialize(x);
    G y;
    return y.deserialize(b.data(), b.size()) == b.size() && y == x;
  } catch (...) { return false; }
}
bool valid(const G1& x, bool nz = false) {
  return canonical(x) && x.isValid() && x.isValidOrder() && (!nz || !x.isZero());
}
bool valid(const G2& x, bool nz = false) {
  return canonical(x) && x.isValid() && x.isValidOrder() && (!nz || !x.isZero());
}
bool valid(const GT& x) { return canonical(x) && mcl::bn::isValidGT(x); }

G1 msm(std::span<const G1> p, std::span<const Fr> s) {
  if (p.size() != s.size()) throw std::invalid_argument("G1 MSM length");
  G1 z; z.clear();
  if (!p.empty()) {
    std::vector<G1> work(p.begin(), p.end());
    G1::mulVec(z, work.data(), s.data(), work.size());
  }
  return z;
}
G2 msm(std::span<const G2> p, std::span<const Fr> s) {
  if (p.size() != s.size()) throw std::invalid_argument("G2 MSM length");
  G2 z; z.clear();
  if (!p.empty()) {
    std::vector<G2> work(p.begin(), p.end());
    G2::mulVec(z, work.data(), s.data(), work.size());
  }
  return z;
}

template<class G> std::vector<G> fold(std::span<const G> x, const Fr& r) {
  if (x.size() < 2 || x.size() % 2) throw std::invalid_argument("fold dimension");
  const std::size_t h = x.size() / 2;
  std::vector<G> out; out.reserve(h);
  for (std::size_t i = 0; i < h; ++i) out.push_back(add(x[i], mul(x[h + i], r)));
  return out;
}
template<class G> std::vector<G> weighted_constant(const G& x, std::span<const Fr> w) {
  std::vector<G> out; out.reserve(w.size());
  for (const auto& a : w) out.push_back(mul(x, a));
  return out;
}

void append(Bytes& out, std::span<const std::uint8_t> x) {
  out.insert(out.end(), x.begin(), x.end());
}
Digest pp_digest(const PublicParameters& p) {
  Bytes b;
  append_frame(b, "bls-agg-bf/pp/v1");
  append_frame(b, "BN254/mcl-v3.00");
  append_frame(b, p.mode == AggregationMode::BasicDistinct ? "basic-distinct" : "augmented");
  append_frame(b, encode_u64(p.k)); append_frame(b, encode_u64(p.d));
  append_frame(b, serialize(p.H));
  for (const auto& x : p.Gamma) append_frame(b, serialize(x));
  for (const auto& x : p.Lambda) append_frame(b, serialize(x));
  append_frame(b, serialize(p.L)); append_frame(b, serialize(p.Lprime));
  return sha256(b);
}
Digest aux_digest(const PublicParameters& p, const Precomputation& a) {
  Bytes b; append_frame(b, "bls-agg-bf/aux/v1"); append_frame(b, p.digest);
  for (const auto& level : a.gamma_chain) for (const auto& x : level) append_frame(b, serialize(x));
  for (const auto& level : a.lambda_chain) for (const auto& x : level) append_frame(b, serialize(x));
  for (const auto* v : {&a.X,&a.delta1L,&a.delta1R,&a.delta2L,&a.delta2R})
    for (const auto& x : *v) append_frame(b, serialize(x));
  for (const auto* c : {&a.g1_round0,&a.g2_round0}) {
    append_frame(b, serialize(c->E)); append_frame(b, serialize(c->F));
    append_frame(b, serialize(c->TL)); append_frame(b, serialize(c->TR));
  }
  return sha256(b);
}
template<class G> G hash_point(std::string_view tag, std::string_view seed, std::size_t i = 0) {
  Bytes b; append_frame(b, tag); append_frame(b, seed); append_frame(b, encode_u64(i));
  G out;
  if constexpr (std::is_same_v<G, G1>) mcl::bn::hashAndMapToG1(out, b.data(), b.size());
  else mcl::bn::hashAndMapToG2(out, b.data(), b.size());
  return out;
}
bool distinct(const std::vector<Bytes>& messages) {
  for (std::size_t i = 0; i < messages.size(); ++i)
    for (std::size_t j = i + 1; j < messages.size(); ++j)
      if (messages[i] == messages[j]) return false;
  return true;
}
bool statement_shape(const PublicParameters& p, const Statement& s) {
  if (s.messages.size() != p.k || s.public_keys.size() != p.k || !valid(s.sigma_agg, true)) return false;
  if (p.mode == AggregationMode::BasicDistinct && !distinct(s.messages)) return false;
  for (const auto& x : s.public_keys) if (!valid(x, true)) return false;
  return true;
}
bool pp_shape(const PublicParameters& p) {
  if (p.d < 1 || p.d >= std::numeric_limits<std::size_t>::digits ||
      p.k != (std::size_t{1} << p.d) || p.Gamma.size() != p.k || p.Lambda.size() != p.k ||
      (p.mode != AggregationMode::BasicDistinct && p.mode != AggregationMode::Augmented)) return false;
  if (!valid(p.H, true) || !valid(p.L, true) || !valid(p.Lprime, true)) return false;
  for (const auto& x : p.Gamma) if (!valid(x, true)) return false;
  for (const auto& x : p.Lambda) if (!valid(x, true)) return false;
  return pp_digest(p) == p.digest;
}
bool targets_equal(const DoryTarget& a, const DoryTarget& b) {
  return a.D0 == b.D0 && a.D1 == b.D1 && a.D2 == b.D2;
}
bool claims_equal(const RexpClaims& a, const RexpClaims& b) {
  return a.E==b.E && a.F==b.F && a.TL==b.TL && a.TR==b.TR;
}
bool aux_equal(const Precomputation& a, const Precomputation& b) {
  if (a.gamma_chain != b.gamma_chain || a.lambda_chain != b.lambda_chain || a.X != b.X ||
      a.delta1L != b.delta1L || a.delta1R != b.delta1R ||
      a.delta2L != b.delta2L || a.delta2R != b.delta2R ||
      !claims_equal(a.g1_round0,b.g1_round0) || !claims_equal(a.g2_round0,b.g2_round0))
    return false;
  return a.digest == b.digest;
}

void absorb_claim(Transcript& tr, std::string_view label, std::size_t q, const RexpClaims& c) {
  std::array<Bytes,4> f{serialize(c.E),serialize(c.F),serialize(c.TL),serialize(c.TR)};
  tr.absorb(label, q, f);
}
void absorb_targets(Transcript& tr, const DoryTarget& a, const DoryTarget& b, const DoryTarget& c) {
  std::array<Bytes,9> f{serialize(a.D0),serialize(a.D1),serialize(a.D2),
    serialize(b.D0),serialize(b.D1),serialize(b.D2),
    serialize(c.D0),serialize(c.D1),serialize(c.D2)};
  tr.absorb("bls-agg-bf/application-targets", 0, f);
}
void absorb_step_a(Transcript& tr, std::size_t level, const DoryStep& s) {
  std::array<Bytes,4> f{serialize(s.A1L),serialize(s.A1R),serialize(s.A2L),serialize(s.A2R)};
  tr.absorb("bls-agg-bf/dory-a", level, f);
}
void absorb_step_w(Transcript& tr, std::size_t level, const DoryStep& s) {
  std::array<Bytes,2> f{serialize(s.W1),serialize(s.W2)};
  tr.absorb("bls-agg-bf/dory-w", level, f);
}

DoryInstance batch(const DoryInstance& old, const DoryInstance& fresh, const GT& u, const Fr& gamma) {
  if (old.witness.Phi.size()!=fresh.witness.Phi.size() ||
      old.witness.Theta.size()!=fresh.witness.Theta.size()) throw std::logic_error("batch dimension");
  DoryInstance out;
  out.target.D0=product({gpow(old.target.D0,fm(gamma,gamma)),gpow(u,gamma),fresh.target.D0});
  out.target.D1=gmul(gpow(old.target.D1,gamma),fresh.target.D1);
  out.target.D2=gmul(gpow(old.target.D2,gamma),fresh.target.D2);
  for (std::size_t i=0;i<old.witness.Phi.size();++i) {
    out.witness.Phi.push_back(add(mul(old.witness.Phi[i],gamma),fresh.witness.Phi[i]));
    out.witness.Theta.push_back(add(mul(old.witness.Theta[i],gamma),fresh.witness.Theta[i]));
  }
  return out;
}
DoryTarget batch_target(const DoryTarget& old, const DoryTarget& fresh, const GT& u, const Fr& gamma) {
  return {product({gpow(old.D0,fm(gamma,gamma)),gpow(u,gamma),fresh.D0}),
          gmul(gpow(old.D1,gamma),fresh.D1), gmul(gpow(old.D2,gamma),fresh.D2)};
}
GT cross(const DoryWitness& a, const DoryWitness& b) {
  return gmul(direct_pairing_product(a.Phi,b.Theta),direct_pairing_product(b.Phi,a.Theta));
}

struct RexpProverData {
  std::vector<DoryInstance> g1, g2;
  std::vector<Fr> r, rp;
};

RexpProverData prove_rexp(const PublicParameters& p, const Precomputation& a,
                          Transcript& tr, Proof& proof) {
  RexpProverData out;
  out.g1.resize(p.d+1); out.g2.resize(p.d+1); out.r.resize(p.d); out.rp.resize(p.d);
  std::vector<G1> h1=p.Gamma; std::vector<G2> h2=p.Lambda;
  GT d1=a.X[0], d2=a.X[0];
  for (std::size_t q=0;q<p.d;++q) {
    const std::size_t m=p.k>>q,h=m/2;
    RexpClaims c1,c2;
    if (q==0) { c1=a.g1_round0; c2=a.g2_round0; }
    else {
      c1={direct_pairing_product(std::span(h1).subspan(h,h),std::span(a.lambda_chain[q]).first(h)),
          direct_pairing_product(std::span(h1).first(h),std::span(a.lambda_chain[q]).subspan(h,h)),
          direct_pairing_product(std::span(h1).first(h),a.lambda_chain[q+1]),
          direct_pairing_product(std::span(h1).subspan(h,h),a.lambda_chain[q+1])};
      c2={direct_pairing_product(std::span(a.gamma_chain[q]).first(h),std::span(h2).subspan(h,h)),
          direct_pairing_product(std::span(a.gamma_chain[q]).subspan(h,h),std::span(h2).first(h)),
          direct_pairing_product(a.gamma_chain[q+1],std::span(h2).first(h)),
          direct_pairing_product(a.gamma_chain[q+1],std::span(h2).subspan(h,h))};
      proof.g1_rexp_claims.push_back(c1); proof.g2_rexp_claims.push_back(c2);
    }
    absorb_claim(tr,"bls-agg-bf/rexp-g1-claim",q,c1);
    absorb_claim(tr,"bls-agg-bf/rexp-g2-claim",q,c2);
    const Fr r=tr.challenge_nonzero("bls-agg-bf/rexp-g1-challenge",q);
    const Fr rp=tr.challenge_nonzero("bls-agg-bf/rexp-g2-challenge",q);
    out.r[p.d-q-1]=r; out.rp[p.d-q-1]=rp;
    auto nh1=fold<G1>(h1,r); auto nh2=fold<G2>(h2,rp);
    DoryInstance i1,i2;
    i1.target={product({d1,gpow(c1.E,r),gpow(c1.F,inv(r))}),
      gmul(c1.TL,gpow(c1.TR,r)),
      gmul(a.delta2L[q],gpow(a.delta2R[q],inv(r)))};
    i1.witness.Phi=nh1;
    i1.witness.Theta=fold<G2>(a.lambda_chain[q],inv(r));
    i2.target={product({d2,gpow(c2.E,rp),gpow(c2.F,inv(rp))}),
      gmul(a.delta1L[q],gpow(a.delta1R[q],inv(rp))),
      gmul(c2.TL,gpow(c2.TR,rp))};
    i2.witness.Phi=fold<G1>(a.gamma_chain[q],inv(rp));
    i2.witness.Theta=nh2;
    out.g1[q+1]=std::move(i1); out.g2[q+1]=std::move(i2);
    h1=std::move(nh1); h2=std::move(nh2);
    d1=out.g1[q+1].target.D1; d2=out.g2[q+1].target.D2;
  }
  proof.R_Gamma=h1.at(0); proof.R_Lambda=h2.at(0);
  std::array<Bytes,2> f{serialize(proof.R_Gamma),serialize(proof.R_Lambda)};
  tr.absorb("bls-agg-bf/rexp-final",0,f);
  return out;
}

std::pair<DoryInstance,DoryStep> fold_dory(const PublicParameters& p,const Precomputation&a,
    const DoryInstance& in,std::size_t level,Transcript& tr) {
  const std::size_t m=in.witness.Phi.size(),h=m/2;
  if (m!=(p.k>>level)||m<2||in.witness.Theta.size()!=m) throw std::logic_error("Dory dimension");
  DoryStep s;
  s.A1L=direct_pairing_product(std::span(in.witness.Phi).first(h),a.lambda_chain[level+1]);
  s.A1R=direct_pairing_product(std::span(in.witness.Phi).subspan(h,h),a.lambda_chain[level+1]);
  s.A2L=direct_pairing_product(a.gamma_chain[level+1],std::span(in.witness.Theta).first(h));
  s.A2R=direct_pairing_product(a.gamma_chain[level+1],std::span(in.witness.Theta).subspan(h,h));
  absorb_step_a(tr,level,s);
  const Fr beta=tr.challenge_nonzero("bls-agg-bf/dory-beta",level), bi=inv(beta);
  std::vector<G1> pc; std::vector<G2> tc; pc.reserve(m);tc.reserve(m);
  for(std::size_t i=0;i<m;++i){pc.push_back(add(in.witness.Phi[i],mul(a.gamma_chain[level][i],beta)));
    tc.push_back(add(in.witness.Theta[i],mul(a.lambda_chain[level][i],bi)));}
  s.W1=direct_pairing_product(std::span(pc).first(h),std::span(tc).subspan(h,h));
  s.W2=direct_pairing_product(std::span(pc).subspan(h,h),std::span(tc).first(h));
  absorb_step_w(tr,level,s);
  const Fr alpha=tr.challenge_nonzero("bls-agg-bf/dory-alpha",level), ai=inv(alpha);
  DoryInstance out;
  for(std::size_t i=0;i<h;++i){out.witness.Phi.push_back(add(mul(pc[i],alpha),pc[h+i]));
    out.witness.Theta.push_back(add(mul(tc[i],ai),tc[h+i]));}
  out.target.D0=product({in.target.D0,a.X[level],gpow(in.target.D1,bi),gpow(in.target.D2,beta),
    gpow(s.W1,alpha),gpow(s.W2,ai)});
  out.target.D1=product({gpow(s.A1L,alpha),s.A1R,gpow(a.delta1L[level],fm(alpha,beta)),gpow(a.delta1R[level],beta)});
  out.target.D2=product({gpow(s.A2L,ai),s.A2R,gpow(a.delta2L[level],fm(ai,bi)),gpow(a.delta2R[level],bi)});
  return {std::move(out),s};
}

DoryTarget fold_target(const Precomputation&a,const DoryTarget&in,const DoryStep&s,
    std::size_t level,Transcript&tr,Fr* beta_out=nullptr,Fr* alpha_out=nullptr){
  absorb_step_a(tr,level,s);const Fr beta=tr.challenge_nonzero("bls-agg-bf/dory-beta",level),bi=inv(beta);
  absorb_step_w(tr,level,s);const Fr alpha=tr.challenge_nonzero("bls-agg-bf/dory-alpha",level),ai=inv(alpha);
  if(beta_out)*beta_out=beta;if(alpha_out)*alpha_out=alpha;
  return {product({in.D0,a.X[level],gpow(in.D1,bi),gpow(in.D2,beta),gpow(s.W1,alpha),gpow(s.W2,ai)}),
    product({gpow(s.A1L,alpha),s.A1R,gpow(a.delta1L[level],fm(alpha,beta)),gpow(a.delta1R[level],beta)}),
    product({gpow(s.A2L,ai),s.A2R,gpow(a.delta2L[level],fm(ai,bi)),gpow(a.delta2R[level],bi)})};
}

bool proof_shape(const PublicParameters&p,const Proof&v){
  if(v.g1_rexp_claims.size()!=p.d-1||v.g2_rexp_claims.size()!=p.d-1||
     v.dory_steps.size()!=p.d||v.insert_g1_u.size()!=p.d||v.insert_g2_u.size()!=p.d)return false;
  for(const auto* x:{&v.cm_M,&v.cm_pk,&v.T,&v.U1,&v.U2})if(!valid(*x))return false;
  if(!valid(v.R_Gamma)||!valid(v.R_Lambda)||!valid(v.Phi_final)||!valid(v.Theta_final))return false;
  for(const auto&c:v.g1_rexp_claims)if(!valid(c.E)||!valid(c.F)||!valid(c.TL)||!valid(c.TR))return false;
  for(const auto&c:v.g2_rexp_claims)if(!valid(c.E)||!valid(c.F)||!valid(c.TL)||!valid(c.TR))return false;
  for(const auto&s:v.dory_steps)for(const auto*x:{&s.A1L,&s.A1R,&s.A2L,&s.A2R,&s.W1,&s.W2})if(!valid(*x))return false;
  for(const auto&x:v.insert_g1_u)if(!valid(x))return false;
  for(const auto&x:v.insert_g2_u)if(!valid(x))return false;
  return true;
}

Digest context_binding(const PublicParameters&p,const Precomputation&a,const Statement&s,
                       std::span<const G1> message_points) {
  Bytes b;append_frame(b,"bls-agg-bf/validated-context/v1");
  append_frame(b,p.digest);append_frame(b,a.digest);append_frame(b,serialize(s.sigma_agg));
  for(const auto&m:s.messages)append_frame(b,m);
  for(const auto&x:s.public_keys)append_frame(b,serialize(x));
  for(const auto&x:message_points)append_frame(b,serialize(x));
  return sha256(b);
}

}

ValidatedVerifierContext::ValidatedVerifierContext(
    PublicParameters p,Precomputation a,Statement s,std::vector<G1> m,Digest b)
    :pp_(std::move(p)),aux_(std::move(a)),statement_(std::move(s)),
     message_points_(std::move(m)),binding_(b){}

GT direct_pairing_product(std::span<const G1> a, std::span<const G2> b) {
  if (a.size()!=b.size()) throw std::invalid_argument("pairing product length");
  GT z;if(a.empty()){z.setOne();return z;}GT ml;
  mcl::bn::millerLoopVec(ml,a.data(),b.data(),a.size(),true);mcl::bn::finalExp(z,ml);return z;
}
std::vector<Fr> tensor_vector(std::span<const Fr> r) {
  std::vector<Fr> w(1);w[0]=1;
  for(const auto&x:r){const auto old=w;w.reserve(old.size()*2);for(const auto&a:old)w.push_back(fm(a,x));}
  return w;
}

SetupResult setup(std::size_t d, AggregationMode mode, std::string_view seed) {
  initialize();
  if(d<1||d>=std::numeric_limits<std::size_t>::digits)throw std::invalid_argument("invalid d");
  SetupResult out;auto&p=out.pp;p.d=d;p.k=std::size_t{1}<<d;p.mode=mode;
  p.H=hash_point<G2>("bls-agg-bf/H",seed);p.L=hash_point<G1>("bls-agg-bf/L",seed);
  p.Lprime=hash_point<G2>("bls-agg-bf/Lprime",seed);
  for(std::size_t i=0;i<p.k;++i){p.Gamma.push_back(hash_point<G1>("bls-agg-bf/Gamma",seed,i));
    p.Lambda.push_back(hash_point<G2>("bls-agg-bf/Lambda",seed,i));}
  p.digest=pp_digest(p);out.aux=precompute(p);return out;
}
Precomputation precompute(const PublicParameters&p){
  if(!pp_shape(p))throw std::invalid_argument("invalid public parameters");
  Precomputation a;a.gamma_chain.push_back(p.Gamma);a.lambda_chain.push_back(p.Lambda);
  for(std::size_t j=1;j<=p.d;++j){const auto n=p.k>>j;
    a.gamma_chain.emplace_back(a.gamma_chain[j-1].begin(),a.gamma_chain[j-1].begin()+n);
    a.lambda_chain.emplace_back(a.lambda_chain[j-1].begin(),a.lambda_chain[j-1].begin()+n);}
  for(std::size_t j=0;j<=p.d;++j)a.X.push_back(direct_pairing_product(a.gamma_chain[j],a.lambda_chain[j]));
  for(std::size_t j=0;j<p.d;++j){const auto m=p.k>>j,h=m/2;
    a.delta1L.push_back(direct_pairing_product(std::span(a.gamma_chain[j]).first(h),a.lambda_chain[j+1]));
    a.delta1R.push_back(direct_pairing_product(std::span(a.gamma_chain[j]).subspan(h,h),a.lambda_chain[j+1]));
    a.delta2L.push_back(direct_pairing_product(a.gamma_chain[j+1],std::span(a.lambda_chain[j]).first(h)));
    a.delta2R.push_back(direct_pairing_product(a.gamma_chain[j+1],std::span(a.lambda_chain[j]).subspan(h,h)));}
  const auto h=p.k/2;
  a.g1_round0={direct_pairing_product(std::span(p.Gamma).subspan(h,h),std::span(p.Lambda).first(h)),
    direct_pairing_product(std::span(p.Gamma).first(h),std::span(p.Lambda).subspan(h,h)),
    a.delta1L[0],a.delta1R[0]};
  a.g2_round0={a.g1_round0.F,a.g1_round0.E,a.delta2L[0],a.delta2R[0]};
  a.digest=aux_digest(p,a);return a;
}
std::vector<G1> hash_messages(const PublicParameters&p,const Statement&s){
  if(s.messages.size()!=p.k||s.public_keys.size()!=p.k)throw std::invalid_argument("statement length");
  std::vector<G1> out;out.reserve(p.k);
  for(std::size_t i=0;i<p.k;++i){Bytes b;
    if(p.mode==AggregationMode::BasicDistinct)b=s.messages[i];
    else{append(b,serialize(s.public_keys[i]));append(b,encode_u64(s.messages[i].size()));append(b,s.messages[i]);}
    G1 x;mcl::bn::hashAndMapToG1(x,b.data(),b.size());out.push_back(x);}
  return out;
}

Proof prove(const PublicParameters&p,const Precomputation&a,const Statement&s){
  if(!pp_shape(p)||!statement_shape(p,s)||!aux_equal(a,precompute(p)))throw std::invalid_argument("invalid prove input");
  const auto M=hash_messages(p,s);Transcript tr(p,s,M);Proof proof;
  proof.cm_M=direct_pairing_product(M,p.Lambda);proof.cm_pk=direct_pairing_product(p.Gamma,s.public_keys);
  proof.T=direct_pairing_product(M,s.public_keys);
  std::array<Bytes,3> claims{serialize(proof.cm_M),serialize(proof.cm_pk),serialize(proof.T)};
  tr.absorb("bls-agg-bf/claims",0,claims);
  auto rdata=prove_rexp(p,a,tr,proof);const auto sw=tensor_vector(rdata.r),tw=tensor_vector(rdata.rp);
  const G1 Y_M=msm(M,sw);const G2 Y_pk=msm(s.public_keys,tw);
  DoryInstance im,ip,it;
  im.target={pair(Y_M,p.Lprime),proof.cm_M,pair(proof.R_Gamma,p.Lprime)};
  im.witness.Phi=M;im.witness.Theta=weighted_constant(p.Lprime,sw);
  ip.target={pair(p.L,Y_pk),pair(p.L,proof.R_Lambda),proof.cm_pk};
  ip.witness.Phi=weighted_constant(p.L,tw);ip.witness.Theta=s.public_keys;
  it.target={proof.T,proof.cm_M,proof.cm_pk};it.witness.Phi=M;it.witness.Theta=s.public_keys;
  absorb_targets(tr,im.target,ip.target,it.target);
  proof.U1=cross(im.witness,ip.witness);tr.absorb("bls-agg-bf/application-u1",0,serialize(proof.U1));
  const Fr eta=tr.challenge_nonzero("bls-agg-bf/application-eta",0);
  auto agg=batch(im,ip,proof.U1,eta);
  proof.U2=cross(agg.witness,it.witness);tr.absorb("bls-agg-bf/application-u2",0,serialize(proof.U2));
  const Fr zeta=tr.challenge_nonzero("bls-agg-bf/application-zeta",0);agg=batch(agg,it,proof.U2,zeta);
  for(std::size_t j=1;j<=p.d;++j){
    auto folded=fold_dory(p,a,agg,j-1,tr);auto bar=std::move(folded.first);proof.dory_steps.push_back(folded.second);
    const GT ug=cross(bar.witness,rdata.g1[j].witness);proof.insert_g1_u.push_back(ug);
    tr.absorb("bls-agg-bf/insert-g1-u",j,serialize(ug));
    const Fr gg=tr.challenge_nonzero("bls-agg-bf/insert-g1-gamma",j);auto hat=batch(bar,rdata.g1[j],ug,gg);
    const GT ul=cross(hat.witness,rdata.g2[j].witness);proof.insert_g2_u.push_back(ul);
    tr.absorb("bls-agg-bf/insert-g2-u",j,serialize(ul));
    const Fr gl=tr.challenge_nonzero("bls-agg-bf/insert-g2-gamma",j);agg=batch(hat,rdata.g2[j],ul,gl);
  }
  proof.Phi_final=agg.witness.Phi.at(0);proof.Theta_final=agg.witness.Theta.at(0);
  std::array<Bytes,2> final{serialize(proof.Phi_final),serialize(proof.Theta_final)};
  tr.absorb("bls-agg-bf/final-opening",0,final);
  return proof;
}

std::optional<ValidatedVerifierContext> prepare_verifier_context(const PublicParameters&p,const Precomputation&a,const Statement&s){
  try{
    if(!pp_shape(p)||!statement_shape(p,s))return std::nullopt;
    const auto expected=precompute(p);if(!aux_equal(a,expected)||a.digest!=aux_digest(p,a))return std::nullopt;
    auto points=hash_messages(p,s);auto binding=context_binding(p,a,s,points);
    return ValidatedVerifierContext(p,a,s,std::move(points),binding);
  }catch(...){return std::nullopt;}
}

static VerificationTrace verify_online_core(const ValidatedVerifierContext&c,const Proof&v,
                                             bool proof_prevalidated,bool parallel_msm){
  VerificationTrace out;
  using Clock=std::chrono::steady_clock;
  const auto elapsed=[](auto x,auto y){return std::chrono::duration<double,std::milli>(y-x).count();};
  const auto total_start=Clock::now();
  try{
    const auto&p=c.parameters();const auto&a=c.precomputation();const auto&s=c.statement();
    const auto message_points=c.message_points();
    if(context_binding(p,a,s,message_points)!=c.binding())return out;
    auto tick=Clock::now();
    if(!proof_prevalidated&&!proof_shape(p,v))return out;
    auto tock=Clock::now();out.proof_validation_ms=elapsed(tick,tock);tick=Clock::now();
    Transcript tr(p,s,message_points);std::array<Bytes,3> claims{serialize(v.cm_M),serialize(v.cm_pk),serialize(v.T)};
    tr.absorb("bls-agg-bf/claims",0,claims);
    std::vector<DoryTarget>g1(p.d+1),g2(p.d+1);std::vector<Fr>rs(p.d),rt(p.d);GT d1=a.X[0],d2=a.X[0];
    out.g1_rexp_challenges.resize(p.d);out.g2_rexp_challenges.resize(p.d);
    double rexp_ms=0;
    for(std::size_t q=0;q<p.d;++q){const auto c1=q? v.g1_rexp_claims[q-1]:a.g1_round0;
      const auto c2=q? v.g2_rexp_claims[q-1]:a.g2_round0;
      absorb_claim(tr,"bls-agg-bf/rexp-g1-claim",q,c1);absorb_claim(tr,"bls-agg-bf/rexp-g2-claim",q,c2);
      const Fr r=tr.challenge_nonzero("bls-agg-bf/rexp-g1-challenge",q),rp=tr.challenge_nonzero("bls-agg-bf/rexp-g2-challenge",q);
      auto phase=Clock::now();
      rs[p.d-q-1]=r;rt[p.d-q-1]=rp;
      out.g1_rexp_challenges[q]=r;out.g2_rexp_challenges[q]=rp;
      g1[q+1]={product({d1,gpow(c1.E,r),gpow(c1.F,inv(r))}),gmul(c1.TL,gpow(c1.TR,r)),
        gmul(a.delta2L[q],gpow(a.delta2R[q],inv(r)))};
      g2[q+1]={product({d2,gpow(c2.E,rp),gpow(c2.F,inv(rp))}),gmul(a.delta1L[q],gpow(a.delta1R[q],inv(rp))),
        gmul(c2.TL,gpow(c2.TR,rp))};d1=g1[q+1].D1;d2=g2[q+1].D2;
      auto phase_end=Clock::now();rexp_ms+=elapsed(phase,phase_end);}
    std::array<Bytes,2> rf{serialize(v.R_Gamma),serialize(v.R_Lambda)};
    tr.absorb("bls-agg-bf/rexp-final",0,rf);
    out.rexp_gt_recurrence_ms=rexp_ms;tick=Clock::now();
    out.g1_rexp_targets=g1;out.g2_rexp_targets=g2;
    const auto sw=tensor_vector(rs),tw=tensor_vector(rt);
    tock=Clock::now();out.tensor_reconstruction_ms=elapsed(tick,tock);
    const auto msm_wall_start=Clock::now();
    if(parallel_msm){
      std::exception_ptr worker_error;
      std::thread g1_worker([&]{
        const auto start=Clock::now();
        try{out.Y_M=msm(message_points,sw);}
        catch(...){worker_error=std::current_exception();}
        out.message_g1_msm_ms=elapsed(start,Clock::now());
      });
      const auto g2_start=Clock::now();
      try{out.Y_pk=msm(s.public_keys,tw);}
      catch(...){g1_worker.join();throw;}
      out.public_key_g2_msm_ms=elapsed(g2_start,Clock::now());
      g1_worker.join();
      if(worker_error)std::rethrow_exception(worker_error);
    }else{
      auto start=Clock::now();out.Y_M=msm(message_points,sw);
      out.message_g1_msm_ms=elapsed(start,Clock::now());start=Clock::now();out.Y_pk=msm(s.public_keys,tw);
      out.public_key_g2_msm_ms=elapsed(start,Clock::now());
    }
    out.public_input_msm_wall_ms=elapsed(msm_wall_start,Clock::now());tick=Clock::now();
    const double transcript_before_application=tr.timing_ms();
    DoryTarget im{pair(out.Y_M,p.Lprime),v.cm_M,pair(v.R_Gamma,p.Lprime)};
    DoryTarget ip{pair(p.L,out.Y_pk),pair(p.L,v.R_Lambda),v.cm_pk};
    DoryTarget it{v.T,v.cm_M,v.cm_pk};absorb_targets(tr,im,ip,it);
    tr.absorb("bls-agg-bf/application-u1",0,serialize(v.U1));const Fr eta=tr.challenge_nonzero("bls-agg-bf/application-eta",0);
    auto agg=batch_target(im,ip,v.U1,eta);tr.absorb("bls-agg-bf/application-u2",0,serialize(v.U2));
    const Fr zeta=tr.challenge_nonzero("bls-agg-bf/application-zeta",0);agg=batch_target(agg,it,v.U2,zeta);
    out.eta=eta;out.zeta=zeta;out.accumulator_targets.push_back(agg);
    tock=Clock::now();out.application_batching_ms=elapsed(tick,tock)-(tr.timing_ms()-transcript_before_application);tick=Clock::now();
    const double transcript_before_dory=tr.timing_ms();
    out.dory_beta.resize(p.d);out.dory_alpha.resize(p.d);out.insert_g1_gamma.resize(p.d);out.insert_g2_gamma.resize(p.d);
    for(std::size_t j=1;j<=p.d;++j){auto bar=fold_target(a,agg,v.dory_steps[j-1],j-1,tr,&out.dory_beta[j-1],&out.dory_alpha[j-1]);
      out.dory_fold_targets.push_back(bar);
      tr.absorb("bls-agg-bf/insert-g1-u",j,serialize(v.insert_g1_u[j-1]));
      const Fr gg=tr.challenge_nonzero("bls-agg-bf/insert-g1-gamma",j);auto hat=batch_target(bar,g1[j],v.insert_g1_u[j-1],gg);
      tr.absorb("bls-agg-bf/insert-g2-u",j,serialize(v.insert_g2_u[j-1]));
      const Fr gl=tr.challenge_nonzero("bls-agg-bf/insert-g2-gamma",j);agg=batch_target(hat,g2[j],v.insert_g2_u[j-1],gl);
      out.insert_g1_gamma[j-1]=gg;out.insert_g2_gamma[j-1]=gl;out.accumulator_targets.push_back(agg);}
    tock=Clock::now();out.integrated_dory_ms=elapsed(tick,tock)-(tr.timing_ms()-transcript_before_dory);tick=Clock::now();
    std::array<Bytes,2> final{serialize(v.Phi_final),serialize(v.Theta_final)};tr.absorb("bls-agg-bf/final-opening",0,final);
    out.final_dory=agg;out.final_g1_rexp=g1[p.d];out.final_g2_rexp=g2[p.d];
    tock=Clock::now();out.final_challenges_ms=0;out.transcript_hashing_ms=tr.timing_ms();tick=Clock::now();
    out.accepted=pair(v.Phi_final,v.Theta_final)==agg.D0&&pair(v.Phi_final,a.lambda_chain[p.d][0])==agg.D1&&
      pair(a.gamma_chain[p.d][0],v.Theta_final)==agg.D2&&pair(v.R_Gamma,a.lambda_chain[p.d][0])==g1[p.d].D1&&
      pair(a.gamma_chain[p.d][0],v.R_Lambda)==g2[p.d].D2&&pair(s.sigma_agg,p.H)==v.T;
    tock=Clock::now();out.terminal_checks_ms=elapsed(tick,tock);
    out.total_online_ms=elapsed(total_start,tock);
    return out;
  }catch(...){return out;}
}

static VerificationTrace verify_online_symbolic_core(const ValidatedVerifierContext&c,
                                                       const ValidatedProof&validated,
                                                       bool capture_intermediates,
                                                       bool split_g2_msm=false){
  VerificationTrace out;
  using Clock=std::chrono::steady_clock;
  const auto elapsed=[](auto x,auto y){return std::chrono::duration<double,std::milli>(y-x).count();};
  const auto total_start=Clock::now();
  try{
    const auto&p=c.parameters();const auto&a=c.precomputation();const auto&s=c.statement();
    if(validated.parameter_digest()!=p.digest)return out;
    const auto&v=validated.proof();const auto message_points=c.message_points();
    if(context_binding(p,a,s,message_points)!=c.binding())return out;
    auto tick=Clock::now();out.proof_validation_ms=0;
    Transcript tr(p,s,message_points);std::array<Bytes,3> claims{serialize(v.cm_M),serialize(v.cm_pk),serialize(v.T)};
    tr.absorb("bls-agg-bf/claims",0,claims);
    SymbolicGTArena arena;
    std::vector<SymbolicTarget>g1(p.d+1),g2(p.d+1);std::vector<Fr>rs(p.d),rt(p.d);
    auto d1=symbolic_atom(arena,a.X[0]),d2=symbolic_atom(arena,a.X[0]);
    out.g1_rexp_challenges.resize(p.d);out.g2_rexp_challenges.resize(p.d);
    if(capture_intermediates){out.g1_rexp_targets.resize(p.d+1);out.g2_rexp_targets.resize(p.d+1);}
    double rexp_ms=0;
    for(std::size_t q=0;q<p.d;++q){const auto c1=q?v.g1_rexp_claims[q-1]:a.g1_round0;
      const auto c2=q?v.g2_rexp_claims[q-1]:a.g2_round0;
      absorb_claim(tr,"bls-agg-bf/rexp-g1-claim",q,c1);absorb_claim(tr,"bls-agg-bf/rexp-g2-claim",q,c2);
      const Fr r=tr.challenge_nonzero("bls-agg-bf/rexp-g1-challenge",q),rp=tr.challenge_nonzero("bls-agg-bf/rexp-g2-challenge",q);
      auto phase=Clock::now();const auto ri=inv(r),rpi=inv(rp);
      rs[p.d-q-1]=r;rt[p.d-q-1]=rp;out.g1_rexp_challenges[q]=r;out.g2_rexp_challenges[q]=rp;
      g1[q+1]={symbolic_sum({d1,symbolic_scale(symbolic_atom(arena,c1.E),r),
                        symbolic_scale(symbolic_atom(arena,c1.F),ri)}),
                symbolic_sum({symbolic_atom(arena,c1.TL),symbolic_scale(symbolic_atom(arena,c1.TR),r)}),
                symbolic_sum({symbolic_atom(arena,a.delta2L[q]),
                              symbolic_scale(symbolic_atom(arena,a.delta2R[q]),ri)})};
      g2[q+1]={symbolic_sum({d2,symbolic_scale(symbolic_atom(arena,c2.E),rp),
                        symbolic_scale(symbolic_atom(arena,c2.F),rpi)}),
                symbolic_sum({symbolic_atom(arena,a.delta1L[q]),
                              symbolic_scale(symbolic_atom(arena,a.delta1R[q]),rpi)}),
                symbolic_sum({symbolic_atom(arena,c2.TL),symbolic_scale(symbolic_atom(arena,c2.TR),rp)})};
      d1=g1[q+1].D1;d2=g2[q+1].D2;
      if(capture_intermediates){out.g1_rexp_targets[q+1]=symbolic_evaluate(arena,g1[q+1]);
        out.g2_rexp_targets[q+1]=symbolic_evaluate(arena,g2[q+1]);}
      rexp_ms+=elapsed(phase,Clock::now());
    }
    std::array<Bytes,2> rf{serialize(v.R_Gamma),serialize(v.R_Lambda)};
    tr.absorb("bls-agg-bf/rexp-final",0,rf);out.rexp_gt_recurrence_ms=rexp_ms;tick=Clock::now();
    const auto sw=tensor_vector(rs),tw=tensor_vector(rt);auto tock=Clock::now();
    out.tensor_reconstruction_ms=elapsed(tick,tock);
    const auto msm_wall_start=Clock::now();
    if(split_g2_msm&&s.public_keys.size()>=2){
      const auto half=s.public_keys.size()/2;G2 g2_left,g2_right;
      std::array<std::exception_ptr,3> errors{};std::vector<std::thread> workers;workers.reserve(2);
      try{
        workers.emplace_back([&]{const auto start=Clock::now();try{out.Y_M=msm(message_points,sw);}
          catch(...){errors[0]=std::current_exception();}out.message_g1_msm_ms=elapsed(start,Clock::now());});
        workers.emplace_back([&]{const auto start=Clock::now();
          try{g2_left=msm(std::span(s.public_keys).first(half),std::span(tw).first(half));}
          catch(...){errors[1]=std::current_exception();}
          out.public_key_g2_left_msm_ms=elapsed(start,Clock::now());});
      }catch(...){for(auto& worker:workers)worker.join();throw;}
      const auto right_start=Clock::now();
      try{g2_right=msm(std::span(s.public_keys).subspan(half),std::span(tw).subspan(half));}
      catch(...){errors[2]=std::current_exception();}
      out.public_key_g2_right_msm_ms=elapsed(right_start,Clock::now());
      for(auto& worker:workers)worker.join();for(const auto& error:errors)if(error)std::rethrow_exception(error);
      out.Y_pk=add(g2_left,g2_right);
      out.public_key_g2_msm_ms=std::max(out.public_key_g2_left_msm_ms,out.public_key_g2_right_msm_ms);
    }else{
      std::exception_ptr worker_error;
      std::thread g1_worker([&]{const auto start=Clock::now();try{out.Y_M=msm(message_points,sw);}
        catch(...){worker_error=std::current_exception();}out.message_g1_msm_ms=elapsed(start,Clock::now());});
      const auto g2_start=Clock::now();try{out.Y_pk=msm(s.public_keys,tw);}
      catch(...){g1_worker.join();throw;}out.public_key_g2_msm_ms=elapsed(g2_start,Clock::now());
      g1_worker.join();if(worker_error)std::rethrow_exception(worker_error);
    }
    out.public_input_msm_wall_ms=elapsed(msm_wall_start,Clock::now());tick=Clock::now();
    const double transcript_before_application=tr.timing_ms();
    DoryTarget im{pair(out.Y_M,p.Lprime),v.cm_M,pair(v.R_Gamma,p.Lprime)};
    DoryTarget ip{pair(p.L,out.Y_pk),pair(p.L,v.R_Lambda),v.cm_pk};
    DoryTarget it{v.T,v.cm_M,v.cm_pk};absorb_targets(tr,im,ip,it);
    SymbolicTarget sim{symbolic_atom(arena,im.D0),symbolic_atom(arena,im.D1),symbolic_atom(arena,im.D2)};
    SymbolicTarget sip{symbolic_atom(arena,ip.D0),symbolic_atom(arena,ip.D1),symbolic_atom(arena,ip.D2)};
    SymbolicTarget sit{symbolic_atom(arena,it.D0),symbolic_atom(arena,it.D1),symbolic_atom(arena,it.D2)};
    tr.absorb("bls-agg-bf/application-u1",0,serialize(v.U1));const Fr eta=tr.challenge_nonzero("bls-agg-bf/application-eta",0);
    auto agg=symbolic_batch(arena,sim,sip,v.U1,eta);tr.absorb("bls-agg-bf/application-u2",0,serialize(v.U2));
    const Fr zeta=tr.challenge_nonzero("bls-agg-bf/application-zeta",0);agg=symbolic_batch(arena,agg,sit,v.U2,zeta);
    out.eta=eta;out.zeta=zeta;if(capture_intermediates)out.accumulator_targets.push_back(symbolic_evaluate(arena,agg));
    tock=Clock::now();out.application_batching_ms=elapsed(tick,tock)-(tr.timing_ms()-transcript_before_application);tick=Clock::now();
    const double transcript_before_dory=tr.timing_ms();
    out.dory_beta.resize(p.d);out.dory_alpha.resize(p.d);out.insert_g1_gamma.resize(p.d);out.insert_g2_gamma.resize(p.d);
    for(std::size_t j=1;j<=p.d;++j){const auto&step=v.dory_steps[j-1];
      absorb_step_a(tr,j-1,step);const Fr beta=tr.challenge_nonzero("bls-agg-bf/dory-beta",j-1);
      absorb_step_w(tr,j-1,step);const Fr alpha=tr.challenge_nonzero("bls-agg-bf/dory-alpha",j-1);
      auto bar=symbolic_fold(arena,a,agg,step,j-1,beta,alpha);
      if(capture_intermediates)out.dory_fold_targets.push_back(symbolic_evaluate(arena,bar));
      tr.absorb("bls-agg-bf/insert-g1-u",j,serialize(v.insert_g1_u[j-1]));
      const Fr gg=tr.challenge_nonzero("bls-agg-bf/insert-g1-gamma",j);
      auto hat=symbolic_batch(arena,bar,g1[j],v.insert_g1_u[j-1],gg);
      tr.absorb("bls-agg-bf/insert-g2-u",j,serialize(v.insert_g2_u[j-1]));
      const Fr gl=tr.challenge_nonzero("bls-agg-bf/insert-g2-gamma",j);
      agg=symbolic_batch(arena,hat,g2[j],v.insert_g2_u[j-1],gl);
      out.dory_beta[j-1]=beta;out.dory_alpha[j-1]=alpha;out.insert_g1_gamma[j-1]=gg;out.insert_g2_gamma[j-1]=gl;
      if(capture_intermediates)out.accumulator_targets.push_back(symbolic_evaluate(arena,agg));
    }
    tock=Clock::now();out.integrated_dory_ms=elapsed(tick,tock)-(tr.timing_ms()-transcript_before_dory);tick=Clock::now();
    std::array<Bytes,2> final{serialize(v.Phi_final),serialize(v.Theta_final)};tr.absorb("bls-agg-bf/final-opening",0,final);
    out.final_challenges_ms=0;out.transcript_hashing_ms=tr.timing_ms();
    const auto multiexp_start=Clock::now();
    if(capture_intermediates){
      out.final_dory=out.accumulator_targets.back();out.final_g1_rexp=out.g1_rexp_targets.back();
      out.final_g2_rexp=out.g2_rexp_targets.back();
    }else{
      std::array<GT,5> values;std::array<std::exception_ptr,5> errors{};
      const std::array<const SymbolicGT*,5> expressions{
        &agg.D0,&agg.D1,&agg.D2,&g1[p.d].D1,&g2[p.d].D2};
      std::vector<std::thread> workers;workers.reserve(4);
      try{for(std::size_t i=0;i<4;++i)workers.emplace_back([&,i]{
          try{values[i]=symbolic_evaluate(arena,*expressions[i]);}
          catch(...){errors[i]=std::current_exception();}
        });}
      catch(...){for(auto& worker:workers)worker.join();throw;}
      try{values[4]=symbolic_evaluate(arena,*expressions[4]);}
      catch(...){errors[4]=std::current_exception();}
      for(auto& worker:workers)worker.join();
      for(const auto& error:errors)if(error)std::rethrow_exception(error);
      out.final_dory={values[0],values[1],values[2]};
      out.final_g1_rexp.D1=values[3];out.final_g2_rexp.D2=values[4];
    }
    out.gt_multiexponentiation_ms=elapsed(multiexp_start,Clock::now());tick=Clock::now();
    out.accepted=pair(v.Phi_final,v.Theta_final)==out.final_dory.D0&&
      pair(v.Phi_final,a.lambda_chain[p.d][0])==out.final_dory.D1&&
      pair(a.gamma_chain[p.d][0],v.Theta_final)==out.final_dory.D2&&
      pair(v.R_Gamma,a.lambda_chain[p.d][0])==out.final_g1_rexp.D1&&
      pair(a.gamma_chain[p.d][0],v.R_Lambda)==out.final_g2_rexp.D2&&pair(s.sigma_agg,p.H)==v.T;
    tock=Clock::now();out.terminal_checks_ms=elapsed(tick,tock);out.total_online_ms=elapsed(total_start,tock);
    return out;
  }catch(...){return out;}
}

VerificationTrace verify_online_diagnostic(const ValidatedVerifierContext&c,const Proof&p){
  return verify_online_core(c,p,false,false);
}
bool verify_online(const ValidatedVerifierContext&c,const Proof&p){return verify_online_diagnostic(c,p).accepted;}
VerificationTrace verify_online_diagnostic(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,false,true);
}
bool verify_online(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_sequential_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  if(p.parameter_digest()!=c.parameters().digest)return {};
  return verify_online_core(c,p.proof(),true,false);
}
bool verify_online_sequential_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_sequential_msm_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_parallel_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  if(p.parameter_digest()!=c.parameters().digest)return {};
  return verify_online_core(c,p.proof(),true,true);
}
bool verify_online_parallel_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_parallel_msm_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_symbolic_gt_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,false,false);
}
bool verify_online_symbolic_gt(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_gt_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_symbolic_gt_differential_trace(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,true);
}
VerificationTrace verify_online_split_g2_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_diagnostic(c,p);
}
bool verify_online_split_g2_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online(c,p);
}
bool verify_safe(const PublicParameters&p,const Precomputation&a,const Statement&s,const Proof&v){
  auto c=prepare_verifier_context(p,a,s);return c&&verify_online(*c,v);
}
bool direct_bls_verify(const PublicParameters&p,const Statement&s){
  try{if(!statement_shape(p,s))return false;return pair(s.sigma_agg,p.H)==direct_pairing_product(hash_messages(p,s),s.public_keys);}
  catch(...){return false;}
}
std::size_t proof_payload_bytes(const Proof&p){
  std::size_t n=serialize(p.cm_M).size()+serialize(p.cm_pk).size()+serialize(p.T).size();
  for(const auto&c:p.g1_rexp_claims)n+=serialize(c.E).size()+serialize(c.F).size()+serialize(c.TL).size()+serialize(c.TR).size();
  for(const auto&c:p.g2_rexp_claims)n+=serialize(c.E).size()+serialize(c.F).size()+serialize(c.TL).size()+serialize(c.TR).size();
  n+=serialize(p.R_Gamma).size()+serialize(p.R_Lambda).size()+serialize(p.U1).size()+serialize(p.U2).size();
  for(const auto&s:p.dory_steps)n+=serialize(s.A1L).size()+serialize(s.A1R).size()+serialize(s.A2L).size()+serialize(s.A2R).size()+serialize(s.W1).size()+serialize(s.W2).size();
  for(const auto&x:p.insert_g1_u)n+=serialize(x).size();for(const auto&x:p.insert_g2_u)n+=serialize(x).size();
  return n+serialize(p.Phi_final).size()+serialize(p.Theta_final).size();
}
}

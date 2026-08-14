#include "blsagg/protocol.hpp"
#include "blsagg/serialization.hpp"
#include "blsagg/transcript.hpp"
#include "internal/crypto.hpp"
#include "internal/validation.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <thread>

namespace blsagg {
namespace {

using internal::add;
using internal::fm;
using internal::gmul;
using internal::gpow;
using internal::inv;
using internal::msm;
using internal::mul;
using internal::pair;
using internal::product;
using internal::valid;

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
  internal::append_raw(out, x);
}
bool targets_equal(const DoryTarget& a, const DoryTarget& b) {
  return a.D0 == b.D0 && a.D1 == b.D1 && a.D2 == b.D2;
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

class VerifierTranscriptReplay {
 public:
  VerifierTranscriptReplay(const PublicParameters& parameters,
                           const Statement& statement,
                           std::span<const G1> message_points,
                           const Proof& proof)
      : transcript_(parameters, statement, message_points) {
    const std::array<Bytes, 3> claims{serialize(proof.cm_M),
                                      serialize(proof.cm_pk),
                                      serialize(proof.T)};
    transcript_.absorb("bls-agg-bf/claims", 0, claims);
  }

  std::pair<Fr, Fr> rexp_round(std::size_t round,
                               const RexpClaims& g1_claim,
                               const RexpClaims& g2_claim) {
    absorb_claim(transcript_, "bls-agg-bf/rexp-g1-claim", round,
                 g1_claim);
    absorb_claim(transcript_, "bls-agg-bf/rexp-g2-claim", round,
                 g2_claim);
    const auto g1_challenge = transcript_.challenge_nonzero(
        "bls-agg-bf/rexp-g1-challenge", round);
    const auto g2_challenge = transcript_.challenge_nonzero(
        "bls-agg-bf/rexp-g2-challenge", round);
    return {g1_challenge, g2_challenge};
  }

  void finish_rexp(const Proof& proof) {
    const std::array<Bytes, 2> final{serialize(proof.R_Gamma),
                                     serialize(proof.R_Lambda)};
    transcript_.absorb("bls-agg-bf/rexp-final", 0, final);
  }

  std::pair<Fr, Fr> application(const DoryTarget& messages,
                                const DoryTarget& public_keys,
                                const DoryTarget& pairing_claim,
                                const Proof& proof) {
    absorb_targets(transcript_, messages, public_keys, pairing_claim);
    transcript_.absorb("bls-agg-bf/application-u1", 0,
                       serialize(proof.U1));
    const auto eta = transcript_.challenge_nonzero(
        "bls-agg-bf/application-eta", 0);
    transcript_.absorb("bls-agg-bf/application-u2", 0,
                       serialize(proof.U2));
    const auto zeta = transcript_.challenge_nonzero(
        "bls-agg-bf/application-zeta", 0);
    return {eta, zeta};
  }

  std::pair<Fr, Fr> dory_round(std::size_t level,
                               const DoryStep& step) {
    absorb_step_a(transcript_, level, step);
    const auto beta =
        transcript_.challenge_nonzero("bls-agg-bf/dory-beta", level);
    absorb_step_w(transcript_, level, step);
    const auto alpha =
        transcript_.challenge_nonzero("bls-agg-bf/dory-alpha", level);
    return {beta, alpha};
  }

  Fr insert_g1(std::size_t round, const GT& cross_term) {
    transcript_.absorb("bls-agg-bf/insert-g1-u", round,
                       serialize(cross_term));
    return transcript_.challenge_nonzero("bls-agg-bf/insert-g1-gamma",
                                         round);
  }

  Fr insert_g2(std::size_t round, const GT& cross_term) {
    transcript_.absorb("bls-agg-bf/insert-g2-u", round,
                       serialize(cross_term));
    return transcript_.challenge_nonzero("bls-agg-bf/insert-g2-gamma",
                                         round);
  }

  void bind_final(const Proof& proof) {
    const std::array<Bytes, 2> final{serialize(proof.Phi_final),
                                     serialize(proof.Theta_final)};
    transcript_.absorb("bls-agg-bf/final-opening", 0, final);
  }

  double timing_ms() const { return transcript_.timing_ms(); }

 private:
  Transcript transcript_;
};

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
    std::size_t level,const Fr&beta,const Fr&alpha){
  const auto bi=inv(beta),ai=inv(alpha);
  return {product({in.D0,a.X[level],gpow(in.D1,bi),gpow(in.D2,beta),gpow(s.W1,alpha),gpow(s.W2,ai)}),
    product({gpow(s.A1L,alpha),s.A1R,gpow(a.delta1L[level],fm(alpha,beta)),gpow(a.delta1R[level],beta)}),
    product({gpow(s.A2L,ai),s.A2R,gpow(a.delta2L[level],fm(ai,bi)),gpow(a.delta2R[level],bi)})};
}

struct ApplicationInstances {
  DoryInstance messages;
  DoryInstance public_keys;
  DoryInstance pairing_claim;
};

void commit_application_claims(const PublicParameters& parameters,
                               const Statement& statement,
                               std::span<const G1> message_points,
                               Transcript& transcript, Proof& proof) {
  proof.cm_M = direct_pairing_product(message_points, parameters.Lambda);
  proof.cm_pk =
      direct_pairing_product(parameters.Gamma, statement.public_keys);
  proof.T = direct_pairing_product(message_points, statement.public_keys);
  const std::array<Bytes, 3> claims{serialize(proof.cm_M),
                                    serialize(proof.cm_pk),
                                    serialize(proof.T)};
  transcript.absorb("bls-agg-bf/claims", 0, claims);
}

ApplicationInstances build_application_instances(
    const PublicParameters& parameters, const Statement& statement,
    std::span<const G1> message_points, const Proof& proof,
    const RexpProverData& rexp) {
  const auto message_weights = tensor_vector(rexp.r);
  const auto key_weights = tensor_vector(rexp.rp);
  const auto reduced_messages = msm(message_points, message_weights);
  const auto reduced_keys = msm(statement.public_keys, key_weights);

  ApplicationInstances result;
  result.messages.target =
      {pair(reduced_messages, parameters.Lprime), proof.cm_M,
       pair(proof.R_Gamma, parameters.Lprime)};
  result.messages.witness.Phi.assign(message_points.begin(),
                                     message_points.end());
  result.messages.witness.Theta =
      weighted_constant(parameters.Lprime, message_weights);

  result.public_keys.target =
      {pair(parameters.L, reduced_keys), pair(parameters.L, proof.R_Lambda),
       proof.cm_pk};
  result.public_keys.witness.Phi =
      weighted_constant(parameters.L, key_weights);
  result.public_keys.witness.Theta = statement.public_keys;

  result.pairing_claim.target = {proof.T, proof.cm_M, proof.cm_pk};
  result.pairing_claim.witness.Phi.assign(message_points.begin(),
                                          message_points.end());
  result.pairing_claim.witness.Theta = statement.public_keys;
  return result;
}

DoryInstance batch_application_instances(const ApplicationInstances& inputs,
                                         Transcript& transcript,
                                         Proof& proof) {
  absorb_targets(transcript, inputs.messages.target, inputs.public_keys.target,
                 inputs.pairing_claim.target);
  proof.U1 = cross(inputs.messages.witness, inputs.public_keys.witness);
  transcript.absorb("bls-agg-bf/application-u1", 0, serialize(proof.U1));
  const auto eta =
      transcript.challenge_nonzero("bls-agg-bf/application-eta", 0);
  auto accumulator =
      batch(inputs.messages, inputs.public_keys, proof.U1, eta);
  proof.U2 = cross(accumulator.witness, inputs.pairing_claim.witness);
  transcript.absorb("bls-agg-bf/application-u2", 0, serialize(proof.U2));
  const auto zeta =
      transcript.challenge_nonzero("bls-agg-bf/application-zeta", 0);
  return batch(accumulator, inputs.pairing_claim, proof.U2, zeta);
}

DoryInstance fold_and_insert_rexp(const PublicParameters& parameters,
                                  const Precomputation& precomputation,
                                  const RexpProverData& rexp,
                                  Transcript& transcript, Proof& proof,
                                  DoryInstance accumulator) {
  proof.dory_steps.reserve(parameters.d);
  proof.insert_g1_u.reserve(parameters.d);
  proof.insert_g2_u.reserve(parameters.d);
  for (std::size_t round = 1; round <= parameters.d; ++round) {
    auto [folded, step] = fold_dory(parameters, precomputation, accumulator,
                                    round - 1, transcript);
    proof.dory_steps.push_back(std::move(step));

    const auto insert_g1 = cross(folded.witness, rexp.g1[round].witness);
    proof.insert_g1_u.push_back(insert_g1);
    transcript.absorb("bls-agg-bf/insert-g1-u", round,
                      serialize(insert_g1));
    const auto gamma_g1 =
        transcript.challenge_nonzero("bls-agg-bf/insert-g1-gamma", round);
    auto with_g1 =
        batch(folded, rexp.g1[round], insert_g1, gamma_g1);

    const auto insert_g2 = cross(with_g1.witness, rexp.g2[round].witness);
    proof.insert_g2_u.push_back(insert_g2);
    transcript.absorb("bls-agg-bf/insert-g2-u", round,
                      serialize(insert_g2));
    const auto gamma_g2 =
        transcript.challenge_nonzero("bls-agg-bf/insert-g2-gamma", round);
    accumulator =
        batch(with_g1, rexp.g2[round], insert_g2, gamma_g2);
  }
  return accumulator;
}

void bind_final_opening(const DoryInstance& accumulator,
                        Transcript& transcript, Proof& proof) {
  if (accumulator.witness.Phi.size() != 1 ||
      accumulator.witness.Theta.size() != 1)
    throw std::logic_error("invalid terminal witness dimension");
  proof.Phi_final = accumulator.witness.Phi.front();
  proof.Theta_final = accumulator.witness.Theta.front();
  const std::array<Bytes, 2> final{serialize(proof.Phi_final),
                                   serialize(proof.Theta_final)};
  transcript.absorb("bls-agg-bf/final-opening", 0, final);
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
  if(!internal::valid_public_parameters(p)||!internal::valid_statement(p,s)||
     !internal::valid_precomputation(p,a))throw std::invalid_argument("invalid prove input");
  const auto message_points=hash_messages(p,s);Transcript tr(p,s,message_points);Proof proof;
  commit_application_claims(p,s,message_points,tr,proof);
  const auto rexp=prove_rexp(p,a,tr,proof);
  const auto applications=build_application_instances(p,s,message_points,proof,rexp);
  auto accumulator=batch_application_instances(applications,tr,proof);
  accumulator=fold_and_insert_rexp(p,a,rexp,tr,proof,std::move(accumulator));
  bind_final_opening(accumulator,tr,proof);
  return proof;
}

std::optional<ValidatedVerifierContext> prepare_verifier_context(const PublicParameters&p,const Precomputation&a,const Statement&s){
  try{
  if(!internal::valid_public_parameters(p)||!internal::valid_statement(p,s)||
     !internal::valid_precomputation(p,a))return std::nullopt;
  auto points=hash_messages(p,s);auto binding=internal::context_binding(p,a,s,points);
    return ValidatedVerifierContext(p,a,s,std::move(points),binding);
  }catch(...){return std::nullopt;}
}

enum class MsmStrategy { Sequential, Parallel, SplitG2 };

bool validated_proof_matches(const PublicParameters& parameters,
                             const ValidatedProof& validated) {
  if (validated.parameter_digest() != parameters.digest) return false;
  try {
    return validated.wire_binding() ==
           internal::validated_proof_binding(parameters, validated.proof());
  } catch (...) {
    return false;
  }
}

static VerificationTrace verify_online_core(const ValidatedVerifierContext&c,const Proof&v,
                                             bool proof_prevalidated,MsmStrategy strategy){
  VerificationTrace out;
  using Clock=std::chrono::steady_clock;
  const auto elapsed=[](auto x,auto y){return std::chrono::duration<double,std::milli>(y-x).count();};
  const auto total_start=Clock::now();
  try{
    const auto&p=c.parameters();const auto&a=c.precomputation();const auto&s=c.statement();
    const auto message_points=c.message_points();
  if(internal::context_binding(p,a,s,message_points)!=c.binding())return out;
    auto tick=Clock::now();
  if(!proof_prevalidated&&!internal::valid_proof(p,v))return out;
    auto tock=Clock::now();out.proof_validation_ms=elapsed(tick,tock);tick=Clock::now();
    VerifierTranscriptReplay replay(p,s,message_points,v);
    std::vector<DoryTarget>g1(p.d+1),g2(p.d+1);std::vector<Fr>rs(p.d),rt(p.d);GT d1=a.X[0],d2=a.X[0];
    out.g1_rexp_challenges.resize(p.d);out.g2_rexp_challenges.resize(p.d);
    double rexp_ms=0;
    for(std::size_t q=0;q<p.d;++q){const auto c1=q? v.g1_rexp_claims[q-1]:a.g1_round0;
      const auto c2=q? v.g2_rexp_claims[q-1]:a.g2_round0;
      const auto [r,rp]=replay.rexp_round(q,c1,c2);
      auto phase=Clock::now();
      rs[p.d-q-1]=r;rt[p.d-q-1]=rp;
      out.g1_rexp_challenges[q]=r;out.g2_rexp_challenges[q]=rp;
      g1[q+1]={product({d1,gpow(c1.E,r),gpow(c1.F,inv(r))}),gmul(c1.TL,gpow(c1.TR,r)),
        gmul(a.delta2L[q],gpow(a.delta2R[q],inv(r)))};
      g2[q+1]={product({d2,gpow(c2.E,rp),gpow(c2.F,inv(rp))}),gmul(a.delta1L[q],gpow(a.delta1R[q],inv(rp))),
        gmul(c2.TL,gpow(c2.TR,rp))};d1=g1[q+1].D1;d2=g2[q+1].D2;
      auto phase_end=Clock::now();rexp_ms+=elapsed(phase,phase_end);}
    replay.finish_rexp(v);
    out.rexp_gt_recurrence_ms=rexp_ms;tick=Clock::now();
    out.g1_rexp_targets=g1;out.g2_rexp_targets=g2;
    const auto sw=tensor_vector(rs),tw=tensor_vector(rt);
    tock=Clock::now();out.tensor_reconstruction_ms=elapsed(tick,tock);
    const auto msm_wall_start=Clock::now();
    if(strategy!=MsmStrategy::Sequential){
      std::exception_ptr worker_error;
      std::jthread g1_worker([&]{
        const auto start=Clock::now();
        try{out.Y_M=msm(message_points,sw);}
        catch(...){worker_error=std::current_exception();}
        out.message_g1_msm_ms=elapsed(start,Clock::now());
      });
      const auto g2_start=Clock::now();
      out.Y_pk=msm(s.public_keys,tw);
      out.public_key_g2_msm_ms=elapsed(g2_start,Clock::now());
      g1_worker.join();
      if(worker_error)std::rethrow_exception(worker_error);
    }else{
      auto start=Clock::now();out.Y_M=msm(message_points,sw);
      out.message_g1_msm_ms=elapsed(start,Clock::now());start=Clock::now();out.Y_pk=msm(s.public_keys,tw);
      out.public_key_g2_msm_ms=elapsed(start,Clock::now());
    }
    out.public_input_msm_wall_ms=elapsed(msm_wall_start,Clock::now());tick=Clock::now();
    const double transcript_before_application=replay.timing_ms();
    DoryTarget im{pair(out.Y_M,p.Lprime),v.cm_M,pair(v.R_Gamma,p.Lprime)};
    DoryTarget ip{pair(p.L,out.Y_pk),pair(p.L,v.R_Lambda),v.cm_pk};
    DoryTarget it{v.T,v.cm_M,v.cm_pk};
    const auto [eta,zeta]=replay.application(im,ip,it,v);
    auto agg=batch_target(im,ip,v.U1,eta);
    agg=batch_target(agg,it,v.U2,zeta);
    out.eta=eta;out.zeta=zeta;out.accumulator_targets.push_back(agg);
    tock=Clock::now();out.application_batching_ms=elapsed(tick,tock)-(replay.timing_ms()-transcript_before_application);tick=Clock::now();
    const double transcript_before_dory=replay.timing_ms();
    out.dory_beta.resize(p.d);out.dory_alpha.resize(p.d);out.insert_g1_gamma.resize(p.d);out.insert_g2_gamma.resize(p.d);
    for(std::size_t j=1;j<=p.d;++j){const auto [beta,alpha]=replay.dory_round(j-1,v.dory_steps[j-1]);
      out.dory_beta[j-1]=beta;out.dory_alpha[j-1]=alpha;
      auto bar=fold_target(a,agg,v.dory_steps[j-1],j-1,beta,alpha);
      out.dory_fold_targets.push_back(bar);
      const Fr gg=replay.insert_g1(j,v.insert_g1_u[j-1]);auto hat=batch_target(bar,g1[j],v.insert_g1_u[j-1],gg);
      const Fr gl=replay.insert_g2(j,v.insert_g2_u[j-1]);agg=batch_target(hat,g2[j],v.insert_g2_u[j-1],gl);
      out.insert_g1_gamma[j-1]=gg;out.insert_g2_gamma[j-1]=gl;out.accumulator_targets.push_back(agg);}
    tock=Clock::now();out.integrated_dory_ms=elapsed(tick,tock)-(replay.timing_ms()-transcript_before_dory);tick=Clock::now();
    replay.bind_final(v);
    out.final_dory=agg;out.final_g1_rexp=g1[p.d];out.final_g2_rexp=g2[p.d];
    tock=Clock::now();out.final_challenges_ms=0;out.transcript_hashing_ms=replay.timing_ms();tick=Clock::now();
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
                                                       MsmStrategy strategy){
  VerificationTrace out;
  using Clock=std::chrono::steady_clock;
  const auto elapsed=[](auto x,auto y){return std::chrono::duration<double,std::milli>(y-x).count();};
  const auto total_start=Clock::now();
  try{
    const auto&p=c.parameters();const auto&a=c.precomputation();const auto&s=c.statement();
    if(!validated_proof_matches(p,validated))return out;
    const auto&v=validated.proof();const auto message_points=c.message_points();
  if(internal::context_binding(p,a,s,message_points)!=c.binding())return out;
    auto tick=Clock::now();out.proof_validation_ms=0;
    VerifierTranscriptReplay replay(p,s,message_points,v);
    SymbolicGTArena arena;
    std::vector<SymbolicTarget>g1(p.d+1),g2(p.d+1);std::vector<Fr>rs(p.d),rt(p.d);
    auto d1=symbolic_atom(arena,a.X[0]),d2=symbolic_atom(arena,a.X[0]);
    out.g1_rexp_challenges.resize(p.d);out.g2_rexp_challenges.resize(p.d);
    if(capture_intermediates){out.g1_rexp_targets.resize(p.d+1);out.g2_rexp_targets.resize(p.d+1);}
    double rexp_ms=0;
    for(std::size_t q=0;q<p.d;++q){const auto c1=q?v.g1_rexp_claims[q-1]:a.g1_round0;
      const auto c2=q?v.g2_rexp_claims[q-1]:a.g2_round0;
      const auto [r,rp]=replay.rexp_round(q,c1,c2);
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
    replay.finish_rexp(v);out.rexp_gt_recurrence_ms=rexp_ms;tick=Clock::now();
    const auto sw=tensor_vector(rs),tw=tensor_vector(rt);auto tock=Clock::now();
    out.tensor_reconstruction_ms=elapsed(tick,tock);
    const auto msm_wall_start=Clock::now();
    if(strategy==MsmStrategy::SplitG2&&s.public_keys.size()>=2){
      const auto half=s.public_keys.size()/2;G2 g2_left,g2_right;
      std::array<std::exception_ptr,3> errors{};std::vector<std::jthread> workers;workers.reserve(2);
      workers.emplace_back([&]{const auto start=Clock::now();try{out.Y_M=msm(message_points,sw);}
          catch(...){errors[0]=std::current_exception();}out.message_g1_msm_ms=elapsed(start,Clock::now());});
      workers.emplace_back([&]{const auto start=Clock::now();
          try{g2_left=msm(std::span(s.public_keys).first(half),std::span(tw).first(half));}
          catch(...){errors[1]=std::current_exception();}
          out.public_key_g2_left_msm_ms=elapsed(start,Clock::now());});
      const auto right_start=Clock::now();
      try{g2_right=msm(std::span(s.public_keys).subspan(half),std::span(tw).subspan(half));}
      catch(...){errors[2]=std::current_exception();}
      out.public_key_g2_right_msm_ms=elapsed(right_start,Clock::now());
      workers.clear();for(const auto& error:errors)if(error)std::rethrow_exception(error);
      out.Y_pk=add(g2_left,g2_right);
      out.public_key_g2_msm_ms=std::max(out.public_key_g2_left_msm_ms,out.public_key_g2_right_msm_ms);
    }else{
      std::exception_ptr worker_error;
      std::jthread g1_worker([&]{const auto start=Clock::now();try{out.Y_M=msm(message_points,sw);}
        catch(...){worker_error=std::current_exception();}out.message_g1_msm_ms=elapsed(start,Clock::now());});
      const auto g2_start=Clock::now();out.Y_pk=msm(s.public_keys,tw);
      out.public_key_g2_msm_ms=elapsed(g2_start,Clock::now());
      g1_worker.join();if(worker_error)std::rethrow_exception(worker_error);
    }
    out.public_input_msm_wall_ms=elapsed(msm_wall_start,Clock::now());tick=Clock::now();
    const double transcript_before_application=replay.timing_ms();
    DoryTarget im{pair(out.Y_M,p.Lprime),v.cm_M,pair(v.R_Gamma,p.Lprime)};
    DoryTarget ip{pair(p.L,out.Y_pk),pair(p.L,v.R_Lambda),v.cm_pk};
    DoryTarget it{v.T,v.cm_M,v.cm_pk};
    SymbolicTarget sim{symbolic_atom(arena,im.D0),symbolic_atom(arena,im.D1),symbolic_atom(arena,im.D2)};
    SymbolicTarget sip{symbolic_atom(arena,ip.D0),symbolic_atom(arena,ip.D1),symbolic_atom(arena,ip.D2)};
    SymbolicTarget sit{symbolic_atom(arena,it.D0),symbolic_atom(arena,it.D1),symbolic_atom(arena,it.D2)};
    const auto [eta,zeta]=replay.application(im,ip,it,v);
    auto agg=symbolic_batch(arena,sim,sip,v.U1,eta);
    agg=symbolic_batch(arena,agg,sit,v.U2,zeta);
    out.eta=eta;out.zeta=zeta;if(capture_intermediates)out.accumulator_targets.push_back(symbolic_evaluate(arena,agg));
    tock=Clock::now();out.application_batching_ms=elapsed(tick,tock)-(replay.timing_ms()-transcript_before_application);tick=Clock::now();
    const double transcript_before_dory=replay.timing_ms();
    out.dory_beta.resize(p.d);out.dory_alpha.resize(p.d);out.insert_g1_gamma.resize(p.d);out.insert_g2_gamma.resize(p.d);
    for(std::size_t j=1;j<=p.d;++j){const auto&step=v.dory_steps[j-1];
      const auto [beta,alpha]=replay.dory_round(j-1,step);
      auto bar=symbolic_fold(arena,a,agg,step,j-1,beta,alpha);
      if(capture_intermediates)out.dory_fold_targets.push_back(symbolic_evaluate(arena,bar));
      const Fr gg=replay.insert_g1(j,v.insert_g1_u[j-1]);
      auto hat=symbolic_batch(arena,bar,g1[j],v.insert_g1_u[j-1],gg);
      const Fr gl=replay.insert_g2(j,v.insert_g2_u[j-1]);
      agg=symbolic_batch(arena,hat,g2[j],v.insert_g2_u[j-1],gl);
      out.dory_beta[j-1]=beta;out.dory_alpha[j-1]=alpha;out.insert_g1_gamma[j-1]=gg;out.insert_g2_gamma[j-1]=gl;
      if(capture_intermediates)out.accumulator_targets.push_back(symbolic_evaluate(arena,agg));
    }
    tock=Clock::now();out.integrated_dory_ms=elapsed(tick,tock)-(replay.timing_ms()-transcript_before_dory);tick=Clock::now();
    replay.bind_final(v);
    out.final_challenges_ms=0;out.transcript_hashing_ms=replay.timing_ms();
    const auto multiexp_start=Clock::now();
    if(capture_intermediates){
      out.final_dory=out.accumulator_targets.back();out.final_g1_rexp=out.g1_rexp_targets.back();
      out.final_g2_rexp=out.g2_rexp_targets.back();
    }else{
      std::array<GT,5> values;std::array<std::exception_ptr,5> errors{};
      const std::array<const SymbolicGT*,5> expressions{
        &agg.D0,&agg.D1,&agg.D2,&g1[p.d].D1,&g2[p.d].D2};
      std::vector<std::jthread> workers;workers.reserve(4);
      for(std::size_t i=0;i<4;++i)workers.emplace_back([&,i]{
          try{values[i]=symbolic_evaluate(arena,*expressions[i]);}
          catch(...){errors[i]=std::current_exception();}
        });
      try{values[4]=symbolic_evaluate(arena,*expressions[4]);}
      catch(...){errors[4]=std::current_exception();}
      workers.clear();
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
  return verify_online_core(c,p,false,MsmStrategy::Sequential);
}
bool verify_online(const ValidatedVerifierContext&c,const Proof&p){return verify_online_diagnostic(c,p).accepted;}
VerificationTrace verify_online_diagnostic(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,false,MsmStrategy::SplitG2);
}
bool verify_online(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_sequential_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  if(!validated_proof_matches(c.parameters(),p))return {};
  return verify_online_core(c,p.proof(),true,MsmStrategy::Sequential);
}
bool verify_online_sequential_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_sequential_msm_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_parallel_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  if(!validated_proof_matches(c.parameters(),p))return {};
  return verify_online_core(c,p.proof(),true,MsmStrategy::Parallel);
}
bool verify_online_parallel_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_parallel_msm_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_symbolic_gt_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,false,MsmStrategy::Parallel);
}
bool verify_online_symbolic_gt(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_gt_diagnostic(c,p).accepted;
}
VerificationTrace verify_online_symbolic_gt_differential_trace(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,true,MsmStrategy::Parallel);
}
VerificationTrace verify_online_split_g2_msm_diagnostic(
    const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_symbolic_core(c,p,false,MsmStrategy::SplitG2);
}
bool verify_online_split_g2_msm(const ValidatedVerifierContext&c,const ValidatedProof&p){
  return verify_online_split_g2_msm_diagnostic(c,p).accepted;
}
bool verify_safe(const PublicParameters&p,const Precomputation&a,const Statement&s,const Proof&v){
  auto c=prepare_verifier_context(p,a,s);return c&&verify_online(*c,v);
}
bool direct_bls_verify(const PublicParameters&p,const Statement&s){
  try{if(!internal::valid_public_parameters(p)||!internal::valid_statement(p,s))return false;return pair(s.sigma_agg,p.H)==direct_pairing_product(hash_messages(p,s),s.public_keys);}
  catch(...){return false;}
}
std::size_t proof_payload_bytes(const Proof&p){
  std::size_t n=0;
  const auto include=[&n](const auto& value){
    std::size_t next{};
    if(!internal::checked_add(n,serialize(value).size(),next))throw std::overflow_error("proof payload size");
    n=next;
  };
  include(p.cm_M);include(p.cm_pk);include(p.T);
  for(const auto&c:p.g1_rexp_claims){include(c.E);include(c.F);include(c.TL);include(c.TR);}
  for(const auto&c:p.g2_rexp_claims){include(c.E);include(c.F);include(c.TL);include(c.TR);}
  include(p.R_Gamma);include(p.R_Lambda);include(p.U1);include(p.U2);
  for(const auto&s:p.dory_steps){include(s.A1L);include(s.A1R);include(s.A2L);include(s.A2R);include(s.W1);include(s.W2);}
  for(const auto&x:p.insert_g1_u)include(x);for(const auto&x:p.insert_g2_u)include(x);
  include(p.Phi_final);include(p.Theta_final);return n;
}
}

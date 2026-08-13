#include "vme_ibf/setup.hpp"
#include "vme_ibf/group_utils.hpp"
#include "internal/crypto.hpp"
#include <cstring>
#include <stdexcept>

namespace vme_ibf {
static void append_raw(Bytes& b,const Bytes& x){b.insert(b.end(),x.begin(),x.end());}
static void append_u64(Bytes& b,std::size_t x){auto e=encode_u64_be(static_cast<std::uint64_t>(x));append_raw(b,e);}
static Fr nonzero(ScalarRng&r){Fr x;do{x=r.random_fr();}while(x.isZero());return x;}
static G1 base_g1(){G1 p;const char*s="VME.BF.G2/BN254/G1-BASE/V2";mcl::bn::hashAndMapToG1(p,s,strlen(s));return p;}
static G2 base_g2(){G2 p;const char*s="VME.BF.G2/BN254/G2-BASE/V2";mcl::bn::hashAndMapToG2(p,s,strlen(s));return p;}
DeterministicRng::DeterministicRng(std::string s):seed_(std::move(s)){}
Fr DeterministicRng::random_fr(){Bytes b;append_frame(b,"VME.BF.G2/DETERMINISTIC-RNG/TEST-ONLY/V2");append_frame(b,seed_);auto c=encode_u64_be(counter_++);b.insert(b.end(),c.begin(),c.end());auto h=sha256(b);Fr x;x.setBigEndianMod(h.data(),h.size());return x;}
Fr SecureRng::random_fr(){Fr x;x.setByCSPRNG();return x;}
Digest compute_crs_digest(const VmeIbfCRS& c){Bytes b;append_frame(b,"VME.BF.G2/CRS/V2");append_frame(b,"BN254/mcl-v3.00");append_frame(b,"Pi_vme.bf:G2;Pi_rexp.bf:G1;e:G1xG2->GT");append_u64(b,c.d);append_u64(b,c.n);append_frame(b,"LEFT-PREFIX-CHAIN/V2");for(auto&p:c.G)append_frame(b,serialize(p));for(auto&p:c.H)append_frame(b,serialize(p));append_frame(b,serialize(c.L));append_frame(b,serialize(c.Lprime));return sha256(b);}
Digest compute_statement_input_digest(const VmeIbfCRS&c,std::span<const Fr>x){Bytes b;append_frame(b,"VME.BF.G2/STATEMENT-INPUT/V2");append_frame(b,c.digest);append_u64(b,c.d);append_u64(b,c.n);for(auto&s:x)append_frame(b,serialize(s));return sha256(b);}
SetupResult setup(std::size_t d,const std::vector<G2>&H,ScalarRng&rng){
 initialize();const size_t n=internal::dimension_size(d);if(H.size()!=n)throw std::invalid_argument("public_H length mismatch");for(auto&p:H)if(!valid_g2(p,true))throw std::invalid_argument("invalid public H");
 SetupResult o;o.crs.d=d;o.crs.n=n;o.crs.H=H;o.crs.G.reserve(n);G2 b2=base_g2();G1 b1=base_g1();for(size_t i=0;i<n;++i)o.crs.G.push_back(g1_pow(b1,nonzero(rng)));o.crs.L=g1_pow(b1,nonzero(rng));o.crs.Lprime=g2_pow(b2,nonzero(rng));if(!valid_g1(o.crs.L,true)||!valid_g2(o.crs.Lprime,true))throw std::runtime_error("invalid generated CRS");for(auto&p:o.crs.G)if(!valid_g1(p,true))throw std::runtime_error("invalid generated auxiliary G");o.crs.digest=compute_crs_digest(o.crs);
 o.precomp.pairing_x.reserve(d+1);for(size_t k=0;k<=d;++k){size_t m=n>>k;o.precomp.pairing_x.push_back(pairing_product(std::span(o.crs.G).first(m),std::span(H).first(m)));}o.precomp.delta1R.reserve(d);o.precomp.delta2R.reserve(d);for(size_t k=0;k<d;++k){size_t m=n>>k,h=m/2;o.precomp.delta1R.push_back(pairing_product(std::span(o.crs.G).subspan(h,h),std::span(H).first(h)));o.precomp.delta2R.push_back(pairing_product(std::span(o.crs.G).first(h),std::span(H).subspan(h,h)));}GT ll;mcl::bn::pairing(ll,o.crs.L,o.crs.Lprime);o.precomp.pairing_LLprime=ll;
 o.statement_input.x.reserve(n);for(size_t i=0;i<n;++i)o.statement_input.x.push_back(rng.random_fr());o.statement_input.digest=compute_statement_input_digest(o.crs,o.statement_input.x);o.prover_input.x=o.statement_input.x;if(!validate_crs(o.crs)||!validate_precomputation_shape(o.crs,o.precomp)||!validate_precomputation_elements(o.precomp)||!audit_precomputation(o.crs,o.precomp))throw std::logic_error("generated setup validation failed");return o;
}
}

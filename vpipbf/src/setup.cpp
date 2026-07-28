#include "vpip_bf/setup.hpp"
#include "vpip_bf/group_utils.hpp"
#include "vpip_bf/transcript.hpp"
#include <cstring>
#include <limits>
#include <stdexcept>

namespace vpip_bf {
namespace {
void append_raw(Bytes& b,const Bytes& x){b.insert(b.end(),x.begin(),x.end());}
void append_u64(Bytes& b,std::size_t x){append_raw(b,encode_u64_be(static_cast<std::uint64_t>(x)));}
Fr nonzero(ScalarRng&r){Fr x;do{x=r.random_fr();}while(x.isZero());return x;}
G1 base_g1(){G1 p;const char*s="vpipbf/fs/v1/BN254/G1";mcl::bn::hashAndMapToG1(p,s,strlen(s));return p;}
G2 base_g2(){G2 p;const char*s="vpipbf/fs/v1/BN254/G2";mcl::bn::hashAndMapToG2(p,s,strlen(s));return p;}
}

DeterministicRng::DeterministicRng(std::string s):seed_(std::move(s)){}
Fr DeterministicRng::random_fr(){Bytes b;append_frame(b,"vpipbf/test-rng/v1");append_frame(b,seed_);append_raw(b,encode_u64_be(counter_++));auto h=sha256(b);Fr x;x.setBigEndianMod(h.data(),h.size());return x;}
Fr SecureRng::random_fr(){Fr x;x.setByCSPRNG();return x;}

Digest compute_crs_digest(const VpipBfCRS&c){Bytes b;append_frame(b,"vpipbf/crs/v1");append_frame(b,"BN254/mcl-v3.00");append_frame(b,"Pi_vpip.bf;Pi_rexp.bf:G1;e:G1xG2->GT");append_u64(b,c.d);append_u64(b,c.n);for(auto&p:c.G)append_frame(b,serialize(p));for(auto&p:c.H)append_frame(b,serialize(p));append_frame(b,serialize(c.Lprime));return sha256(b);}
Digest compute_precomputation_digest(const VpipBfCRS&c,const VpipBfPrecomputation&p){Bytes b;append_frame(b,"vpipbf/precomputation/v1");append_frame(b,c.digest);for(auto&x:p.pairing_x)append_frame(b,serialize(x));for(auto&x:p.delta1R)append_frame(b,serialize(x));for(auto&x:p.delta2R)append_frame(b,serialize(x));return sha256(b);}
Digest compute_statement_input_digest(const VpipBfCRS&c,std::span<const G1>X){Bytes b;append_frame(b,"vpipbf/statement-input/v1");append_frame(b,c.digest);for(auto&x:X)append_frame(b,serialize(x));return sha256(b);}

SetupResult setup_parameters(std::size_t d,const std::vector<G2>&lambda,const std::vector<G1>&X,ScalarRng&rng){
 initialize();if(d<1||d>=std::numeric_limits<std::size_t>::digits)throw std::invalid_argument("invalid d");size_t n=size_t{1}<<d;if(lambda.size()!=n||X.size()!=n)throw std::invalid_argument("public vector length mismatch");for(auto&p:lambda)if(!valid_g2(p,true))throw std::invalid_argument("invalid Lambda");for(auto&p:X)if(!valid_g1(p))throw std::invalid_argument("invalid X");
 SetupResult o;o.crs.d=d;o.crs.n=n;o.crs.H=lambda;o.crs.G.reserve(n);G1 b1=base_g1();for(size_t i=0;i<n;++i)o.crs.G.push_back(g1_pow(b1,nonzero(rng)));o.crs.Lprime=g2_pow(base_g2(),nonzero(rng));o.crs.digest=compute_crs_digest(o.crs);
 o.statement_input.X=X;o.statement_input.digest=compute_statement_input_digest(o.crs,X);o.prover_input.X=X;return o;
}
VpipBfPrecomputation precompute(const VpipBfCRS&c){if(c.d<1||c.n!=(size_t{1}<<c.d)||c.G.size()!=c.n||c.H.size()!=c.n||compute_crs_digest(c)!=c.digest)throw std::invalid_argument("invalid CRS");VpipBfPrecomputation p;for(size_t k=0;k<=c.d;++k){size_t m=c.n>>k;p.pairing_x.push_back(pairing_product(std::span(c.G).first(m),std::span(c.H).first(m)));}for(size_t k=0;k<c.d;++k){size_t m=c.n>>k,h=m/2;p.delta1R.push_back(pairing_product(std::span(c.G).subspan(h,h),std::span(c.H).first(h)));p.delta2R.push_back(pairing_product(std::span(c.G).first(h),std::span(c.H).subspan(h,h)));}p.digest=compute_precomputation_digest(c,p);return p;}
SetupResult setup(std::size_t d,const std::vector<G2>&lambda,const std::vector<G1>&X,ScalarRng&rng){auto o=setup_parameters(d,lambda,X,rng);o.precomp=precompute(o.crs);return o;}
SetupResult setup(std::size_t d,const std::vector<G2>&lambda,ScalarRng&rng){initialize();if(d>=std::numeric_limits<size_t>::digits)throw std::invalid_argument("invalid d");size_t n=size_t{1}<<d;G1 b=base_g1();std::vector<G1>X;X.reserve(n);for(size_t i=0;i<n;++i)X.push_back(g1_pow(b,nonzero(rng)));return setup(d,lambda,X,rng);}
}

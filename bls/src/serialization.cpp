#include "blsagg/serialization.hpp"
#include "blsagg/transcript.hpp"
#include "internal/crypto.hpp"
#include "internal/validation.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>

namespace blsagg {
namespace {
constexpr std::uint16_t kVersion=1;
constexpr std::size_t kHeaderBytes=24;

void seterr(DecodeError*e,DecodeError x){if(e)*e=x;}
void raw(Bytes&o,std::span<const std::uint8_t>x){internal::append_raw(o,x);}
void u16(Bytes&o,std::uint16_t x){o.push_back(std::uint8_t(x>>8));o.push_back(std::uint8_t(x));}
void u32(Bytes&o,std::uint32_t x){for(int s=24;s>=0;s-=8)o.push_back(std::uint8_t(x>>s));}
std::uint64_t as_u64(std::size_t x){
  if constexpr(sizeof(std::size_t)>sizeof(std::uint64_t))
    if(x>std::numeric_limits<std::uint64_t>::max())throw std::overflow_error("wire length");
  return static_cast<std::uint64_t>(x);
}
void u64(Bytes&o,std::uint64_t x){raw(o,encode_u64(x));}
void frame(Bytes&o,std::span<const std::uint8_t>x){u64(o,as_u64(x.size()));raw(o,x);}
template<class T>void element(Bytes&o,const T&x){auto b=serialize(x);frame(o,b);}
std::uint8_t mode_byte(AggregationMode m){return static_cast<std::uint8_t>(m);}

Bytes wrap(const char magic[9],AggregationMode mode,std::size_t d,const Bytes&body){
  if(d>std::numeric_limits<std::uint32_t>::max())throw std::overflow_error("dimension");
  std::size_t wire_size{};
  if(!internal::checked_add(kHeaderBytes,body.size(),wire_size))throw std::overflow_error("wire length");
  Bytes o;o.reserve(wire_size);
  o.insert(o.end(),magic,magic+8);u16(o,kVersion);o.push_back(mode_byte(mode));o.push_back(0);
  u32(o,static_cast<std::uint32_t>(d));u64(o,as_u64(body.size()));raw(o,body);return o;
}

class Reader{
 public:
  Reader(std::span<const std::uint8_t>b,DecodeError*e):b_(b),e_(e){}
  bool bytes(std::size_t n,std::span<const std::uint8_t>&x){
    if(pos_>b_.size()||n>b_.size()-pos_){fail(DecodeError::Truncated);return false;}x=b_.subspan(pos_,n);pos_+=n;return true;
  }
  bool u64v(std::uint64_t&x){std::span<const std::uint8_t>q;if(!bytes(8,q))return false;x=0;for(auto c:q)x=(x<<8)|c;return true;}
  bool count(std::size_t&x){
    std::uint64_t n;if(!u64v(n))return false;
    if(n>std::numeric_limits<std::size_t>::max()){fail(DecodeError::IntegerOverflow);return false;}
    x=static_cast<std::size_t>(n);
    if(x>(remaining()/8)+1){fail(DecodeError::InvalidLength);return false;}
    return true;
  }
  bool framed(std::span<const std::uint8_t>&x){
    std::uint64_t n;if(!u64v(n))return false;
    if(n>std::numeric_limits<std::size_t>::max()){fail(DecodeError::IntegerOverflow);return false;}
    return bytes(static_cast<std::size_t>(n),x);
  }
  std::size_t remaining()const{return pos_<=b_.size()?b_.size()-pos_:0;}
  bool done(){if(pos_!=b_.size()){fail(DecodeError::TrailingBytes);return false;}return true;}
  void fail(DecodeError x){if(ok_){ok_=false;seterr(e_,x);}}
 private:
  std::span<const std::uint8_t>b_;std::size_t pos_{};DecodeError*e_{};bool ok_{true};
};

bool header(std::span<const std::uint8_t>b,const char magic[9],AggregationMode expected,
            bool require_mode,std::size_t&d,std::span<const std::uint8_t>&body,DecodeError*e){
  seterr(e,DecodeError::None);
  if(b.size()<kHeaderBytes){seterr(e,DecodeError::Truncated);return false;}
  if(!std::equal(b.begin(),b.begin()+8,reinterpret_cast<const std::uint8_t*>(magic))){seterr(e,DecodeError::WrongMagic);return false;}
  const auto ver=std::uint16_t(std::uint16_t(b[8])<<8|b[9]);if(ver!=kVersion){seterr(e,DecodeError::WrongVersion);return false;}
  const auto mb=b[10];if((mb!=1&&mb!=2)||(require_mode&&mb!=mode_byte(expected))){seterr(e,DecodeError::WrongMode);return false;}
  if(b[11]!=0){seterr(e,DecodeError::WrongVersion);return false;}
  d=(std::uint32_t(b[12])<<24)|(std::uint32_t(b[13])<<16)|(std::uint32_t(b[14])<<8)|b[15];
  if(d<1||d>=std::numeric_limits<std::size_t>::digits){seterr(e,DecodeError::InvalidDimension);return false;}
  std::uint64_t n=0;for(std::size_t i=16;i<24;++i)n=(n<<8)|b[i];
  if(n>std::numeric_limits<std::size_t>::max()){seterr(e,DecodeError::IntegerOverflow);return false;}
  const auto want=static_cast<std::size_t>(n);
  if(want>b.size()-kHeaderBytes){seterr(e,DecodeError::Truncated);return false;}
  if(want<b.size()-kHeaderBytes){seterr(e,DecodeError::TrailingBytes);return false;}
  body=b.subspan(kHeaderBytes,want);return true;
}

template<class T>bool decode_element(Reader&r,T&out,DecodeError bad,bool nonidentity){
  std::span<const std::uint8_t>b;if(!r.framed(b))return false;T x;
  if(x.deserialize(b.data(),b.size())!=b.size()){r.fail(bad);return false;}
  if(serialize(x)!=Bytes(b.begin(),b.end())){r.fail(DecodeError::NonCanonical);return false;}
  bool good=true;
  if constexpr(std::is_same_v<T,G1>||std::is_same_v<T,G2>)good=x.isValid()&&x.isValidOrder();
  else good=mcl::bn::isValidGT(x);
  if(!good){r.fail(bad);return false;}
  if(nonidentity&&x.isZero()){r.fail(DecodeError::IdentityNotAllowed);return false;}
  out=x;return true;
}
bool digest(Reader&r,Digest&d){std::span<const std::uint8_t>x;if(!r.framed(x))return false;if(x.size()!=d.size()){r.fail(DecodeError::InvalidLength);return false;}std::copy(x.begin(),x.end(),d.begin());return true;}
template<class T>void vec_elements(Bytes&b,const std::vector<T>&v){u64(b,v.size());for(const auto&x:v)element(b,x);}
template<class T>bool read_vec(Reader&r,std::vector<T>&v,DecodeError bad,bool nz,std::optional<std::size_t>want={}){
  std::size_t n;if(!r.count(n))return false;if(want&&n!=*want){r.fail(DecodeError::InvalidShape);return false;}
  v.clear();v.reserve(n);for(std::size_t i=0;i<n;++i){T x;if(!decode_element(r,x,bad,nz))return false;v.push_back(x);}return true;
}
void claims(Bytes&b,const RexpClaims&c){element(b,c.E);element(b,c.F);element(b,c.TL);element(b,c.TR);}
bool read_claims(Reader&r,RexpClaims&c){return decode_element(r,c.E,DecodeError::InvalidGT,false)&&decode_element(r,c.F,DecodeError::InvalidGT,false)&&decode_element(r,c.TL,DecodeError::InvalidGT,false)&&decode_element(r,c.TR,DecodeError::InvalidGT,false);}
void target_vec(Bytes&b,const std::vector<GT>&v){vec_elements(b,v);}
bool mode_valid(AggregationMode m){return m==AggregationMode::BasicDistinct||m==AggregationMode::Augmented;}
}

ValidatedProof::ValidatedProof(Proof proof,Digest parameter_digest,Digest wire_binding)
    :proof_(std::move(proof)),parameter_digest_(parameter_digest),wire_binding_(wire_binding){}

Bytes serialize_public_parameters(const PublicParameters&p){
  if(!internal::valid_public_parameters(p))throw std::invalid_argument("invalid public parameters");
  Bytes b;element(b,p.H);vec_elements(b,p.Gamma);vec_elements(b,p.Lambda);element(b,p.L);element(b,p.Lprime);frame(b,p.digest);
  return wrap("BLSAPP01",p.mode,p.d,b);
}
bool deserialize_public_parameters(std::span<const std::uint8_t>in,PublicParameters&o,DecodeError*e){
  o={};seterr(e,DecodeError::None);
  try{
    std::size_t d=0;std::span<const std::uint8_t>body;if(!header(in,"BLSAPP01",AggregationMode::BasicDistinct,false,d,body,e))return false;
    const auto mode=static_cast<AggregationMode>(in[10]);if(!mode_valid(mode)){seterr(e,DecodeError::WrongMode);return false;}
    const std::size_t k=std::size_t{1}<<d;Reader r(body,e);PublicParameters p;p.d=d;p.k=k;p.mode=mode;
    if(!decode_element(r,p.H,DecodeError::InvalidG2,true)||!read_vec(r,p.Gamma,DecodeError::InvalidG1,true,k)||
       !read_vec(r,p.Lambda,DecodeError::InvalidG2,true,k)||!decode_element(r,p.L,DecodeError::InvalidG1,true)||
       !decode_element(r,p.Lprime,DecodeError::InvalidG2,true)||!digest(r,p.digest)||!r.done())return false;
    if(!internal::valid_public_parameters(p)){seterr(e,DecodeError::InvalidDigest);return false;}o=std::move(p);return true;
  }catch(...){seterr(e,DecodeError::InvalidShape);return false;}
}

Bytes serialize_precomputation(const PublicParameters&p,const Precomputation&a){
  if(!internal::valid_precomputation(p,a))throw std::invalid_argument("invalid precomputation");
  Bytes b;frame(b,p.digest);frame(b,a.digest);u64(b,a.gamma_chain.size());for(const auto&v:a.gamma_chain)vec_elements(b,v);
  u64(b,a.lambda_chain.size());for(const auto&v:a.lambda_chain)vec_elements(b,v);
  target_vec(b,a.X);target_vec(b,a.delta1L);target_vec(b,a.delta1R);target_vec(b,a.delta2L);target_vec(b,a.delta2R);
  claims(b,a.g1_round0);claims(b,a.g2_round0);return wrap("BLSAUX01",p.mode,p.d,b);
}
bool deserialize_precomputation(std::span<const std::uint8_t>in,const PublicParameters&p,Precomputation&o,DecodeError*e){
  o={};seterr(e,DecodeError::None);
  try{
    if(!internal::valid_public_parameters(p)){seterr(e,DecodeError::InvalidDigest);return false;}
    std::size_t d=0;std::span<const std::uint8_t>body;if(!header(in,"BLSAUX01",p.mode,true,d,body,e))return false;
    if(d!=p.d){seterr(e,DecodeError::InvalidDimension);return false;}
    Reader r(body,e);Digest pd{};Precomputation a;if(!digest(r,pd))return false;if(pd!=p.digest){seterr(e,DecodeError::InvalidDigest);return false;}if(!digest(r,a.digest))return false;
    std::size_t levels;if(!r.count(levels))return false;if(levels!=p.d+1){seterr(e,DecodeError::InvalidShape);return false;}
    for(std::size_t j=0;j<levels;++j){std::vector<G1>v;if(!read_vec(r,v,DecodeError::InvalidG1,true,p.k>>j))return false;a.gamma_chain.push_back(std::move(v));}
    if(!r.count(levels))return false;if(levels!=p.d+1){seterr(e,DecodeError::InvalidShape);return false;}
    for(std::size_t j=0;j<levels;++j){std::vector<G2>v;if(!read_vec(r,v,DecodeError::InvalidG2,true,p.k>>j))return false;a.lambda_chain.push_back(std::move(v));}
    if(!read_vec(r,a.X,DecodeError::InvalidGT,false,p.d+1)||!read_vec(r,a.delta1L,DecodeError::InvalidGT,false,p.d)||
       !read_vec(r,a.delta1R,DecodeError::InvalidGT,false,p.d)||!read_vec(r,a.delta2L,DecodeError::InvalidGT,false,p.d)||
       !read_vec(r,a.delta2R,DecodeError::InvalidGT,false,p.d)||!read_claims(r,a.g1_round0)||!read_claims(r,a.g2_round0)||!r.done())return false;
    if(!internal::valid_precomputation(p,a)){seterr(e,DecodeError::InvalidDigest);return false;}
    o=std::move(a);return true;
  }catch(...){seterr(e,DecodeError::InvalidShape);return false;}
}

Bytes serialize_statement(const PublicParameters&p,const Statement&s){
  if(!internal::valid_public_parameters(p)||!internal::valid_statement(p,s))throw std::invalid_argument("invalid statement");
  Bytes b;element(b,s.sigma_agg);u64(b,s.messages.size());for(const auto&m:s.messages)frame(b,m);vec_elements(b,s.public_keys);
  return wrap("BLSAST01",p.mode,p.d,b);
}
bool deserialize_statement(std::span<const std::uint8_t>in,const PublicParameters&p,Statement&o,DecodeError*e){
  o={};seterr(e,DecodeError::None);
  try{
    if(!internal::valid_public_parameters(p)){seterr(e,DecodeError::InvalidDigest);return false;}
    std::size_t d=0;std::span<const std::uint8_t>body;if(!header(in,"BLSAST01",p.mode,true,d,body,e))return false;
    if(d!=p.d){seterr(e,DecodeError::InvalidDimension);return false;}
    Reader r(body,e);Statement s;if(!decode_element(r,s.sigma_agg,DecodeError::InvalidG1,true))return false;
    std::size_t n;if(!r.count(n))return false;if(n!=p.k){seterr(e,DecodeError::InvalidShape);return false;}s.messages.reserve(n);
    for(std::size_t i=0;i<n;++i){std::span<const std::uint8_t>x;if(!r.framed(x))return false;s.messages.emplace_back(x.begin(),x.end());}
    if(!read_vec(r,s.public_keys,DecodeError::InvalidG2,true,p.k)||!r.done())return false;
    if(!internal::valid_statement(p,s)){seterr(e,DecodeError::InvalidShape);return false;}
    o=std::move(s);return true;
  }catch(...){seterr(e,DecodeError::InvalidShape);return false;}
}

Bytes serialize_proof(const PublicParameters&p,const Proof&v){
  if(!internal::valid_public_parameters(p)||!internal::valid_proof(p,v))throw std::invalid_argument("invalid proof");
  Bytes b;element(b,v.cm_M);element(b,v.cm_pk);element(b,v.T);
  u64(b,v.g1_rexp_claims.size());for(const auto&x:v.g1_rexp_claims)claims(b,x);
  u64(b,v.g2_rexp_claims.size());for(const auto&x:v.g2_rexp_claims)claims(b,x);
  element(b,v.R_Gamma);element(b,v.R_Lambda);element(b,v.U1);element(b,v.U2);
  u64(b,v.dory_steps.size());for(const auto&x:v.dory_steps){element(b,x.A1L);element(b,x.A1R);element(b,x.A2L);element(b,x.A2R);element(b,x.W1);element(b,x.W2);}
  vec_elements(b,v.insert_g1_u);vec_elements(b,v.insert_g2_u);element(b,v.Phi_final);element(b,v.Theta_final);
  return wrap("BLSAPF01",p.mode,p.d,b);
}
bool deserialize_proof(std::span<const std::uint8_t>in,const PublicParameters&p,Proof&o,DecodeError*e){
  o={};seterr(e,DecodeError::None);
  try{
    if(!internal::valid_public_parameters(p)){seterr(e,DecodeError::InvalidDigest);return false;}
    std::size_t d=0;std::span<const std::uint8_t>body;if(!header(in,"BLSAPF01",p.mode,true,d,body,e))return false;
    if(d!=p.d){seterr(e,DecodeError::InvalidDimension);return false;}
    Reader r(body,e);Proof v;if(!decode_element(r,v.cm_M,DecodeError::InvalidGT,false)||!decode_element(r,v.cm_pk,DecodeError::InvalidGT,false)||!decode_element(r,v.T,DecodeError::InvalidGT,false))return false;
    std::size_t n;if(!r.count(n))return false;if(n!=p.d-1){seterr(e,DecodeError::InvalidShape);return false;}for(std::size_t i=0;i<n;++i){RexpClaims c;if(!read_claims(r,c))return false;v.g1_rexp_claims.push_back(c);}
    if(!r.count(n))return false;if(n!=p.d-1){seterr(e,DecodeError::InvalidShape);return false;}for(std::size_t i=0;i<n;++i){RexpClaims c;if(!read_claims(r,c))return false;v.g2_rexp_claims.push_back(c);}
    if(!decode_element(r,v.R_Gamma,DecodeError::InvalidG1,false)||!decode_element(r,v.R_Lambda,DecodeError::InvalidG2,false)||
       !decode_element(r,v.U1,DecodeError::InvalidGT,false)||!decode_element(r,v.U2,DecodeError::InvalidGT,false))return false;
    if(!r.count(n))return false;if(n!=p.d){seterr(e,DecodeError::InvalidShape);return false;}for(std::size_t i=0;i<n;++i){DoryStep x;if(!decode_element(r,x.A1L,DecodeError::InvalidGT,false)||!decode_element(r,x.A1R,DecodeError::InvalidGT,false)||!decode_element(r,x.A2L,DecodeError::InvalidGT,false)||!decode_element(r,x.A2R,DecodeError::InvalidGT,false)||!decode_element(r,x.W1,DecodeError::InvalidGT,false)||!decode_element(r,x.W2,DecodeError::InvalidGT,false))return false;v.dory_steps.push_back(x);}
    if(!read_vec(r,v.insert_g1_u,DecodeError::InvalidGT,false,p.d)||!read_vec(r,v.insert_g2_u,DecodeError::InvalidGT,false,p.d)||
       !decode_element(r,v.Phi_final,DecodeError::InvalidG1,false)||!decode_element(r,v.Theta_final,DecodeError::InvalidG2,false)||!r.done())return false;
    if(!internal::valid_proof(p,v)){seterr(e,DecodeError::InvalidShape);return false;}o=std::move(v);return true;
  }catch(...){seterr(e,DecodeError::InvalidShape);return false;}
}

std::optional<ValidatedProof> deserialize_and_validate_proof(
    std::span<const std::uint8_t>in,const PublicParameters&p,DecodeError*e){
  Proof proof;if(!deserialize_proof(in,p,proof,e))return std::nullopt;
  const auto binding=internal::validated_proof_binding(p,proof);
  return ValidatedProof(std::move(proof),p.digest,binding);
}

std::size_t proof_mathematical_payload_bytes(const Proof&p){return proof_payload_bytes(p);}
std::size_t proof_wire_bytes(const PublicParameters&p,const Proof&v){return serialize_proof(p,v).size();}
}

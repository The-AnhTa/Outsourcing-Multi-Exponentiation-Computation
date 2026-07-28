#include "vme_ibf/group_utils.hpp"
#include <mcl/fp.hpp>
#include <iomanip>
#include <mutex>
#include <sstream>
#include <stdexcept>

namespace vme_ibf {
void initialize(){ static std::once_flag once; std::call_once(once,[]{mcl::bn::initPairing(mcl::BN254);}); }
Bytes encode_u64_be(std::uint64_t v){ Bytes b; for(int s=56;s>=0;s-=8)b.push_back(static_cast<std::uint8_t>(v>>s)); return b; }
void append_frame(Bytes& out,std::span<const std::uint8_t> f){ auto n=encode_u64_be(f.size());out.insert(out.end(),n.begin(),n.end());out.insert(out.end(),f.begin(),f.end()); }
void append_frame(Bytes& out,std::string_view f){append_frame(out,{reinterpret_cast<const std::uint8_t*>(f.data()),f.size()});}
Digest sha256(std::span<const std::uint8_t> in){if(in.size()>UINT32_MAX)throw std::length_error("SHA-256 input too large");Digest d{};auto n=mcl::fp::sha256(d.data(),static_cast<std::uint32_t>(d.size()),in.data(),static_cast<std::uint32_t>(in.size()));if(n!=d.size())throw std::runtime_error("SHA-256 failed");return d;}
template<class T> static Bytes enc(const T& v){Bytes b(2048);auto n=v.serialize(b.data(),b.size());if(!n)throw std::runtime_error("mcl serialization failed");b.resize(n);return b;}
Bytes serialize(const Fr& v){return enc(v);} Bytes serialize(const G1& v){return enc(v);} Bytes serialize(const G2& v){return enc(v);} Bytes serialize(const GT& v){return enc(v);}
bool valid_g1(const G1& p,bool nz){return p.isValid()&&p.isValidOrder()&&(!nz||!p.isZero());} bool valid_g2(const G2& p,bool nz){return p.isValid()&&p.isValidOrder()&&(!nz||!p.isZero());}
GT pairing_product(std::span<const G1>a,std::span<const G2>b){if(a.size()!=b.size())throw std::invalid_argument("pairing length mismatch");GT z;if(a.empty()){z.setOne();return z;}GT ml;mcl::bn::millerLoopVec(ml,a.data(),b.data(),a.size(),true);mcl::bn::finalExp(z,ml);return z;}
G1 g1_add(const G1&a,const G1&b){G1 z;G1::add(z,a,b);return z;} G2 g2_add(const G2&a,const G2&b){G2 z;G2::add(z,a,b);return z;}
G1 g1_pow(const G1&p,const Fr&s){G1 z;G1::mul(z,p,s);return z;} G2 g2_pow(const G2&p,const Fr&s){G2 z;G2::mul(z,p,s);return z;}
GT gt_mul(const GT&a,const GT&b){GT z;GT::mul(z,a,b);return z;} GT gt_pow(const GT&a,const Fr&s){GT z;GT::pow(z,a,s);return z;}
G1 g1_multiexp_reference(std::span<const G1>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("MSM length mismatch");G1 z;z.clear();for(size_t i=0;i<a.size();++i)z=g1_add(z,g1_pow(a[i],s[i]));return z;}
G1 g1_multiexp(std::span<const G1>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("MSM length mismatch");G1 z;z.clear();if(!a.empty()){std::vector<G1> work(a.begin(),a.end());G1::mulVec(z,work.data(),s.data(),a.size());}return z;}
G2 g2_multiexp_reference(std::span<const G2>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("MSM length mismatch");G2 z;z.clear();for(size_t i=0;i<a.size();++i)z=g2_add(z,g2_pow(a[i],s[i]));return z;}
G2 g2_multiexp(std::span<const G2>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("MSM length mismatch");G2 z;z.clear();if(!a.empty()){std::vector<G2> work(a.begin(),a.end());G2::mulVec(z,work.data(),s.data(),a.size());}return z;}
G2 g2_multiexp_protocol(std::span<const G2>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("MSM length mismatch");if(a.empty()){G2 z;z.clear();return z;}std::vector<G2>level;level.reserve(a.size());for(size_t i=0;i<a.size();++i)level.push_back(g2_pow(a[i],s[i]));while(level.size()>1){std::vector<G2>next;next.reserve((level.size()+1)/2);size_t i=0;for(;i+1<level.size();i+=2)next.push_back(g2_add(level[i],level[i+1]));if(i<level.size())next.push_back(level[i]);level=std::move(next);}return level[0];}
template<class G>static std::vector<G> cm(std::span<const G>a,std::span<const G>b){if(a.size()!=b.size())throw std::invalid_argument("component length mismatch");std::vector<G>z;z.reserve(a.size());for(size_t i=0;i<a.size();++i){G q;G::add(q,a[i],b[i]);z.push_back(q);}return z;}
std::vector<G1> component_mul(std::span<const G1>a,std::span<const G1>b){return cm(a,b);} std::vector<G2> component_mul(std::span<const G2>a,std::span<const G2>b){return cm(a,b);}
template<class G>static std::vector<G> cp(std::span<const G>a,std::span<const Fr>s){if(a.size()!=s.size())throw std::invalid_argument("component length mismatch");std::vector<G>z;z.reserve(a.size());for(size_t i=0;i<a.size();++i){G q;G::mul(q,a[i],s[i]);z.push_back(q);}return z;}
std::vector<G1> component_pow(std::span<const G1>a,std::span<const Fr>s){return cp(a,s);} std::vector<G2> component_pow(std::span<const G2>a,std::span<const Fr>s){return cp(a,s);}
Fr inner_product(std::span<const Fr>a,std::span<const Fr>b){if(a.size()!=b.size())throw std::invalid_argument("inner product length mismatch");Fr z;z.clear();for(size_t i=0;i<a.size();++i){Fr q;Fr::mul(q,a[i],b[i]);Fr::add(z,z,q);}return z;}
Fr inverse_nonzero(const Fr&a){if(a.isZero())throw std::invalid_argument("cannot invert zero");Fr z;Fr::inv(z,a);return z;}
std::string hex(std::span<const std::uint8_t>b){std::ostringstream o;o<<std::hex<<std::setfill('0');for(auto x:b)o<<std::setw(2)<<unsigned(x);return o.str();}
}

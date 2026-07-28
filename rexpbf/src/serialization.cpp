#include "rexpbf/serialization.hpp"
#include "rexpbf/pairing.hpp"
#include "rexpbf/prove.hpp"
#include "rexpbf/setup.hpp"
#include <cybozu/sha2.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

namespace rexpbf {
namespace {
using Bytes = std::vector<std::uint8_t>;
void append(Bytes& out, std::span<const std::uint8_t> in) { out.insert(out.end(), in.begin(), in.end()); }
void frame_text(Bytes& out, std::string_view s) {
    append_frame(out, {reinterpret_cast<const std::uint8_t*>(s.data()), s.size()});
}
template<class T> Bytes serialize(const T& value) {
    Bytes out(2048);
    const auto n = value.serialize(out.data(), out.size());
    if (!n) throw std::runtime_error("mcl binary serialization failed");
    out.resize(n);
    return out;
}
template<class T> T deserialize(std::span<const std::uint8_t> bytes) {
    if (bytes.empty()) throw std::invalid_argument("empty mcl encoding");
    T value;
    const auto consumed = value.deserialize(bytes.data(), bytes.size());
    if (consumed != bytes.size()) throw std::invalid_argument("malformed or trailing mcl encoding");
    if (serialize(value) != Bytes(bytes.begin(), bytes.end()))
        throw std::invalid_argument("non-canonical mcl encoding");
    return value;
}
void framed(Bytes& out, const Bytes& bytes) { append_frame(out, bytes); }
class Reader {
public:
    explicit Reader(std::span<const std::uint8_t> bytes) : bytes_(bytes) {}
    std::uint64_t u64() {
        if (bytes_.size() - pos_ < 8) throw std::invalid_argument("truncated u64");
        std::uint64_t v=0;for(int i=0;i<8;++i)v=(v<<8)|bytes_[pos_++];return v;
    }
    std::span<const std::uint8_t> frame() {
        const auto n=u64();if(n>bytes_.size()-pos_)throw std::invalid_argument("truncated frame");
        auto out=bytes_.subspan(pos_,static_cast<std::size_t>(n));pos_+=static_cast<std::size_t>(n);return out;
    }
    void expect(std::string_view value) {
        auto got=frame();if(got.size()!=value.size()||!std::equal(got.begin(),got.end(),value.begin()))
            throw std::invalid_argument("wire domain mismatch");
    }
    template<class T> T element(){return deserialize<T>(frame());}
    Digest32 digest(){auto f=frame();if(f.size()!=32)throw std::invalid_argument("digest width mismatch");Digest32 d{};std::copy(f.begin(),f.end(),d.begin());return d;}
    void finish(){if(pos_!=bytes_.size())throw std::invalid_argument("trailing wire bytes");}
private: std::span<const std::uint8_t> bytes_;std::size_t pos_{};
};
template<class P> void require_wire_point(const P&p,std::string_view name){
    if(p.isZero()||!p.isValid()||!p.isValidOrder())throw std::invalid_argument(std::string("invalid ")+std::string(name));
}
}

std::array<std::uint8_t, 8> encode_u64_be(std::uint64_t value) {
    std::array<std::uint8_t, 8> out{};
    for (std::size_t i = 0; i < 8; ++i) out[i] = static_cast<std::uint8_t>(value >> (56 - 8 * i));
    return out;
}
void append_frame(Bytes& out, std::span<const std::uint8_t> bytes) {
    if (bytes.size() > std::numeric_limits<std::uint64_t>::max())
        throw std::length_error("frame is too large");
    const auto size = encode_u64_be(static_cast<std::uint64_t>(bytes.size()));
    append(out, size);
    append(out, bytes);
}
Bytes serialize_g1(const G1& p) { return serialize(p); }
Bytes serialize_g2(const G2& p) { return serialize(p); }
Bytes serialize_gt(const GT& p) { return serialize(p); }
Bytes serialize_fr(const Fr& p) { return serialize(p); }
Bytes serialize_crs(const CRS& crs) {
    initialize_bn254();
    Bytes out;
    append(out, encode_u64_be(crs.d));
    append(out, encode_u64_be(crs.n));
    append(out, encode_u64_be(crs.gamma.size()));
    for (const auto& point : crs.gamma) append(out, serialize_g1(point));
    append(out, encode_u64_be(crs.lambda.size()));
    for (const auto& point : crs.lambda) append(out, serialize_g2(point));
    append(out, crs.digest);
    return out;
}
Bytes serialize_crs_wire(const CRS& crs) {
    if(!validate_crs(crs))throw std::invalid_argument("cannot serialize invalid CRS");Bytes out;
    frame_text(out,"REXP-BF-G1-CRS-WIRE-BN254-V1");frame_text(out,"BN254");frame_text(out,"G1");
    append(out,encode_u64_be(crs.d));append(out,encode_u64_be(crs.n));append(out,encode_u64_be(crs.gamma.size()));
    for(const auto&p:crs.gamma)framed(out,serialize_g1(p));append(out,encode_u64_be(crs.lambda.size()));
    for(const auto&p:crs.lambda)framed(out,serialize_g2(p));append_frame(out,crs.digest);return out;
}
Bytes serialize_statement_wire(const CRS&crs,const Statement&s){
    if(!validate_statement_shape(crs,s))throw std::invalid_argument("cannot serialize invalid statement");Bytes out;
    frame_text(out,"REXP-BF-G1-STATEMENT-WIRE-BN254-V1");frame_text(out,"BN254");frame_text(out,"G1");
    append(out,encode_u64_be(crs.d));append(out,encode_u64_be(crs.n));append_frame(out,s.crs_digest);append(out,encode_u64_be(s.h.size()));
    for(const auto&p:s.h)framed(out,serialize_g1(p));for(const auto*x:{&s.d1_initial,&s.e0,&s.f0,&s.t_left0,&s.t_right0})framed(out,serialize_gt(*x));
    append_frame(out,s.digest);append(out,encode_u64_be(s.pairing_terms));return out;
}
G1 deserialize_g1(std::span<const std::uint8_t> b) { return deserialize<G1>(b); }
G2 deserialize_g2(std::span<const std::uint8_t> b) { return deserialize<G2>(b); }
GT deserialize_gt(std::span<const std::uint8_t> b) { return deserialize<GT>(b); }
CRS deserialize_crs_wire(std::span<const std::uint8_t>b){initialize_bn254();Reader r(b);r.expect("REXP-BF-G1-CRS-WIRE-BN254-V1");r.expect("BN254");r.expect("G1");CRS c;
    c.d=static_cast<std::size_t>(r.u64());c.n=static_cast<std::size_t>(r.u64());if(c.d==0||c.d>=std::numeric_limits<std::size_t>::digits||c.n!=(std::size_t{1}<<c.d))throw std::invalid_argument("invalid CRS dimensions");
    if(r.u64()!=c.n)throw std::invalid_argument("G1 count mismatch");c.gamma.reserve(c.n);for(std::size_t i=0;i<c.n;++i){auto p=r.element<G1>();require_wire_point(p,"CRS G1");c.gamma.push_back(p);}if(r.u64()!=c.n)throw std::invalid_argument("G2 count mismatch");c.lambda.reserve(c.n);for(std::size_t i=0;i<c.n;++i){auto p=r.element<G2>();require_wire_point(p,"CRS G2");c.lambda.push_back(p);}c.digest=r.digest();r.finish();if(!validate_crs(c))throw std::invalid_argument("CRS digest/validation failure");return c;
}
Statement deserialize_statement_wire(std::span<const std::uint8_t>b,const CRS&c){initialize_bn254();Reader r(b);r.expect("REXP-BF-G1-STATEMENT-WIRE-BN254-V1");r.expect("BN254");r.expect("G1");if(r.u64()!=c.d||r.u64()!=c.n)throw std::invalid_argument("statement dimensions mismatch");Statement s;s.crs_digest=r.digest();if(r.u64()!=c.n)throw std::invalid_argument("H count mismatch");s.h.reserve(c.n);for(std::size_t i=0;i<c.n;++i){auto p=r.element<G1>();require_wire_point(p,"statement H");s.h.push_back(p);}s.d1_initial=r.element<GT>();s.e0=r.element<GT>();s.f0=r.element<GT>();s.t_left0=r.element<GT>();s.t_right0=r.element<GT>();s.digest=r.digest();s.pairing_terms=r.u64();r.finish();if(!validate_statement_shape(c,s))throw std::invalid_argument("statement digest/validation failure");for(const auto*x:{&s.d1_initial,&s.e0,&s.f0,&s.t_left0,&s.t_right0})if(!mcl::bn::isValidGT(*x))throw std::invalid_argument("invalid statement GT");return s;}

Bytes serialize_proof_wire(const Proof&p,std::size_t d,std::size_t n){if(d==0||d>=std::numeric_limits<std::size_t>::digits||n!=(std::size_t{1}<<d)||p.steps.size()!=d-1)throw std::invalid_argument("proof dimensions mismatch");Bytes out;frame_text(out,"REXP-BF-G1-PROOF-WIRE-BN254-V1");frame_text(out,"BN254");frame_text(out,"G1");append(out,encode_u64_be(d));append(out,encode_u64_be(n));append(out,encode_u64_be(p.steps.size()));auto elem=[&](const auto&x){framed(out,serialize(x));};for(const auto&s:p.steps){for(const auto*x:{&s.dory_fold.d1_left,&s.dory_fold.d1_right,&s.dory_fold.d2_left,&s.dory_fold.d2_right,&s.dory_fold.w1,&s.dory_fold.w2,&s.rexp_round.e,&s.rexp_round.f,&s.rexp_round.t_left,&s.rexp_round.t_right,&s.u})elem(*x);}elem(p.phi_final);elem(p.theta_final);elem(p.r_final);return out;}
Proof deserialize_proof_wire(std::span<const std::uint8_t>b,std::size_t d,std::size_t n){initialize_bn254();Reader r(b);r.expect("REXP-BF-G1-PROOF-WIRE-BN254-V1");r.expect("BN254");r.expect("G1");if(r.u64()!=d||r.u64()!=n||d==0||d>=std::numeric_limits<std::size_t>::digits||n!=(std::size_t{1}<<d))throw std::invalid_argument("proof dimensions mismatch");if(r.u64()!=d-1)throw std::invalid_argument("proof step count mismatch");Proof p;p.steps.resize(d-1);for(auto&s:p.steps){GT* fields[]={&s.dory_fold.d1_left,&s.dory_fold.d1_right,&s.dory_fold.d2_left,&s.dory_fold.d2_right,&s.dory_fold.w1,&s.dory_fold.w2,&s.rexp_round.e,&s.rexp_round.f,&s.rexp_round.t_left,&s.rexp_round.t_right,&s.u};for(auto*x:fields){*x=r.element<GT>();if(!mcl::bn::isValidGT(*x))throw std::invalid_argument("invalid proof GT");}}p.phi_final=r.element<G1>();p.theta_final=r.element<G2>();p.r_final=r.element<G1>();r.finish();require_wire_point(p.phi_final,"proof Phi");require_wire_point(p.theta_final,"proof Theta");require_wire_point(p.r_final,"proof R");return p;}
Digest32 sha256(std::span<const std::uint8_t> bytes) {
    Digest32 out{};
    cybozu::Sha256 hash;
    if (hash.digest(out.data(), out.size(), bytes.data(), bytes.size()) != out.size())
        throw std::runtime_error("SHA-256 failed");
    return out;
}
Digest32 compute_crs_digest(const CRS& crs) {
    initialize_bn254();
    Bytes out;
    frame_text(out, "REXP-BF-G1-CRS-BN254-V1");
    frame_text(out, "BN254");
    frame_text(out, "G1-G2-GT-PAIRING");
    append(out, encode_u64_be(crs.d));
    append(out, encode_u64_be(crs.n));
    frame_text(out, "LEFT-PREFIX-GENERATOR-CHAIN-V1");
    for (const auto& p : crs.gamma) framed(out, serialize_g1(p));
    for (const auto& p : crs.lambda) framed(out, serialize_g2(p));
    return sha256(out);
}
Digest32 compute_statement_digest(const CRS& crs, const Statement& s) {
    initialize_bn254();
    Bytes out;
    frame_text(out, "REXP-BF-G1-STATEMENT-BN254-V1");
    frame_text(out, "G1");
    append_frame(out, s.crs_digest);
    append(out, encode_u64_be(crs.d));
    append(out, encode_u64_be(crs.n));
    for (const auto& p : s.h) framed(out, serialize_g1(p));
    framed(out, serialize_gt(s.d1_initial));
    framed(out, serialize_gt(s.e0));
    framed(out, serialize_gt(s.f0));
    framed(out, serialize_gt(s.t_left0));
    framed(out, serialize_gt(s.t_right0));
    return sha256(out);
}
}

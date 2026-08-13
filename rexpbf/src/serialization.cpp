#include "rexpbf/serialization.hpp"
#include "rexpbf/pairing.hpp"
#include "rexpbf/prove.hpp"
#include "rexpbf/setup.hpp"
#include "internal/crypto.hpp"
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
        std::uint64_t value = 0;
        for (int i = 0; i < 8; ++i) value = (value << 8) | bytes_[pos_++];
        return value;
    }
    std::span<const std::uint8_t> frame() {
        const auto size = u64();
        if (size > remaining()) throw std::invalid_argument("truncated frame");
        const auto n = static_cast<std::size_t>(size);
        auto out = bytes_.subspan(pos_, n);
        pos_ += n;
        return out;
    }
    void expect(std::string_view value) {
        const auto got = frame();
        if (got.size()!=value.size()||!std::equal(got.begin(),got.end(),value.begin()))
            throw std::invalid_argument("wire domain mismatch");
    }
    template<class T> T element(){return deserialize<T>(frame());}
    Digest32 digest(){const auto f=frame();if(f.size()!=32)throw std::invalid_argument("digest width mismatch");Digest32 d{};std::copy(f.begin(),f.end(),d.begin());return d;}
    std::size_t size() {
        const auto value = u64();
        if (value > std::numeric_limits<std::size_t>::max())
            throw std::invalid_argument("wire integer exceeds platform size");
        return static_cast<std::size_t>(value);
    }
    std::size_t remaining() const { return bytes_.size() - pos_; }
    void require_framed_elements(std::size_t count, std::size_t trailing_bytes = 0) const {
        constexpr std::size_t minimum_frame_bytes = 9;
        if (trailing_bytes > remaining()
            || count > (remaining() - trailing_bytes) / minimum_frame_bytes)
            throw std::invalid_argument("wire element count exceeds remaining input");
    }
    void finish(){if(pos_!=bytes_.size())throw std::invalid_argument("trailing wire bytes");}
private: std::span<const std::uint8_t> bytes_;std::size_t pos_{};
};
template<class P> void require_wire_point(const P&p,std::string_view name){
    if(p.isZero()||!p.isValid()||!p.isValidOrder())throw std::invalid_argument(std::string("invalid ")+std::string(name));
}
void require_wire_gt(const GT& value, std::string_view name) {
    if (!mcl::bn::isValidGT(value))
        throw std::invalid_argument(std::string("invalid ") + std::string(name));
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
    if(!validate_crs(crs)||!validate_statement_shape(crs,s)||!validate_statement_elements(s)
       ||!validate_statement_digest(crs,s))throw std::invalid_argument("cannot serialize invalid statement");Bytes out;
    frame_text(out,"REXP-BF-G1-STATEMENT-WIRE-BN254-V1");frame_text(out,"BN254");frame_text(out,"G1");
    append(out,encode_u64_be(crs.d));append(out,encode_u64_be(crs.n));append_frame(out,s.crs_digest);append(out,encode_u64_be(s.h.size()));
    for(const auto&p:s.h)framed(out,serialize_g1(p));for(const auto*x:{&s.d1_initial,&s.e0,&s.f0,&s.t_left0,&s.t_right0})framed(out,serialize_gt(*x));
    append_frame(out,s.digest);append(out,encode_u64_be(s.pairing_terms));return out;
}
G1 deserialize_g1(std::span<const std::uint8_t> b) { return deserialize<G1>(b); }
G2 deserialize_g2(std::span<const std::uint8_t> b) { return deserialize<G2>(b); }
GT deserialize_gt(std::span<const std::uint8_t> b) { return deserialize<GT>(b); }
CRS deserialize_crs_wire(std::span<const std::uint8_t> bytes) {
    initialize_bn254();
    Reader reader(bytes);
    reader.expect("REXP-BF-G1-CRS-WIRE-BN254-V1");
    reader.expect("BN254");
    reader.expect("G1");
    CRS crs;
    crs.d = reader.size();
    crs.n = reader.size();
    internal::validate_dimension(crs.d);
    if (crs.n != (std::size_t{1} << crs.d))
        throw std::invalid_argument("invalid CRS dimensions");
    if (reader.size() != crs.n) throw std::invalid_argument("G1 count mismatch");
    reader.require_framed_elements(crs.n, 8);
    crs.gamma.reserve(crs.n);
    for (std::size_t i = 0; i < crs.n; ++i) {
        auto point = reader.element<G1>();
        require_wire_point(point, "CRS G1");
        crs.gamma.push_back(point);
    }
    if (reader.size() != crs.n) throw std::invalid_argument("G2 count mismatch");
    reader.require_framed_elements(crs.n, 8 + 32);
    crs.lambda.reserve(crs.n);
    for (std::size_t i = 0; i < crs.n; ++i) {
        auto point = reader.element<G2>();
        require_wire_point(point, "CRS G2");
        crs.lambda.push_back(point);
    }
    crs.digest = reader.digest();
    reader.finish();
    if (!validate_crs(crs)) throw std::invalid_argument("CRS digest/validation failure");
    return crs;
}
Statement deserialize_statement_wire(std::span<const std::uint8_t> bytes,const CRS& crs) {
    initialize_bn254();
    if (!validate_crs(crs)) throw std::invalid_argument("invalid CRS");
    Reader reader(bytes);
    reader.expect("REXP-BF-G1-STATEMENT-WIRE-BN254-V1");
    reader.expect("BN254");
    reader.expect("G1");
    if (reader.size()!=crs.d || reader.size()!=crs.n)
        throw std::invalid_argument("statement dimensions mismatch");
    Statement statement;
    statement.crs_digest = reader.digest();
    if (reader.size()!=crs.n) throw std::invalid_argument("H count mismatch");
    reader.require_framed_elements(crs.n, 5 * 9 + 8 + 32 + 8);
    statement.h.reserve(crs.n);
    for (std::size_t i=0; i<crs.n; ++i) {
        auto point=reader.element<G1>();
        require_wire_point(point,"statement H");
        statement.h.push_back(point);
    }
    GT* fields[]{&statement.d1_initial, &statement.e0, &statement.f0,
                 &statement.t_left0, &statement.t_right0};
    for (GT* field : fields) {
        *field = reader.element<GT>();
        require_wire_gt(*field, "statement GT");
    }
    statement.digest=reader.digest();
    statement.pairing_terms=reader.u64();
    reader.finish();
    if (!validate_statement_shape(crs,statement)
        || !validate_statement_elements(statement)
        || !validate_statement_digest(crs,statement))
        throw std::invalid_argument("statement digest/validation failure");
    return statement;
}

Bytes serialize_proof_wire(const Proof&p,std::size_t d,std::size_t n){if(d==0||d>=std::numeric_limits<std::size_t>::digits||n!=(std::size_t{1}<<d)||p.steps.size()!=d-1)throw std::invalid_argument("proof dimensions mismatch");Bytes out;frame_text(out,"REXP-BF-G1-PROOF-WIRE-BN254-V1");frame_text(out,"BN254");frame_text(out,"G1");append(out,encode_u64_be(d));append(out,encode_u64_be(n));append(out,encode_u64_be(p.steps.size()));auto elem=[&](const auto&x){framed(out,serialize(x));};for(const auto&s:p.steps){for(const auto*x:{&s.dory_fold.d1_left,&s.dory_fold.d1_right,&s.dory_fold.d2_left,&s.dory_fold.d2_right,&s.dory_fold.w1,&s.dory_fold.w2,&s.rexp_round.e,&s.rexp_round.f,&s.rexp_round.t_left,&s.rexp_round.t_right,&s.u})elem(*x);}elem(p.phi_final);elem(p.theta_final);elem(p.r_final);return out;}
Proof deserialize_proof_wire(std::span<const std::uint8_t> bytes,std::size_t d,std::size_t n) {
    initialize_bn254();
    internal::validate_dimension(d);
    if (n != (std::size_t{1} << d)) throw std::invalid_argument("proof dimensions mismatch");
    Reader reader(bytes);
    reader.expect("REXP-BF-G1-PROOF-WIRE-BN254-V1");
    reader.expect("BN254");
    reader.expect("G1");
    if (reader.size()!=d || reader.size()!=n)
        throw std::invalid_argument("proof dimensions mismatch");
    const auto step_count = reader.size();
    if (step_count != d-1) throw std::invalid_argument("proof step count mismatch");
    constexpr std::size_t fields_per_step = 11;
    reader.require_framed_elements(step_count * fields_per_step + 3);
    Proof proof;
    proof.steps.resize(step_count);
    for (auto& step : proof.steps) {
        GT* fields[]{&step.dory_fold.d1_left,&step.dory_fold.d1_right,
            &step.dory_fold.d2_left,&step.dory_fold.d2_right,&step.dory_fold.w1,
            &step.dory_fold.w2,&step.rexp_round.e,&step.rexp_round.f,
            &step.rexp_round.t_left,&step.rexp_round.t_right,&step.u};
        for (GT* field : fields) {
            *field=reader.element<GT>();
            require_wire_gt(*field,"proof GT");
        }
    }
    proof.phi_final=reader.element<G1>();
    proof.theta_final=reader.element<G2>();
    proof.r_final=reader.element<G1>();
    reader.finish();
    require_wire_point(proof.phi_final,"proof Phi");
    require_wire_point(proof.theta_final,"proof Theta");
    require_wire_point(proof.r_final,"proof R");
    return proof;
}
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

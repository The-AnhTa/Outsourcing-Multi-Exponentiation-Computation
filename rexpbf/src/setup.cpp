#include "rexpbf/setup.hpp"
#include "rexpbf/pairing.hpp"
#include "rexpbf/serialization.hpp"
#include "internal/crypto.hpp"
#include <algorithm>
#include <chrono>
#include <stdexcept>
#include <string_view>

namespace rexpbf {
namespace {
using Bytes = std::vector<std::uint8_t>;
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b-a).count();
}
Bytes derivation_message(std::string_view domain, std::span<const std::uint8_t> seed, std::size_t i) {
    Bytes out(reinterpret_cast<const std::uint8_t*>(domain.data()),
              reinterpret_cast<const std::uint8_t*>(domain.data()) + domain.size());
    append_frame(out, seed);
    const auto index = encode_u64_be(static_cast<std::uint64_t>(i));
    out.insert(out.end(), index.begin(), index.end());
    return out;
}
G1 derive_g1(std::string_view domain, std::span<const std::uint8_t> seed, std::size_t i) {
    const auto msg = derivation_message(domain, seed, i);
    G1 p; mcl::bn::hashAndMapToG1(p, msg.data(), msg.size());
    internal::require_point(p, true, "derived G1 point"); return p;
}
G2 derive_g2(std::string_view domain, std::span<const std::uint8_t> seed, std::size_t i) {
    const auto msg = derivation_message(domain, seed, i);
    G2 p; mcl::bn::hashAndMapToG2(p, msg.data(), msg.size());
    internal::require_point(p, true, "derived G2 point"); return p;
}
bool digest_equal(const Digest32& a, const Digest32& b) { return a == b; }
Precomputation make_precomputation(const CRS& crs) {
    Precomputation p; p.crs_digest = crs.digest;
    p.x.reserve(crs.d + 1);
    p.delta1_right.reserve(crs.d); p.delta2_right.reserve(crs.d);
    for (std::size_t k = 0; k <= crs.d; ++k) {
        auto g = gamma_level(crs, k); auto l = lambda_level(crs, k);
        p.x.push_back(pairing_product(g, l)); p.pairing_terms += g.size();
    }
    for (std::size_t k = 0; k < crs.d; ++k) {
        auto g = gamma_level(crs, k); auto l = lambda_level(crs, k);
        const auto half = g.size() / 2;
        p.delta1_right.push_back(pairing_product(g.subspan(half), l.first(half)));
        p.delta2_right.push_back(pairing_product(g.first(half), l.subspan(half)));
        p.pairing_terms += static_cast<std::uint64_t>(2 * half);
    }
    if (p.pairing_terms != 4 * static_cast<std::uint64_t>(crs.n) - 3)
        throw std::logic_error("precomputation pairing count invariant failed");
    return p;
}
Statement make_statement(const CRS& crs, std::span<const std::uint8_t> seed, SetupBreakdown* profile) {
    Statement s; s.crs_digest = crs.digest; s.h.reserve(crs.n);
    auto phase = Clock::now();
    for (std::size_t i = 0; i < crs.n; ++i)
        s.h.push_back(derive_g1("REXP-BF-G1-PUBLIC-H-V1", seed, i));
    if (profile) profile->h_generation_ms = elapsed(phase, Clock::now());
    phase = Clock::now();
    const auto lambda = std::span<const G2>(crs.lambda);
    const auto h = std::span<const G1>(s.h);
    const auto half = crs.n / 2;
    s.d1_initial = pairing_product(h, lambda); s.pairing_terms += crs.n;
    s.e0 = pairing_product(h.subspan(half), lambda.first(half)); s.pairing_terms += half;
    s.f0 = pairing_product(h.first(half), lambda.subspan(half)); s.pairing_terms += half;
    s.t_left0 = pairing_product(h.first(half), lambda.first(half)); s.pairing_terms += half;
    s.t_right0 = pairing_product(h.subspan(half), lambda.first(half)); s.pairing_terms += half;
    if (s.pairing_terms != 3 * static_cast<std::uint64_t>(crs.n))
        throw std::logic_error("statement pairing count invariant failed");
    if (profile) profile->statement_gt_ms = elapsed(phase, Clock::now());
    phase = Clock::now();
    s.digest = compute_statement_digest(crs, s);
    if (profile) profile->statement_digest_ms = elapsed(phase, Clock::now());
    return s;
}
SetupResult setup_impl(const SetupConfig& config, SetupBreakdown* profile) {
    auto total = Clock::now();
    initialize_bn254(); internal::validate_dimension(config.d);
    SetupResult out;
    out.crs.d = config.d; out.crs.n = std::size_t{1} << config.d;
    out.crs.gamma.reserve(out.crs.n); out.crs.lambda.reserve(out.crs.n);
    auto phase = Clock::now();
    for (std::size_t i = 0; i < out.crs.n; ++i) {
        out.crs.gamma.push_back(derive_g1("REXP-BF-G1-CRS-GAMMA-V1", config.crs_seed, i));
        out.crs.lambda.push_back(derive_g2("REXP-BF-G1-CRS-LAMBDA-V1", config.crs_seed, i));
    }
    if (profile) profile->crs_generation_ms = elapsed(phase, Clock::now());
    phase = Clock::now(); out.crs.digest = compute_crs_digest(out.crs);
    if (profile) profile->crs_digest_ms = elapsed(phase, Clock::now());
    phase = Clock::now(); out.precomputation = make_precomputation(out.crs);
    if (profile) profile->precomputation_ms = elapsed(phase, Clock::now());
    out.statement = make_statement(out.crs, config.instance_seed, profile);
    out.prover_input.h = out.statement.h;
    if (!validate_crs(out.crs)
        || !validate_precomputation_shape(out.crs, out.precomputation)
        || !validate_precomputation_elements(out.precomputation)
        || !validate_statement_shape(out.crs, out.statement)
        || !validate_statement_elements(out.statement)
        || !validate_statement_digest(out.crs, out.statement))
        throw std::logic_error("Setup validation failed");
    if (profile) profile->total_ms = elapsed(total, Clock::now());
    return out;
}
}

std::span<const G1> gamma_level(const CRS& crs, std::size_t level) {
    if (level > crs.d) throw std::out_of_range("gamma level exceeds d");
    return {crs.gamma.data(), crs.n >> level};
}
std::span<const G2> lambda_level(const CRS& crs, std::size_t level) {
    if (level > crs.d) throw std::out_of_range("lambda level exceeds d");
    return {crs.lambda.data(), crs.n >> level};
}
SetupResult setup(const SetupConfig& config) {
    return setup_impl(config, nullptr);
}
SetupResult setup_with_breakdown(const SetupConfig& config, SetupBreakdown& breakdown) {
    breakdown = {}; return setup_impl(config, &breakdown);
}
bool validate_crs(const CRS& crs) {
    return validate_crs_shape(crs) && validate_crs_points(crs)
        && validate_crs_digest(crs);
}
bool validate_crs_shape(const CRS& crs) {
    try {
        internal::validate_dimension(crs.d);
        return crs.n == (std::size_t{1} << crs.d)
            && crs.gamma.size() == crs.n && crs.lambda.size() == crs.n;
    } catch (...) { return false; }
}
bool validate_crs_points(const CRS& crs) {
    try {
        for (const auto& p : crs.gamma)
            internal::require_point(p, true, "CRS G1 point");
        for (const auto& p : crs.lambda)
            internal::require_point(p, true, "CRS G2 point");
        return true;
    } catch (...) { return false; }
}
bool validate_crs_digest(const CRS& crs) {
    try { return digest_equal(crs.digest, compute_crs_digest(crs)); }
    catch (...) { return false; }
}
bool validate_precomputation_shape(const CRS& c, const Precomputation& p) {
    return validate_crs_shape(c) && p.crs_digest == c.digest && p.x.size() == c.d + 1
        && p.delta1_right.size() == c.d && p.delta2_right.size() == c.d
        && p.pairing_terms == 4 * static_cast<std::uint64_t>(c.n) - 3;
}
bool validate_precomputation_elements(const Precomputation& p) {
    const auto valid = [](const GT& value) { return mcl::bn::isValidGT(value); };
    return std::all_of(p.x.begin(), p.x.end(), valid)
        && std::all_of(p.delta1_right.begin(), p.delta1_right.end(), valid)
        && std::all_of(p.delta2_right.begin(), p.delta2_right.end(), valid);
}
bool validate_statement_shape(const CRS& c, const Statement& s) {
    return validate_crs_shape(c) && s.crs_digest == c.digest && s.h.size() == c.n
        && s.pairing_terms == 3 * static_cast<std::uint64_t>(c.n);
}
bool validate_statement_elements(const Statement& s) {
    try {
        for (const auto& point : s.h)
            internal::require_point(point, true, "statement G1 point");
        return mcl::bn::isValidGT(s.d1_initial) && mcl::bn::isValidGT(s.e0)
            && mcl::bn::isValidGT(s.f0) && mcl::bn::isValidGT(s.t_left0)
            && mcl::bn::isValidGT(s.t_right0);
    } catch (...) { return false; }
}
bool validate_statement_digest(const CRS& c, const Statement& s) {
    try { return s.digest == compute_statement_digest(c, s); }
    catch (...) { return false; }
}
bool audit_precomputation(const CRS& c, const Precomputation& p) {
    if (!validate_crs(c) || !validate_precomputation_shape(c, p)
        || !validate_precomputation_elements(p)) return false;
    const auto expected = make_precomputation(c);
    return expected.x == p.x && expected.delta1_right == p.delta1_right
        && expected.delta2_right == p.delta2_right;
}
bool audit_statement(const CRS& c, const Statement& s) {
    if (!validate_crs(c) || !validate_statement_shape(c, s)
        || !validate_statement_elements(s) || !validate_statement_digest(c, s)) return false;
    const auto l = std::span<const G2>(c.lambda); const auto h = std::span<const G1>(s.h);
    const auto half = c.n / 2;
    return s.d1_initial == pairing_product(h, l)
        && s.e0 == pairing_product(h.subspan(half), l.first(half))
        && s.f0 == pairing_product(h.first(half), l.subspan(half))
        && s.t_left0 == pairing_product(h.first(half), l.first(half))
        && s.t_right0 == pairing_product(h.subspan(half), l.first(half));
}
}

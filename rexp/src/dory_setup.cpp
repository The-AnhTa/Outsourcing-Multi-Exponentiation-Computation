#include "rexp/dory_setup.hpp"
#include "internal/crypto.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace rexp {
namespace {

using internal::Bytes;
using internal::append;
using internal::append_fixed_element;
using internal::append_frame;
using internal::append_framed_element;
using internal::append_u64_be;
using internal::pairing_product;
using internal::serialize_element;
using internal::sha256;

constexpr std::string_view kCrsDomain = "DORY-CRS-BN254-V1";
constexpr std::string_view kCurve = "BN254";
constexpr std::string_view kChain = "LEFT-PREFIX-GENERATOR-CHAIN-V1";
constexpr std::string_view kG1Domain = "REXP-CRS-G1-V1";
constexpr std::string_view kG2Domain = "REXP-CRS-G2-V1";
constexpr std::string_view kCrsWire = "DORY-CRS-WIRE-BN254-V1";
constexpr std::string_view kPrecompWire = "DORY-PRECOMP-WIRE-BN254-V1";
constexpr std::string_view kStatementWire = "DORY-STATEMENT-WIRE-BN254-V1";

G1 hash_g1(std::string_view domain, std::string_view seed, std::size_t index) {
    Bytes message;
    append(message, domain);
    append_frame(message, seed);
    append_u64_be(message, static_cast<std::uint64_t>(index));
    G1 point;
    mcl::bn::hashAndMapToG1(point, message.data(), message.size());
    return point;
}

G2 hash_g2(std::string_view domain, std::string_view seed, std::size_t index) {
    Bytes message;
    append(message, domain);
    append_frame(message, seed);
    append_u64_be(message, static_cast<std::uint64_t>(index));
    G2 point;
    mcl::bn::hashAndMapToG2(point, message.data(), message.size());
    return point;
}

template<class Point>
void require_point(const Point& point, const char* name) {
    if (point.isZero() || !point.isValid() || !point.isValidOrder()) {
        throw std::runtime_error(std::string("invalid generated ") + name);
    }
}

Fr default_scalar() {
    Fr value;
    do {
        value.setByCSPRNG();
    } while (value.isZero());
    return value;
}

Fr sample_nonzero(const NonzeroScalarSampler& sampler) {
    Fr value = sampler ? sampler() : default_scalar();
    if (value.isZero()) throw std::invalid_argument("RNG returned a zero scalar");
    return value;
}

std::pair<DoryStatement, DoryWitness> generate_instance(
    const DoryCRS& crs,
    std::string_view crs_seed,
    const NonzeroScalarSampler& sampler) {
    DoryWitness witness;
    witness.Phi.reserve(crs.n);
    witness.Theta.reserve(crs.n);
    G1 g1_base = hash_g1("DORY-WITNESS-G1-BASE-V1", crs_seed, 0);
    G2 g2_base = hash_g2("DORY-WITNESS-G2-BASE-V1", crs_seed, 0);
    require_point(g1_base, "witness G1 base");
    require_point(g2_base, "witness G2 base");
    for (std::size_t i = 0; i < crs.n; ++i) {
        G1 phi;
        G2 theta;
        G1::mul(phi, g1_base, sample_nonzero(sampler));
        G2::mul(theta, g2_base, sample_nonzero(sampler));
        witness.Phi.push_back(phi);
        witness.Theta.push_back(theta);
    }

    DoryStatement statement;
    statement.D0 = pairing_product(
        witness.Phi, 0, witness.Theta, 0, crs.n);
    statement.D1 = pairing_product(
        witness.Phi, 0, crs.Lambda, 0, crs.n);
    statement.D2 = pairing_product(
        crs.Gamma, 0, witness.Theta, 0, crs.n);
    return {std::move(statement), std::move(witness)};
}

bool same_digest(const Digest& a, const Digest& b) {
    return std::equal(a.begin(), a.end(), b.begin());
}

void validate_d(std::size_t d) {
    if (d >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument("d is too large for size_t");
    }

    const std::size_t n = std::size_t{1} << d;
    if (n > (std::numeric_limits<std::size_t>::max() - 3U) / 4U) {
        throw std::invalid_argument("d makes the Setup pairing count overflow");
    }
}

}

void initialize() {
    static std::once_flag once;
    std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}

Digest ComputeDoryCRSDigest(const DoryCRS& crs) {
    initialize();
    if (crs.n == 0 || crs.Gamma.size() != crs.n || crs.Lambda.size() != crs.n) {
        throw std::invalid_argument("cannot digest malformed CRS");
    }
    Bytes input;
    append_frame(input, kCrsDomain);
    append_frame(input, kCurve);
    append_u64_be(input, static_cast<std::uint64_t>(crs.d));
    append_u64_be(input, static_cast<std::uint64_t>(crs.n));
    append_frame(input, kChain);
    for (const G1& point : crs.Gamma) append_framed_element(input, point);
    for (const G2& point : crs.Lambda) append_framed_element(input, point);
    return sha256(input);
}

SetupResult SetupDory(
    std::size_t d,
    std::string_view crs_seed,
    NonzeroScalarSampler sampler) {
    BatchSetupResult batch = SetupDoryBatch(d, 1, crs_seed, std::move(sampler));
    SetupResult out;
    out.crs = std::move(batch.crs);
    out.precomp = std::move(batch.precomp);
    out.statement = std::move(batch.statements.front());
    out.witness = std::move(batch.witnesses.front());
    return out;
}

DoryCRS GenerateDoryCRS(std::size_t d, std::string_view crs_seed) {
    initialize();
    validate_d(d);
    const std::size_t n = std::size_t{1} << d;
    DoryCRS out;
    out.d = d;
    out.n = n;
    out.Gamma.reserve(n);
    out.Lambda.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        G1 gamma = hash_g1(kG1Domain, crs_seed, i);
        G2 lambda = hash_g2(kG2Domain, crs_seed, i);
        require_point(gamma, "Gamma");
        require_point(lambda, "Lambda");
        out.Gamma.push_back(gamma);
        out.Lambda.push_back(lambda);
    }
    out.digest = ComputeDoryCRSDigest(out);
    return out;
}

DoryPrecomputation PrepareDoryPrecomputation(const DoryCRS& crs) {
    if (!ValidateCRS(crs)) {
        throw std::invalid_argument("cannot prepare malformed Dory CRS");
    }
    DoryPrecomputation out;
    out.X.reserve(crs.d + 1);
    for (std::size_t k = 0; k <= crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        out.X.push_back(pairing_product(crs.Gamma, 0, crs.Lambda, 0, m));
        out.pairing_product_terms += m;
    }
    out.Delta1R.reserve(crs.d);
    out.Delta2R.reserve(crs.d);
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        const std::size_t h = m / 2;
        out.Delta1R.push_back(
            pairing_product(crs.Gamma, h, crs.Lambda, 0, h));
        out.Delta2R.push_back(
            pairing_product(crs.Gamma, 0, crs.Lambda, h, h));
        out.pairing_product_terms += 2 * h;
    }
    out.crs_digest = crs.digest;
    return out;
}

BatchSetupResult SetupDoryBatch(
    std::size_t d,
    std::size_t batch_size,
    std::string_view crs_seed,
    NonzeroScalarSampler sampler) {
    if (batch_size == 0) {
        throw std::invalid_argument("batch size B must be at least one");
    }

    BatchSetupResult out;
    out.crs = GenerateDoryCRS(d, crs_seed);
    out.precomp = PrepareDoryPrecomputation(out.crs);

    out.statements.reserve(batch_size);
    out.witnesses.reserve(batch_size);
    for (std::size_t instance = 0; instance < batch_size; ++instance) {
        auto generated = generate_instance(out.crs, crs_seed, sampler);
        out.statements.push_back(std::move(generated.first));
        out.witnesses.push_back(std::move(generated.second));
    }
    return out;
}

bool ValidateCRS(const DoryCRS& crs) {
    try {
        validate_d(crs.d);
        if (crs.n != (std::size_t{1} << crs.d)
            || crs.Gamma.size() != crs.n || crs.Lambda.size() != crs.n) return false;
        for (const G1& point : crs.Gamma) require_point(point, "Gamma");
        for (const G2& point : crs.Lambda) require_point(point, "Lambda");
        return same_digest(crs.digest, ComputeDoryCRSDigest(crs));
    } catch (...) {
        return false;
    }
}

bool ValidatePrecomputation(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp) {
    if (!ValidateCRS(crs) || precomp.X.size() != crs.d + 1
        || precomp.Delta1R.size() != crs.d
        || precomp.Delta2R.size() != crs.d
        || !same_digest(precomp.crs_digest, crs.digest)
        || precomp.pairing_product_terms != 4 * crs.n - 3) return false;
    for (std::size_t k = 0; k <= crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        if (precomp.X[k] != pairing_product(crs.Gamma, 0, crs.Lambda, 0, m)) {
            return false;
        }
    }
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        const std::size_t h = m / 2;
        if (precomp.Delta1R[k]
                != pairing_product(crs.Gamma, h, crs.Lambda, 0, h)
            || precomp.Delta2R[k]
                != pairing_product(crs.Gamma, 0, crs.Lambda, h, h)) return false;
    }
    return true;
}

std::vector<std::uint8_t> SerializeCRS(const DoryCRS& crs) {
    if (!ValidateCRS(crs)) throw std::invalid_argument("cannot serialize invalid CRS");
    Bytes out;
    append_frame(out, kCrsWire);
    append_frame(out, kCurve);
    append_u64_be(out, crs.d);
    append_u64_be(out, crs.n);
    append_frame(out, kChain);
    append_u64_be(out, crs.Gamma.size());



    const std::size_t g1_size = serialize_element(crs.Gamma.front()).size();
    const std::size_t g2_size = serialize_element(crs.Lambda.front()).size();
    for (const G1& point : crs.Gamma) {
        append_fixed_element(out, point, g1_size, "G1");
    }
    append_u64_be(out, crs.Lambda.size());
    for (const G2& point : crs.Lambda) {
        append_fixed_element(out, point, g2_size, "G2");
    }
    append_frame(out, crs.digest.data(), crs.digest.size());
    return out;
}

std::vector<std::uint8_t> SerializePrecomputation(
    const DoryPrecomputation& precomp) {
    Bytes out;
    append_frame(out, kPrecompWire);
    append_frame(out, precomp.crs_digest.data(), precomp.crs_digest.size());
    append_u64_be(out, precomp.X.size());
    for (const GT& value : precomp.X) append_framed_element(out, value);
    append_u64_be(out, precomp.Delta1R.size());
    for (const GT& value : precomp.Delta1R) append_framed_element(out, value);
    append_u64_be(out, precomp.Delta2R.size());
    for (const GT& value : precomp.Delta2R) append_framed_element(out, value);
    append_u64_be(out, precomp.pairing_product_terms);
    return out;
}

std::vector<std::uint8_t> SerializeStatement(const DoryStatement& statement) {
    Bytes out;
    append_frame(out, kStatementWire);
    append_frame(out, kCurve);
    append_framed_element(out, statement.D0);
    append_framed_element(out, statement.D1);
    append_framed_element(out, statement.D2);
    return out;
}

SetupSizes MeasureSetupSizes(const SetupResult& setup) {
    return {
        SerializeCRS(setup.crs).size(),
        SerializePrecomputation(setup.precomp).size(),
        SerializeStatement(setup.statement).size()};
}

SetupSizes MeasureSetupSizes(const BatchSetupResult& setup) {
    std::size_t statements = 0;
    for (const DoryStatement& statement : setup.statements) {
        statements += SerializeStatement(statement).size();
    }
    return {
        SerializeCRS(setup.crs).size(),
        SerializePrecomputation(setup.precomp).size(),
        statements};
}

std::string DigestHex(const Digest& digest) {
    std::ostringstream out;
    out << std::hex << std::setfill('0');
    for (std::uint8_t byte : digest) out << std::setw(2) << unsigned(byte);
    return out.str();
}

}

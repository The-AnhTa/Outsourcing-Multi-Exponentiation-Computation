#include "rexp/serialization.hpp"

#include "internal/crypto.hpp"
#include "internal/rexp_helpers.hpp"
#include "internal/rexp_transcript.hpp"

#include <chrono>
#include <limits>
#include <stdexcept>
#include <utility>

namespace rexp {
namespace {
using Clock = std::chrono::steady_clock;
}

std::vector<std::uint8_t> SerializeRexpCRS(const RawRexpCRS& crs) {
    internal::check_dimension(crs.d, crs.n, "Rexp");
    internal::Bytes bytes;
    internal::append_u64(bytes, crs.d);
    internal::append_u64(bytes, crs.n);
    for (const auto& point : crs.Gamma) internal::frame_element(bytes, point);
    for (const auto& point : crs.Lambda) internal::frame_element(bytes, point);
    const Digest digest = internal::rexp_crs_digest(crs);
    internal::append(bytes, digest.data(), digest.size());
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpPrecomputation(
    const PreparedPublicParameters& params) {
    internal::Bytes bytes;
    for (const auto& value : params.X()) internal::frame_element(bytes, value);
    for (const auto& value : params.Delta1R()) internal::frame_element(bytes, value);
    for (const auto& value : params.Delta2R()) internal::frame_element(bytes, value);
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpStatement(
    const RawRexpStatement& statement) {
    internal::Bytes bytes;
    for (const auto& point : statement.H) internal::frame_element(bytes, point);
    return bytes;
}

std::vector<std::uint8_t> SerializePreparedStatement(
    const PreparedStatement& statement) {
    internal::Bytes bytes;
    for (const auto& point : statement.H()) internal::frame_element(bytes, point);
    internal::frame_element(bytes, statement.D1Initial());
    internal::frame_element(bytes, statement.E0());
    internal::frame_element(bytes, statement.F0());
    internal::frame_element(bytes, statement.TL0());
    internal::frame_element(bytes, statement.TR0());
    internal::append(bytes, statement.digest().data(), statement.digest().size());
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpProof(
    const RexpProof& proof, std::size_t d) {
    if (d >= std::numeric_limits<std::size_t>::digits
        || proof.doryProofs.size() != d
        || proof.dynamicRoundMessages.size() != (d ? d - 1 : 0)) {
        throw std::invalid_argument("wrong proof shape");
    }
    internal::Bytes bytes;
    for (std::size_t round = 0; round < d; ++round) {
        if (proof.doryProofs[round].rounds.size() != d - round - 1) {
            throw std::invalid_argument("wrong embedded Dory depth");
        }
        if (round) {
            const auto& message = proof.dynamicRoundMessages[round - 1];
            internal::frame_element(bytes, message.E);
            internal::frame_element(bytes, message.F);
            internal::frame_element(bytes, message.TL);
            internal::frame_element(bytes, message.TR);
        }
        for (const auto& message : proof.doryProofs[round].rounds) {
            internal::frame_element(bytes, message.D1L);
            internal::frame_element(bytes, message.D1R);
            internal::frame_element(bytes, message.D2L);
            internal::frame_element(bytes, message.D2R);
            internal::frame_element(bytes, message.W1);
            internal::frame_element(bytes, message.W2);
        }
        internal::frame_element(bytes, proof.doryProofs[round].PhiFinal);
        internal::frame_element(bytes, proof.doryProofs[round].ThetaFinal);
    }
    internal::frame_element(bytes, proof.R);
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpCRSWire(const RawRexpCRS& crs) {
    internal::check_dimension(crs.d, crs.n, "Rexp");
    internal::Bytes bytes;
    internal::frame(bytes, "REXP-G1-CRS-WIRE-BN254-V1");
    internal::append_u64(bytes, crs.d);
    internal::append_u64(bytes, crs.n);
    internal::append_u64(bytes, crs.Gamma.size());
    for (const auto& point : crs.Gamma) internal::frame_element(bytes, point);
    internal::append_u64(bytes, crs.Lambda.size());
    for (const auto& point : crs.Lambda) internal::frame_element(bytes, point);
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpStatementWire(
    const RawRexpStatement& statement, std::size_t n) {
    if (statement.H.size() != n) {
        throw std::invalid_argument("statement wire length mismatch");
    }
    internal::Bytes bytes;
    internal::frame(bytes, "REXP-G1-RAW-STATEMENT-WIRE-BN254-V1");
    internal::append_u64(bytes, n);
    internal::append_u64(bytes, statement.H.size());
    for (const auto& point : statement.H) internal::frame_element(bytes, point);
    return bytes;
}

std::vector<std::uint8_t> SerializeRexpProofWire(
    const RexpProof& proof, std::size_t d) {
    (void)SerializeRexpProof(proof, d);
    internal::Bytes bytes;
    internal::frame(bytes, "REXP-G1-PROOF-WIRE-BN254-V1");
    internal::append_u64(bytes, d);
    internal::append_u64(bytes, proof.dynamicRoundMessages.size());
    internal::append_u64(bytes, proof.doryProofs.size());
    for (std::size_t round = 0; round < d; ++round) {
        if (round) {
            const auto& message = proof.dynamicRoundMessages[round - 1];
            internal::frame_element(bytes, message.E);
            internal::frame_element(bytes, message.F);
            internal::frame_element(bytes, message.TL);
            internal::frame_element(bytes, message.TR);
        }
        internal::append_u64(bytes, proof.doryProofs[round].rounds.size());
        for (const auto& message : proof.doryProofs[round].rounds) {
            internal::frame_element(bytes, message.D1L);
            internal::frame_element(bytes, message.D1R);
            internal::frame_element(bytes, message.D2L);
            internal::frame_element(bytes, message.D2R);
            internal::frame_element(bytes, message.W1);
            internal::frame_element(bytes, message.W2);
        }
        internal::frame_element(bytes, proof.doryProofs[round].PhiFinal);
        internal::frame_element(bytes, proof.doryProofs[round].ThetaFinal);
    }
    internal::frame_element(bytes, proof.R);
    return bytes;
}

RawRexpCRS DeserializeRexpCRSWire(const internal::Bytes& bytes) {
    internal::Reader reader(bytes);
    reader.expect("REXP-G1-CRS-WIRE-BN254-V1");
    RawRexpCRS crs;
    crs.d = reader.readSize();
    crs.n = reader.readSize();
    internal::check_dimension(crs.d, crs.n, "Rexp");
    if (reader.readU64() != crs.n) throw std::invalid_argument("G1 count mismatch");
    reader.requireFramedElementsFit(crs.n);
    crs.Gamma.reserve(crs.n);
    for (std::size_t i = 0; i < crs.n; ++i) crs.Gamma.push_back(reader.element<G1>());
    if (reader.readU64() != crs.n) throw std::invalid_argument("G2 count mismatch");
    reader.requireFramedElementsFit(crs.n);
    crs.Lambda.reserve(crs.n);
    for (std::size_t i = 0; i < crs.n; ++i) crs.Lambda.push_back(reader.element<G2>());
    reader.finish();
    (void)PreparePublicParameters(crs);
    return crs;
}

RawRexpStatement DeserializeRexpStatementWire(
    const internal::Bytes& bytes, std::size_t expected_n) {
    internal::Reader reader(bytes);
    reader.expect("REXP-G1-RAW-STATEMENT-WIRE-BN254-V1");
    if (reader.readU64() != expected_n || reader.readU64() != expected_n) {
        throw std::invalid_argument("statement count mismatch");
    }
    reader.requireFramedElementsFit(expected_n);
    RawRexpStatement statement;
    statement.H.reserve(expected_n);
    for (std::size_t i = 0; i < expected_n; ++i) {
        statement.H.push_back(reader.element<G1>());
    }
    reader.finish();
    for (const auto& point : statement.H) {
        internal::require_point(point, false, "statement G1 point");
    }
    return statement;
}

RexpProof DeserializeRexpProofWire(
    const internal::Bytes& bytes,
    std::size_t expected_d,
    RexpProofValidationMetrics* metrics) {
    if (expected_d >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument("proof d too large");
    }
    const auto decode_start = Clock::now();
    internal::Reader reader(bytes);
    reader.expect("REXP-G1-PROOF-WIRE-BN254-V1");
    if (reader.readU64() != expected_d) throw std::invalid_argument("proof d mismatch");
    const std::size_t dynamic = expected_d ? expected_d - 1 : 0;
    if (reader.readU64() != dynamic || reader.readU64() != expected_d) {
        throw std::invalid_argument("proof count mismatch");
    }
    RexpProof proof;
    proof.dynamicRoundMessages.reserve(dynamic);
    proof.doryProofs.reserve(expected_d);
    for (std::size_t round = 0; round < expected_d; ++round) {
        if (round) proof.dynamicRoundMessages.push_back({
            reader.element<GT>(), reader.element<GT>(),
            reader.element<GT>(), reader.element<GT>()});
        const std::size_t depth = expected_d - round - 1;
        if (reader.readU64() != depth) throw std::invalid_argument("Dory depth mismatch");
        DoryProof dory;
        dory.rounds.reserve(depth);
        for (std::size_t i = 0; i < depth; ++i) dory.rounds.push_back({
            reader.element<GT>(), reader.element<GT>(), reader.element<GT>(),
            reader.element<GT>(), reader.element<GT>(), reader.element<GT>()});
        dory.PhiFinal = reader.element<G1>();
        dory.ThetaFinal = reader.element<G2>();
        proof.doryProofs.push_back(std::move(dory));
    }
    proof.R = reader.element<G1>();
    reader.finish();
    if (!proof.R.isValid() || !proof.R.isValidOrder()) {
        throw std::invalid_argument("invalid proof R G1 point");
    }
    for (const auto& dory : proof.doryProofs) {
        if (!dory.PhiFinal.isValid() || !dory.PhiFinal.isValidOrder()
            || !dory.ThetaFinal.isValid() || !dory.ThetaFinal.isValidOrder()) {
            throw std::invalid_argument("invalid embedded Dory final point");
        }
    }
    RexpProofValidationMetrics local;
    internal::validate_rexp_proof_gt(proof, &local);
    if (metrics) {
        *metrics = local;
        metrics->canonical_decode_ms =
            std::chrono::duration<double, std::milli>(Clock::now() - decode_start)
                .count() - local.gt_subgroup_validation_ms;
    }
    return proof;
}

ValidatedRexpProof DeserializeValidatedRexpProofWire(
    const internal::Bytes& bytes,
    std::size_t d,
    RexpProofValidationMetrics* metrics) {
    RexpProof proof = DeserializeRexpProofWire(bytes, d, metrics);
    // Deserialization has already performed shape, point, canonical encoding,
    // and GT subgroup checks. Move the validated representation into its
    // wrapper so large proofs are not copied.
    ValidatedRexpProof result;
    result.proof_ = std::move(proof);
    result.d_ = d;
    return result;
}

} // namespace rexp

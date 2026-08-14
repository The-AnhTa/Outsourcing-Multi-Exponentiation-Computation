#include "rexp/rexp.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Bytes = std::vector<std::uint8_t>;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Function>
void check_throws(Function&& function, const std::string& message) {
    try {
        function();
    } catch (const std::exception&) {
        return;
    }
    throw std::runtime_error(message);
}

std::string hex(const std::uint8_t* data, std::size_t size) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string out;
    out.reserve(2 * size);
    for (std::size_t i = 0; i < size; ++i) {
        out.push_back(digits[data[i] >> 4]);
        out.push_back(digits[data[i] & 0x0f]);
    }
    return out;
}

template<class T>
std::string encoded_hex(const T& value) {
    Bytes bytes(2048);
    const std::size_t written = value.serialize(bytes.data(), bytes.size());
    check(written != 0, "element serialization failed in test");
    return hex(bytes.data(), written);
}

std::string sha256_hex(const Bytes& bytes) {
    rexp::Digest digest{};
    const auto written = mcl::fp::sha256(
        digest.data(), static_cast<std::uint32_t>(digest.size()),
        bytes.data(), static_cast<std::uint32_t>(bytes.size()));
    check(written == digest.size(), "SHA-256 failed in test");
    return hex(digest.data(), digest.size());
}

void append_u64(Bytes& bytes, std::uint64_t value) {
    for (int shift = 56; shift >= 0; shift -= 8) {
        bytes.push_back(static_cast<std::uint8_t>(value >> shift));
    }
}

void append_frame(Bytes& bytes, const std::string& text) {
    append_u64(bytes, text.size());
    bytes.insert(bytes.end(), text.begin(), text.end());
}

rexp::NonzeroScalarSampler sequential_sampler() {
    auto counter = std::make_shared<std::uint64_t>(0);
    return [counter] {
        rexp::Fr value;
        value = ++*counter;
        return value;
    };
}

void test_challenge_vector() {
    rexp::Digest digest{};
    for (std::size_t i = 0; i < digest.size(); ++i) {
        digest[i] = static_cast<std::uint8_t>(i);
    }
    const rexp::Fr challenge =
        rexp::ChallengeNonzeroFr(digest, "DORY-BETA-V1", 2);
    const std::string actual = encoded_hex(challenge);
    const std::string expected =
        "29f2dfab6428b5846a94660058c573665b006b25984c7d58288e52f9342cc70f";
    check(actual == expected, "challenge vector mismatch: " + actual);
}

void test_dory_setup_and_proofs() {
    const auto first = rexp::SetupDory(2, "dory-test-seed", sequential_sampler());
    const auto second = rexp::SetupDory(2, "dory-test-seed", sequential_sampler());
    check(first.crs.digest == second.crs.digest, "CRS generation is not deterministic");
    check(first.crs.Gamma == second.crs.Gamma, "G1 CRS generation differs");
    check(first.crs.Lambda == second.crs.Lambda, "G2 CRS generation differs");
    check(rexp::ValidateCRS(first.crs), "generated CRS is invalid");
    check(rexp::ValidatePrecomputation(first.crs, first.precomp),
          "generated precomputation is invalid");

    auto malformed_crs = first.crs;
    malformed_crs.Gamma.pop_back();
    check(!rexp::ValidateCRS(malformed_crs), "malformed CRS was accepted");

    const auto proof = rexp::Prove(first.crs, first.statement, first.witness);
    rexp::VerifyMetrics metrics;
    check(rexp::Verify(first.crs, first.precomp, first.statement, proof, &metrics),
          "ordinary Dory proof was rejected");
    check(metrics.terminal_pairings == 1, "ordinary verification metric mismatch");

    const auto prepared = rexp::PrepareVerifier(first.crs, first.precomp);
    check(prepared.has_value(), "Dory verifier preparation failed");
    check(rexp::VerifyPrepared(*prepared, first.statement, proof),
          "prepared Dory verification failed");
    rexp::VerifyAuditTrace trace;
    check(rexp::VerifyPreparedAudit(*prepared, first.statement, proof, metrics, trace),
          "audited Dory verification failed");
    check(trace.evaluated_lhs == trace.terminal_pairing_rhs,
          "Dory audit sides differ");

    auto malformed_proof = proof;
    malformed_proof.rounds.pop_back();
    check(!rexp::Verify(first.crs, first.precomp, first.statement, malformed_proof),
          "wrong Dory proof shape was accepted");

    rexp::Digest transcript{};
    transcript[0] = 0x42;
    rexp::Digest prover_end{}, verifier_end{}, reference_end{};
    const auto embedded = rexp::ProveEmbedded(
        first.crs, first.witness, transcript, &prover_end);
    check(rexp::VerifyEmbeddedDeferred(
              first.crs, first.precomp, first.statement, embedded,
              transcript, &verifier_end),
          "embedded deferred Dory verification failed");
    check(rexp::VerifyEmbeddedReference(
              first.crs, first.precomp, first.statement, embedded,
              transcript, &reference_end),
          "embedded reference Dory verification failed");
    check(prover_end == verifier_end && prover_end == reference_end,
          "embedded transcript end differs");

    const std::string proof_hash = sha256_hex(rexp::SerializeProof(proof, 2));
    const std::string expected_hash =
        "46bab83a6720c8353004194239f2313a1efe342ef97c23af83e9b18de12ccceb";
    check(proof_hash == expected_hash, "Dory proof wire vector mismatch: " + proof_hash);
}

void test_dory_batch() {
    const auto setup =
        rexp::SetupDoryBatch(2, 3, "dory-batch-test-seed", sequential_sampler());
    const auto proof = rexp::ProveBatch(
        setup.crs, setup.statements, setup.witnesses);
    rexp::VerifyMetrics metrics;
    check(rexp::VerifyBatch(
              setup.crs, setup.precomp, setup.statements, proof, &metrics),
          "batch Dory proof was rejected");
    const auto prepared = rexp::PrepareVerifier(setup.crs, setup.precomp);
    check(prepared.has_value(), "batch verifier preparation failed");
    check(rexp::VerifyBatchPrepared(*prepared, setup.statements, proof),
          "prepared batch verification failed");
    rexp::VerifyAuditTrace trace;
    check(rexp::VerifyBatchPreparedAudit(
              *prepared, setup.statements, proof, metrics, trace),
          "audited batch verification failed");
    check(rexp::DeriveBatchGammas(setup.crs, setup.statements, proof).size() == 2,
          "batch challenge count mismatch");
    check(rexp::DeriveBatchChallenges(setup.crs, setup.statements, proof)
              .dory.beta.size() == 2,
          "batch Dory challenge count mismatch");

    auto malformed = proof;
    malformed.batchCrossTerms.pop_back();
    check(!rexp::VerifyBatch(
              setup.crs, setup.precomp, setup.statements, malformed),
          "malformed batch proof was accepted");
}

void test_rexp_dimension(std::size_t d) {
    const auto setup = rexp::Setup(d, "rexp-test-seed");
    const auto proof = rexp::Prove(setup.params, setup.statement, setup.proverInput);
    static constexpr std::size_t expected_proof_bytes[] = {40, 152, 4184};
    static constexpr std::size_t expected_crs_bytes[] = {160, 272, 496};
    if (d < 3) {
        check(rexp::SerializeRexpProof(proof, d).size()
                  == expected_proof_bytes[d],
              "Rexp proof size regression");
        check(rexp::SerializeRexpCRS(setup.rawCRS).size()
                  == expected_crs_bytes[d],
              "Rexp CRS size regression");
    }
    rexp::RexpVerifyMetrics metrics;
    check(rexp::VerifyPrepared(setup.params, setup.statement, proof, &metrics),
          "prepared Rexp proof was rejected");
    check(rexp::VerifyReference(setup.params, setup.statement, proof),
          "reference Rexp proof was rejected");
    const auto validated = rexp::ValidateRexpProof(proof, d);
    check(rexp::VerifyValidatedProof(setup.params, setup.statement, validated),
          "validated Rexp proof was rejected");
    check(rexp::VerifyOptimized(setup.params, setup.statement, validated),
          "optimized Rexp proof was rejected");
    auto movable_proof = proof;
    const auto move_validated =
        rexp::ValidateRexpProof(std::move(movable_proof), d);
    check(rexp::VerifyOptimized(setup.params, setup.statement, move_validated),
          "move-validated Rexp proof was rejected");
    check(rexp::Verify(setup.rawCRS, setup.rawStatement, proof),
          "raw Rexp verification failed");

    const Bytes crs_wire = rexp::SerializeRexpCRSWire(setup.rawCRS);
    const auto decoded_crs = rexp::DeserializeRexpCRSWire(crs_wire);
    check(decoded_crs.d == setup.rawCRS.d && decoded_crs.n == setup.rawCRS.n
              && decoded_crs.Gamma == setup.rawCRS.Gamma
              && decoded_crs.Lambda == setup.rawCRS.Lambda,
          "CRS wire round trip differs");

    const Bytes statement_wire =
        rexp::SerializeRexpStatementWire(setup.rawStatement, setup.params.n());
    const auto decoded_statement =
        rexp::DeserializeRexpStatementWire(statement_wire, setup.params.n());
    check(decoded_statement.H == setup.rawStatement.H,
          "statement wire round trip differs");

    const Bytes proof_wire = rexp::SerializeRexpProofWire(proof, d);
    rexp::RexpProofValidationMetrics validation_metrics;
    const auto decoded_proof =
        rexp::DeserializeValidatedRexpProofWire(proof_wire, d, &validation_metrics);
    check(rexp::VerifyOptimized(setup.params, setup.statement, decoded_proof),
          "decoded Rexp proof was rejected");

    Bytes truncated = proof_wire;
    truncated.pop_back();
    check_throws(
        [&] { rexp::DeserializeRexpProofWire(truncated, d); },
        "truncated proof wire was accepted");
    Bytes trailing = proof_wire;
    trailing.push_back(0);
    check_throws(
        [&] { rexp::DeserializeRexpProofWire(trailing, d); },
        "proof wire with trailing byte was accepted");
    Bytes wrong_domain = proof_wire;
    if (wrong_domain.size() > 8) wrong_domain[8] ^= 1;
    check_throws(
        [&] { rexp::DeserializeRexpProofWire(wrong_domain, d); },
        "proof wire with wrong domain was accepted");
    check_throws(
        [&] { rexp::DeserializeRexpProofWire(proof_wire, d + 1); },
        "proof wire with wrong expected dimension was accepted");

    auto tampered = proof;
    rexp::G1::add(tampered.R, tampered.R, setup.params.Gamma().front());
    check(!rexp::VerifyPrepared(setup.params, setup.statement, tampered),
          "tampered Rexp proof was accepted");
}

void test_rexp_validation_failures() {
    const auto setup = rexp::Setup(2, "rexp-validation-seed");
    const auto proof = rexp::Prove(setup.params, setup.statement, setup.proverInput);

    auto malformed = proof;
    malformed.doryProofs.pop_back();
    check_throws(
        [&] { rexp::ValidateRexpProof(malformed, 2); },
        "malformed Rexp proof was validated");

    auto wrong_statement = setup.rawStatement;
    wrong_statement.H.pop_back();
    check_throws(
        [&] { rexp::PrepareStatement(setup.params, wrong_statement); },
        "wrong statement length was accepted");

    auto wrong_crs = setup.rawCRS;
    wrong_crs.Lambda.pop_back();
    check_throws(
        [&] { rexp::PreparePublicParameters(wrong_crs); },
        "wrong CRS length was accepted");
}

void test_wire_allocation_bounds() {
    Bytes crs_wire;
    append_frame(crs_wire, "REXP-G1-CRS-WIRE-BN254-V1");
    append_u64(crs_wire, 30);
    append_u64(crs_wire, std::uint64_t{1} << 30);
    append_u64(crs_wire, std::uint64_t{1} << 30);
    check_throws(
        [&] { rexp::DeserializeRexpCRSWire(crs_wire); },
        "impossible CRS element count reached allocation");

    Bytes statement_wire;
    append_frame(
        statement_wire, "REXP-G1-RAW-STATEMENT-WIRE-BN254-V1");
    append_u64(statement_wire, std::uint64_t{1} << 30);
    append_u64(statement_wire, std::uint64_t{1} << 30);
    check_throws(
        [&] {
            rexp::DeserializeRexpStatementWire(
                statement_wire, std::size_t{1} << 30);
        },
        "impossible statement element count reached allocation");
}

} 

int main() {
    try {
        rexp::initialize();
        test_challenge_vector();
        test_dory_setup_and_proofs();
        test_dory_batch();
        test_rexp_dimension(0);
        test_rexp_dimension(1);
        test_rexp_dimension(2);
        test_rexp_validation_failures();
        test_wire_allocation_bounds();
        std::cout << "All rexp tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rexp test failure: " << error.what() << '\n';
        return 1;
    }
}

#include "rexpbf/rexpbf.hpp"

#include <cstdint>
#include <iostream>
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

Bytes bytes(const std::string& text) {
    return {text.begin(), text.end()};
}

std::uint64_t read_u64(const Bytes& input, std::size_t offset) {
    check(offset <= input.size() && input.size() - offset >= 8,
          "test wire u64 is truncated");
    std::uint64_t value = 0;
    for (std::size_t i = 0; i < 8; ++i)
        value = (value << 8) | input[offset + i];
    return value;
}

void write_u64(Bytes& output, std::size_t offset, std::uint64_t value) {
    check(offset <= output.size() && output.size() - offset >= 8,
          "test wire u64 is truncated");
    for (std::size_t i = 0; i < 8; ++i)
        output[offset + i] = static_cast<std::uint8_t>(value >> (56 - 8 * i));
}

void skip_frame(const Bytes& input, std::size_t& offset) {
    const auto size = read_u64(input, offset);
    offset += 8;
    check(size <= input.size() - offset, "test wire frame is truncated");
    offset += static_cast<std::size_t>(size);
}

std::string hex(const rexpbf::Digest32& digest) {
    static constexpr char digits[] = "0123456789abcdef";
    std::string result;
    result.reserve(64);
    for (std::uint8_t byte : digest) {
        result.push_back(digits[byte >> 4]);
        result.push_back(digits[byte & 0x0f]);
    }
    return result;
}

rexpbf::SetupConfig config(std::size_t d) {
    return {
        d,
        bytes("REXPBF-TEST-V1|CRS"),
        bytes("REXPBF-TEST-V1|STATEMENT")};
}

void test_setup_and_proof_fixture() {
    const auto first = rexpbf::setup(config(3));
    const auto second = rexpbf::setup(config(3));
    check(first.crs.gamma == second.crs.gamma, "CRS G1 is not deterministic");
    check(first.crs.lambda == second.crs.lambda, "CRS G2 is not deterministic");
    check(first.crs.digest == second.crs.digest, "CRS digest is not deterministic");
    check(first.statement.h == second.statement.h, "statement is not deterministic");
    check(first.statement.digest == second.statement.digest,
          "statement digest is not deterministic");
    check(rexpbf::validate_crs(first.crs), "generated CRS is invalid");
    check(rexpbf::validate_crs_shape(first.crs), "generated CRS shape is invalid");
    check(rexpbf::validate_crs_points(first.crs), "generated CRS points are invalid");
    check(rexpbf::validate_crs_digest(first.crs), "generated CRS digest is invalid");
    check(rexpbf::validate_precomputation_shape(
              first.crs, first.precomputation),
          "generated precomputation shape is invalid");
    check(rexpbf::validate_precomputation_elements(first.precomputation),
          "generated precomputation elements are invalid");
    check(rexpbf::validate_statement_shape(first.crs, first.statement),
          "generated statement shape is invalid");
    check(rexpbf::validate_statement_elements(first.statement),
          "generated statement elements are invalid");
    check(rexpbf::validate_statement_digest(first.crs, first.statement),
          "generated statement digest is invalid");
    check(rexpbf::audit_precomputation(first.crs, first.precomputation),
          "generated precomputation audit failed");
    check(rexpbf::audit_statement(first.crs, first.statement),
          "generated statement audit failed");

    const auto proof = rexpbf::prove(
        first.crs, first.precomputation, first.statement,
        first.prover_input).proof;
    check(rexpbf::verify_reference(
              first.crs, first.precomputation, first.statement, proof),
          "reference verification failed");
    check(rexpbf::verify_online(
              first.crs, first.precomputation, first.statement, proof),
          "online verification failed");
    check(rexpbf::verify(
              first.crs, first.precomputation, first.statement, proof),
          "validated verification failed");
    const auto validated = rexpbf::validate_verification_inputs(
        first.crs, first.precomputation, first.statement, proof);
    check(validated.has_value(), "input validation failed");
    check(rexpbf::verify_prevalidated(*validated),
          "prevalidated verification failed");

    rexpbf::VerifyDiagnostics diagnostics;
    check(rexpbf::verify_with_diagnostics(
              first.crs, first.precomputation, first.statement, proof,
              diagnostics),
          "diagnostic verification failed");
    check(diagnostics.gt_terms_after_coalescing
              == diagnostics.gt_nonzero_scalar_terms,
          "diagnostic nonzero term count is inconsistent");
    check(diagnostics.gt_bases.size()
              == diagnostics.gt_terms_after_coalescing
              && diagnostics.gt_scalars.size()
                  == diagnostics.gt_terms_after_coalescing,
          "diagnostic GT vectors and term counts are inconsistent");
    check(diagnostics.gt_terms_after_coalescing
              <= diagnostics.gt_terms_before_coalescing,
          "diagnostic coalescing increased the term count");
    rexpbf::VerifyCoreBreakdown breakdown;
    check(rexpbf::verify_online_with_breakdown(
              first.crs, first.precomputation, first.statement, proof,
              breakdown),
          "profiled verification failed");
    check(breakdown.normalization_calls == 4 * (first.crs.d - 1) + 2,
          "profiled normalization count is inconsistent");
    check(breakdown.total_terms_before_normalization
              >= breakdown.total_terms_after_normalization
                   + breakdown.zero_terms_removed,
          "profiled normalization terms are inconsistent");

    const auto trace = rexpbf::replay_challenges(
        first.crs, first.statement, proof);
    check(trace.rho.size() == 3 && trace.beta.size() == 2
              && trace.alpha.size() == 2 && trace.gamma.size() == 2,
          "challenge trace shape mismatch");

    const Bytes proof_wire =
        rexpbf::serialize_proof_wire(proof, first.crs.d, first.crs.n);
    const std::string proof_hash = hex(rexpbf::sha256(proof_wire));
    const std::string expected_proof_hash =
        "18de22b27cabebe7d386afb4686d7464f890a6ee56771305c59825394069d50c";
    check(proof_hash == expected_proof_hash,
          "proof wire hash mismatch: " + proof_hash);
    const std::string trace_hash = hex(trace.final_digest);
    const std::string expected_trace_hash =
        "d449170e2c25632b72c4339d58547fa6451ab4d8ee3fac4cec9596a8bc0e810b";
    check(trace_hash == expected_trace_hash,
          "transcript digest mismatch: " + trace_hash);
}

void test_wire_round_trips() {
    const auto instance = rexpbf::setup(config(2));
    const auto proof = rexpbf::prove(
        instance.crs, instance.precomputation, instance.statement,
        instance.prover_input).proof;

    const Bytes crs_wire = rexpbf::serialize_crs_wire(instance.crs);
    const auto crs = rexpbf::deserialize_crs_wire(crs_wire);
    check(crs.gamma == instance.crs.gamma && crs.lambda == instance.crs.lambda
              && crs.digest == instance.crs.digest,
          "CRS wire round trip differs");

    const Bytes statement_wire =
        rexpbf::serialize_statement_wire(instance.crs, instance.statement);
    const auto statement =
        rexpbf::deserialize_statement_wire(statement_wire, instance.crs);
    check(statement.h == instance.statement.h
              && statement.digest == instance.statement.digest,
          "statement wire round trip differs");

    const Bytes proof_wire =
        rexpbf::serialize_proof_wire(proof, instance.crs.d, instance.crs.n);
    const auto decoded = rexpbf::deserialize_proof_wire(
        proof_wire, instance.crs.d, instance.crs.n);
    check(rexpbf::verify(
              instance.crs, instance.precomputation, instance.statement,
              decoded),
          "decoded proof was rejected");

    Bytes truncated = proof_wire;
    truncated.pop_back();
    check_throws(
        [&] {
            rexpbf::deserialize_proof_wire(
                truncated, instance.crs.d, instance.crs.n);
        },
        "truncated proof wire was accepted");
    Bytes trailing = proof_wire;
    trailing.push_back(0);
    check_throws(
        [&] {
            rexpbf::deserialize_proof_wire(
                trailing, instance.crs.d, instance.crs.n);
        },
        "proof wire with trailing bytes was accepted");
    Bytes wrong_domain = proof_wire;
    wrong_domain[8] ^= 1;
    check_throws(
        [&] {
            rexpbf::deserialize_proof_wire(
                wrong_domain, instance.crs.d, instance.crs.n);
        },
        "proof wire with wrong domain was accepted");

    Bytes truncated_crs = crs_wire;
    truncated_crs.pop_back();
    check_throws(
        [&] { rexpbf::deserialize_crs_wire(truncated_crs); },
        "truncated CRS wire was accepted");

    Bytes truncated_statement = statement_wire;
    truncated_statement.pop_back();
    check_throws(
        [&] {
            rexpbf::deserialize_statement_wire(
                truncated_statement, instance.crs);
        },
        "truncated statement wire was accepted");

    Bytes hostile_crs = crs_wire;
    std::size_t dimensions_offset = 0;
    skip_frame(hostile_crs, dimensions_offset);
    skip_frame(hostile_crs, dimensions_offset);
    skip_frame(hostile_crs, dimensions_offset);
    constexpr std::uint64_t hostile_d = 30;
    constexpr std::uint64_t hostile_n = std::uint64_t{1} << hostile_d;
    write_u64(hostile_crs, dimensions_offset, hostile_d);
    write_u64(hostile_crs, dimensions_offset + 8, hostile_n);
    write_u64(hostile_crs, dimensions_offset + 16, hostile_n);
    check_throws(
        [&] { rexpbf::deserialize_crs_wire(hostile_crs); },
        "CRS wire with an oversized element count was accepted");
}

void test_tampering() {
    const auto instance = rexpbf::setup(config(2));
    const auto proof = rexpbf::prove(
        instance.crs, instance.precomputation, instance.statement,
        instance.prover_input).proof;
    auto tampered = proof;
    rexpbf::GT::mul(
        tampered.steps[0].u, tampered.steps[0].u,
        instance.precomputation.x[0]);
    check(!rexpbf::verify(
              instance.crs, instance.precomputation, instance.statement,
              tampered),
          "tampered batch U was accepted");

    tampered = proof;
    rexpbf::G1::add(
        tampered.r_final, tampered.r_final, instance.crs.gamma[0]);
    check(!rexpbf::verify(
              instance.crs, instance.precomputation, instance.statement,
              tampered),
          "tampered final R was accepted");

    auto wrong_precomputation = instance.precomputation;
    wrong_precomputation.crs_digest[0] ^= 1;
    check(!rexpbf::verify(
              instance.crs, wrong_precomputation, instance.statement, proof),
          "wrong precomputation binding was accepted");

    auto owned_instance = rexpbf::setup(config(2));
    auto owned_proof = rexpbf::prove(
        owned_instance.crs, owned_instance.precomputation,
        owned_instance.statement, owned_instance.prover_input).proof;
    const auto token = rexpbf::validate_verification_inputs(
        owned_instance.crs, owned_instance.precomputation,
        owned_instance.statement, owned_proof);
    check(token.has_value(), "owning validation token was not created");
    owned_instance.crs.digest[0] ^= 1;
    owned_instance.precomputation.crs_digest[0] ^= 1;
    rexpbf::G1::add(
        owned_proof.r_final, owned_proof.r_final,
        owned_instance.crs.gamma[0]);
    check(rexpbf::verify_prevalidated(*token),
          "validated token changed after its source inputs were mutated");

    auto invalid_statement = instance.statement;
    invalid_statement.digest[0] ^= 1;
    check_throws(
        [&] {
            (void)rexpbf::replay_challenges(
                instance.crs, invalid_statement, proof);
        },
        "public challenge replay accepted an invalid statement");
}

void test_arithmetic_helpers() {
    const auto instance = rexpbf::setup(config(2));
    std::vector<rexpbf::Fr> scalars(3);
    scalars[0] = 2;
    scalars[1] = 3;
    scalars[2] = 5;
    std::vector<rexpbf::GT> bases {
        instance.precomputation.x[0],
        instance.precomputation.x[1],
        instance.precomputation.delta1_right[0]};
    check(rexpbf::gt_multiexp_pippenger(bases, scalars)
              == rexpbf::gt_multiexp_naive(bases, scalars),
          "GT multiexp implementations disagree");

    const auto inverses = rexpbf::batch_invert_nonzero(scalars);
    rexpbf::Fr one;
    one = 1;
    for (std::size_t i = 0; i < scalars.size(); ++i) {
        rexpbf::Fr product;
        rexpbf::Fr::mul(product, scalars[i], inverses[i]);
        check(product == one, "batch inversion is incorrect");
    }

    const auto gamma = std::span<const rexpbf::G1>(instance.crs.gamma).first(2);
    const auto lambda = std::span<const rexpbf::G2>(instance.crs.lambda).first(2);
    const auto product_pairing = rexpbf::pairing_product(gamma, lambda);
    rexpbf::GT first, second, expected;
    mcl::bn::pairing(first, gamma[0], lambda[0]);
    mcl::bn::pairing(second, gamma[1], lambda[1]);
    rexpbf::GT::mul(expected, first, second);
    check(product_pairing == expected, "pairing product differs");
}

void test_baseline_sizes() {
    static constexpr std::size_t proof_sizes[] = {237, 4549, 8861, 13173};
    static constexpr std::size_t crs_sizes[] = {355, 579, 1027, 1923};
    for (std::size_t d = 1; d <= 4; ++d) {
        const auto instance = rexpbf::setup(config(d));
        const auto proof = rexpbf::prove(
            instance.crs, instance.precomputation, instance.statement,
            instance.prover_input).proof;
        check(rexpbf::serialize_proof_wire(proof, d, instance.crs.n).size()
                  == proof_sizes[d - 1],
              "proof size regression");
        check(rexpbf::serialize_crs_wire(instance.crs).size()
                  == crs_sizes[d - 1],
              "CRS size regression");
    }
}

} // namespace

int main() {
    try {
        rexpbf::initialize_bn254();
        test_setup_and_proof_fixture();
        test_wire_round_trips();
        test_tampering();
        test_arithmetic_helpers();
        test_baseline_sizes();
        std::cout << "All rexpbf tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rexpbf test failure: " << error.what() << '\n';
        return 1;
    }
}

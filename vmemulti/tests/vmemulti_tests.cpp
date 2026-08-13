#include "vme_ibf/serialization.hpp"
#include "vme_ibf/vmemulti.hpp"
#include "vme_ibf/verify_online.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace vme_ibf;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Function>
void check_throws(Function&& function, const std::string& message) {
    try { function(); }
    catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}

std::vector<G2> public_h(std::size_t d) {
    const std::size_t n = std::size_t{1} << d;
    std::vector<G2> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto seed = "VMEMULTI/TEST/H/" + std::to_string(d)
            + "/" + std::to_string(i);
        G2 point;
        mcl::bn::hashAndMapToG2(point, seed.data(), seed.size());
        result.push_back(point);
    }
    return result;
}

VmeMultiStatement make_statement(const VmeIbfCRS& crs, std::size_t m,
                                 DeterministicRng& rng) {
    VmeMultiStatement result;
    result.x_instances.resize(m);
    result.X_instances.reserve(m);
    for (auto& instance : result.x_instances) {
        instance.reserve(crs.n);
        for (std::size_t j = 0; j < crs.n; ++j)
            instance.push_back(rng.random_fr());
        result.X_instances.push_back(g2_multiexp(crs.H, instance));
    }
    return result;
}

struct Fixture {
    SetupResult setup;
    VmeMultiStatement statement;
    VmeMultiProverTrace trace;
    VmeIbfProof proof;
};

Fixture fixture(std::size_t d, std::size_t m) {
    DeterministicRng setup_rng(
        "VMEMULTI/TEST/SETUP/d=" + std::to_string(d));
    Fixture result;
    result.setup = setup(d, public_h(d), setup_rng);
    DeterministicRng statement_rng(
        "VMEMULTI/TEST/STATEMENT/d=" + std::to_string(d)
        + "/m=" + std::to_string(m));
    result.statement = make_statement(result.setup.crs, m, statement_rng);
    result.trace = prove_vmemulti_with_trace(
        result.setup.crs, result.setup.precomp, result.statement);
    result.proof = assemble_public_proof(result.trace.phase1,
                                         result.trace.phase2);
    return result;
}

Digest aggregate_digest(const AggregatedInstance& aggregate) {
    Bytes bytes;
    append_frame(bytes, "VMEMULTI/TEST/AGGREGATE/V1");
    for (const auto& scalar : aggregate.z)
        append_frame(bytes, serialize(scalar));
    append_frame(bytes, serialize(aggregate.Y));
    return sha256(bytes);
}

void test_golden_protocol() {
    const auto first = fixture(3, 4);
    const auto second = fixture(3, 4);
    const auto binding = compute_vmemulti_binding_digest(
        first.setup.crs, first.setup.precomp, first.statement);
    Digest after_gamma;
    const auto gamma = derive_vmemulti_aggregate_gamma(
        first.setup.crs, first.setup.precomp, first.statement, &after_gamma);

    check(first.setup.crs.digest == second.setup.crs.digest,
          "CRS is not deterministic");
    check(binding == compute_vmemulti_binding_digest(
              second.setup.crs, second.setup.precomp, second.statement),
          "binding digest is not deterministic");
    check(serialize(gamma) == serialize(second.trace.aggregate_gamma),
          "aggregate challenge is not deterministic");
    check(aggregate_digest(first.trace.aggregate)
              == aggregate_digest(second.trace.aggregate),
          "aggregate is not deterministic");
    check(serialize_proof(3, first.proof)
              == serialize_proof(3, second.proof),
          "proof is not deterministic");

    check(verify_vmemulti_diagnostic(
              first.setup.crs, first.setup.precomp,
              first.statement, first.proof).accepted,
          "reference verifier rejected fixture");
    const auto combined = verify_vmemulti_combined_diagnostic(
        first.setup.crs, first.setup.precomp,
        first.statement, first.proof);
    check(combined.accepted
              && combined.evaluated_combined_reference
                  == combined.evaluated_combined_residual
              && combined.gt_terms_after_normalize
                  <= combined.gt_terms_before_normalize
              && combined.pairing_terms_after_normalize
                  <= combined.pairing_terms_before_normalize
              && combined.gt_multiexp_calls == 1
              && combined.pairing_product_calls == 1
              && combined.miller_loop_batches <= 1
              && combined.final_exponentiations <= 1,
          "combined reference/optimized evaluation disagrees");
    check(verify_vmemulti(first.setup.crs, first.setup.precomp,
                          first.statement, first.proof),
          "online verifier rejected fixture");
    const auto validated = validate_vmemulti_inputs(
        first.setup.crs, first.setup.precomp,
        first.statement, first.proof);
    check(validated.has_value(), "online token creation failed");
    const auto online = verify_vmemulti_online_with_trace(*validated);
    check(online.accepted && online.total_ms >= 0
              && online.transcript_prefix_ms >= 0
              && online.aggregate_gamma_ms >= 0
              && online.aggregate_z_ms >= 0
              && online.aggregate_Y_ms >= 0
              && online.core.total_ms >= 0,
          "multi online timing trace is inconsistent");

    const auto crs_hash = hex(sha256(serialize_crs(first.setup.crs)));
    const auto precomputation_hash = hex(sha256(serialize_precomputation(
        first.setup.crs, first.setup.precomp)));
    const auto binding_hash = hex(binding);
    const auto gamma_hash = hex(sha256(serialize(gamma)));
    const auto after_gamma_hash = hex(after_gamma);
    const auto aggregate_hash = hex(aggregate_digest(first.trace.aggregate));
    const auto proof_hash = hex(sha256(serialize_proof(3, first.proof)));
    const auto transcript_hash = hex(first.trace.phase2.final_transcript_digest);
    check(crs_hash
              == "385c2c3f29eac1d3979501ef8c8477dfc012f4527d8a9f151e09935205967839",
          "CRS golden changed: " + crs_hash);
    check(precomputation_hash
              == "3b5fc269720f81c45fc17149a08e1313e99acb740cab575214f88262d4ac3464",
          "precomputation golden changed: " + precomputation_hash);
    check(binding_hash
              == "9e494bcb3661472153ff39e5e1ff06699b406fbe73d1da7cd9d0f2ec155bbb61",
          "binding golden changed: " + binding_hash);
    check(gamma_hash
              == "0fbb1787eb77081c0a83671dc5febfbae5f63e8db3c35e317df8ddc695fa6bb7",
          "aggregate challenge golden changed: " + gamma_hash);
    check(after_gamma_hash
              == "7580ecad031acecac14390bc4bfc3f2b34efde8beaff8278dcb21de174e23010",
          "post-aggregation transcript golden changed: " + after_gamma_hash);
    check(aggregate_hash
              == "cf0d8783555adf741193ae74cdcd7891f46cc23f3efff195aebaccd9bd3ae142",
          "aggregate golden changed: " + aggregate_hash);
    check(proof_hash
              == "83a11dd012c0047bd6f7b3c40921920e71d57e01ef23611828214b1c1c0068f7",
          "proof golden changed: " + proof_hash);
    check(transcript_hash
              == "930bd71952e5f69ad6f6e3557c0964cf36e6c4089d53d1ebfa033bc78b8b8b9d",
          "transcript golden changed: " + transcript_hash);
}

void test_size_and_instance_baselines() {
    constexpr std::size_t proof_sizes[]{2832, 7056, 11280};
    constexpr std::size_t crs_sizes[]{304, 496, 880};
    for (std::size_t d = 1; d <= 3; ++d) {
        const auto one = fixture(d, 1);
        const auto many = fixture(d, 4);
        check(serialize_proof(d, one.proof).size() == proof_sizes[d - 1]
                  && serialize_proof(d, many.proof).size()
                      == proof_sizes[d - 1],
              "proof size baseline changed");
        check(serialize_crs(one.setup.crs).size() == crs_sizes[d - 1],
              "CRS size baseline changed");
        check(verify_vmemulti(one.setup.crs, one.setup.precomp,
                              one.statement, one.proof)
                  && verify_vmemulti(many.setup.crs, many.setup.precomp,
                                     many.statement, many.proof),
              "dimension/instance baseline rejected");
    }
}

void test_owning_validation_boundary() {
    auto value = fixture(2, 3);
    const auto validated = validate_vmemulti_inputs(
        value.setup.crs, value.setup.precomp, value.statement, value.proof);
    check(validated.has_value(), "validated multi token was not created");
    value.setup.crs.digest[0] ^= 1;
    value.setup.precomp.pairing_x.clear();
    value.statement.x_instances.clear();
    value.proof.batch_U.clear();
    check(verify_vmemulti_online(*validated),
          "validated multi token changed after source mutation");
}

void test_object_validation() {
    const auto value = fixture(2, 2);
    check(validate_crs(value.setup.crs)
              && validate_precomputation_shape(
                  value.setup.crs, value.setup.precomp)
              && validate_precomputation_elements(value.setup.precomp)
              && audit_precomputation(value.setup.crs, value.setup.precomp),
          "valid setup objects were rejected");
    auto crs = value.setup.crs;
    crs.digest[0] ^= 1;
    check(!validate_crs(crs), "bad CRS digest was accepted");
    auto precomputation = value.setup.precomp;
    GT::mul(precomputation.pairing_x[0], precomputation.pairing_x[0],
            value.setup.precomp.pairing_x[1]);
    check(validate_precomputation_elements(precomputation)
              && !audit_precomputation(value.setup.crs, precomputation),
          "inconsistent precomputation was accepted");
}

void test_phase2_validation() {
    const auto value = fixture(2, 2);

    auto bad_precomputation = value.setup.precomp;
    GT::mul(bad_precomputation.pairing_x[0],
            bad_precomputation.pairing_x[0],
            value.setup.precomp.pairing_x[1]);
    check_throws([&] {
        (void)prove_phase2(
            value.setup.crs, bad_precomputation, value.trace.phase1);
    }, "phase 2 accepted inconsistent precomputation");

    auto bad_phase1 = value.trace.phase1;
    bad_phase1.fresh[1].D0.clear();
    check_throws([&] {
        (void)prove_phase2(
            value.setup.crs, value.setup.precomp, bad_phase1);
    }, "phase 2 accepted an invalid fresh instance");

    bad_phase1 = value.trace.phase1;
    bad_phase1.dynamic_claims[0].E.clear();
    check_throws([&] {
        (void)prove_phase2(
            value.setup.crs, value.setup.precomp, bad_phase1);
    }, "phase 2 accepted an invalid phase-1 claim");

    bad_phase1 = value.trace.phase1;
    bad_phase1.transcript_start[0] ^= 1;
    check_throws([&] {
        (void)prove_phase2(
            value.setup.crs, value.setup.precomp, bad_phase1);
    }, "phase 2 accepted a mismatched transcript start state");
}

void test_aggregation_reference() {
    const auto value = fixture(2, 4);
    const Fr gamma = value.trace.aggregate_gamma;
    std::vector<Fr> expected_z(value.setup.crs.n);
    for (auto& scalar : expected_z) scalar.clear();
    std::vector<Fr> powers(value.statement.x_instances.size());
    Fr power;
    power = 1;
    for (std::size_t i = 0; i < powers.size(); ++i) {
        powers[i] = power;
        for (std::size_t j = 0; j < expected_z.size(); ++j) {
            Fr term;
            Fr::mul(term, power, value.statement.x_instances[i][j]);
            Fr::add(expected_z[j], expected_z[j], term);
        }
        Fr::mul(power, power, gamma);
    }
    check(value.trace.aggregate.z == expected_z
              && value.trace.aggregate.Y
                  == g2_multiexp_reference(
                      value.statement.X_instances, powers),
          "optimized aggregation disagrees with reference");

    auto invalid = value.statement;
    invalid.x_instances[0].pop_back();
    check(!validate_vmemulti_statement(value.setup.crs, invalid),
          "ragged statement matrix was accepted");
    invalid = value.statement;
    invalid.x_instances.resize(kMaxVmeMultiInstances + 1,
                               invalid.x_instances.front());
    invalid.X_instances.resize(kMaxVmeMultiInstances + 1,
                               invalid.X_instances.front());
    check(!validate_vmemulti_statement(value.setup.crs, invalid),
          "oversized instance set was accepted");
}

void test_wire_and_tampering() {
    const auto value = fixture(2, 3);
    DecodeError error = DecodeError::None;
    VmeIbfCRS crs;
    const auto crs_wire = serialize_crs(value.setup.crs);
    check(deserialize_crs(crs_wire, crs, &error)
              && crs.digest == value.setup.crs.digest,
          "CRS round trip failed");
    VmeIbfPrecomputation precomputation;
    const auto pre_wire = serialize_precomputation(
        value.setup.crs, value.setup.precomp);
    check(deserialize_precomputation(
              pre_wire, value.setup.crs, precomputation, &error)
              && audit_precomputation(value.setup.crs, precomputation),
          "precomputation round trip failed");
    VmeIbfProof proof;
    const auto proof_wire = serialize_proof(2, value.proof);
    check(deserialize_proof(proof_wire, value.setup.crs, proof, &error),
          "proof round trip failed");
    check(verify_vmemulti(value.setup.crs, value.setup.precomp,
                          value.statement, proof),
          "decoded proof failed verification");

    auto truncated = proof_wire;
    truncated.pop_back();
    proof = value.proof;
    check(!deserialize_proof(truncated, value.setup.crs, proof, &error)
              && error == DecodeError::Truncated
              && serialize_proof(2, proof) == proof_wire,
          "truncated decode was accepted or modified output");
    auto trailing = proof_wire;
    trailing.push_back(0);
    check(!deserialize_proof(trailing, value.setup.crs, proof, &error)
              && error == DecodeError::TrailingBytes,
          "trailing bytes were accepted");
    auto wrong_magic = proof_wire;
    wrong_magic[0] ^= 1;
    check(!deserialize_proof(wrong_magic, value.setup.crs, proof, &error)
              && error == DecodeError::WrongMagic,
          "wrong magic was accepted");
    auto wrong_version = proof_wire;
    wrong_version[9] ^= 1;
    check(!deserialize_proof(wrong_version, value.setup.crs, proof, &error)
              && error == DecodeError::UnsupportedVersion,
          "unsupported version was accepted");
    auto wrong_curve = proof_wire;
    wrong_curve[11] ^= 1;
    check(!deserialize_proof(wrong_curve, value.setup.crs, proof, &error)
              && error == DecodeError::WrongCurve,
          "wrong curve was accepted");
    auto wrong_dimension = proof_wire;
    wrong_dimension[15] = 3;
    check(!deserialize_proof(
              wrong_dimension, value.setup.crs, proof, &error)
              && error == DecodeError::InvalidDimension,
          "wrong dimension was accepted");
    auto invalid_crs = value.setup.crs;
    invalid_crs.digest[0] ^= 1;
    precomputation = value.setup.precomp;
    check(!deserialize_precomputation(
              pre_wire, invalid_crs, precomputation, &error)
              && error == DecodeError::InvalidCrs
              && precomputation.pairing_x
                  == value.setup.precomp.pairing_x,
          "precomputation decoder accepted an invalid CRS or modified output");
    VmeIbfStatement core_statement = value.trace.phase1.statement;
    const auto statement_wire = serialize_statement(
        value.setup.crs, core_statement);
    check(!deserialize_statement(
              statement_wire, invalid_crs, core_statement, &error)
              && error == DecodeError::InvalidCrs
              && core_statement.digest
                  == value.trace.phase1.statement.digest,
          "statement decoder accepted an invalid CRS or modified output");
    check(!deserialize_proof(proof_wire, invalid_crs, proof, &error)
              && error == DecodeError::InvalidCrs,
          "proof decoder accepted an invalid CRS");
    check_throws([&] { (void)serialize_crs(invalid_crs); },
                 "serializer accepted an invalid CRS");

    proof = value.proof;
    GT::mul(proof.batch_U[0], proof.batch_U[0],
            value.setup.precomp.pairing_x[0]);
    check(!verify_vmemulti_diagnostic(
              value.setup.crs, value.setup.precomp,
              value.statement, proof).accepted
              && !verify_vmemulti(
                  value.setup.crs, value.setup.precomp,
                  value.statement, proof),
          "tampered proof was accepted");
}

} // namespace

int main() {
    try {
        initialize();
        test_golden_protocol();
        test_size_and_instance_baselines();
        test_owning_validation_boundary();
        test_object_validation();
        test_phase2_validation();
        test_aggregation_reference();
        test_wire_and_tampering();
        std::cout << "All vmemulti tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "vmemulti test failure: " << error.what() << '\n';
        return 1;
    }
}

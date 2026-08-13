#include "vme_ibf/batch_inversion.hpp"
#include "vme_ibf/gt_multiexp.hpp"
#include "vme_ibf/phase1.hpp"
#include "vme_ibf/phase2.hpp"
#include "vme_ibf/proof.hpp"
#include "vme_ibf/serialization.hpp"
#include "vme_ibf/setup.hpp"
#include "vme_ibf/verify_combined.hpp"
#include "vme_ibf/verify_deferred.hpp"
#include "vme_ibf/verify_online.hpp"
#include "vme_ibf/verify_reference.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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
        const auto seed = "VME.BF.G2/TEST/H/" + std::to_string(d)
            + "/" + std::to_string(i);
        G2 point;
        mcl::bn::hashAndMapToG2(point, seed.data(), seed.size());
        result.push_back(point);
    }
    return result;
}

struct Fixture {
    SetupResult setup;
    Phase1Result phase1;
    Phase2Result phase2;
    VmeIbfProof proof;
};

Fixture fixture(std::size_t d) {
    DeterministicRng rng("VME.BF.G2/TEST/V1/d=" + std::to_string(d));
    Fixture result;
    result.setup = vme_ibf::setup(d, public_h(d), rng);
    result.phase1 = prove_phase1(
        result.setup.crs, result.setup.precomp, result.setup.statement_input);
    result.phase2 = prove_phase2(
        result.setup.crs, result.setup.precomp, result.phase1);
    result.proof = assemble_public_proof(result.phase1, result.phase2);
    return result;
}

void test_determinism_and_protocol() {
    const auto first = fixture(3);
    const auto second = fixture(3);
    check(first.setup.crs.digest == second.setup.crs.digest,
          "CRS digest is not deterministic");
    check(validate_crs_shape(first.setup.crs)
              && validate_crs_elements(first.setup.crs)
              && validate_crs_digest(first.setup.crs)
              && validate_crs(first.setup.crs),
          "generated CRS validation failed");
    check(validate_precomputation_shape(first.setup.crs, first.setup.precomp)
              && validate_precomputation_elements(first.setup.precomp)
              && audit_precomputation(first.setup.crs, first.setup.precomp),
          "generated precomputation validation failed");
    check(first.phase1.statement.digest == second.phase1.statement.digest,
          "statement digest is not deterministic");
    check(validate_statement_shape(first.setup.crs, first.phase1.statement)
              && validate_statement_elements(first.phase1.statement)
              && validate_statement_digest(first.setup.crs, first.phase1.statement),
          "generated statement validation failed");
    check(serialize_proof(3, first.proof) == serialize_proof(3, second.proof),
          "proof is not deterministic");
    check(first.phase2.final_transcript_digest
              == second.phase2.final_transcript_digest,
          "transcript is not deterministic");

    check(verify_reference(first.setup.crs, first.setup.precomp,
                           first.phase1.statement, first.proof),
          "reference verifier rejected fixture");
    check(verify_deferred(first.setup.crs, first.setup.precomp,
                          first.phase1.statement, first.proof),
          "deferred verifier rejected fixture");
    check(verify_deferred_combined(first.setup.crs, first.setup.precomp,
                                   first.phase1.statement, first.proof),
          "combined verifier rejected fixture");
    const auto deferred_trace = verify_deferred_with_trace(
        first.setup.crs, first.setup.precomp,
        first.phase1.statement, first.proof);
    check(deferred_trace.accepted
              && deferred_trace.dory_residual_reference
                  == deferred_trace.dory_residual_pippenger
              && deferred_trace.rexp_residual_reference
                  == deferred_trace.rexp_residual_pippenger,
          "deferred reference and optimized traces disagree");
    check(deferred_trace.gt_multiexp_calls == 2
              && deferred_trace.pairing_product_calls == 2,
          "deferred operation counters are inconsistent");
    const auto combined_trace = verify_deferred_combined_with_trace(
        first.setup.crs, first.setup.precomp,
        first.phase1.statement, first.proof);
    check(combined_trace.accepted
              && combined_trace.evaluated_combined_reference
                  == combined_trace.evaluated_combined_residual,
          "combined reference and optimized traces disagree");
    check(combined_trace.gt_terms_after_normalize
              <= combined_trace.gt_terms_before_normalize
              && combined_trace.pairing_terms_after_normalize
                  <= combined_trace.pairing_terms_before_normalize
              && combined_trace.gt_multiexp_calls == 1
              && combined_trace.pairing_product_calls == 1
              && combined_trace.miller_loop_batches <= 1
              && combined_trace.final_exponentiations <= 1,
          "combined trace counters are inconsistent");
    const auto inputs = validate_verification_inputs(
        first.setup.crs, first.setup.precomp,
        first.phase1.statement, first.proof);
    check(inputs.has_value(),
          "online input preparation failed");
    check(verify_online(*inputs), "online verifier rejected fixture");
    const auto online_trace = verify_online_with_trace(*inputs);
    check(online_trace.accepted && online_trace.total_ms >= 0
              && online_trace.gt_msm_pairing_ms >= 0,
          "online trace is inconsistent");

    const auto proof_hash = hex(sha256(serialize_proof(3, first.proof)));
    const auto crs_hash = hex(sha256(serialize_crs(first.setup.crs)));
    const auto precomputation_hash = hex(sha256(
        serialize_precomputation(first.setup.crs, first.setup.precomp)));
    const auto statement_hash = hex(sha256(
        serialize_statement(first.setup.crs, first.phase1.statement)));
    const auto transcript_hash = hex(first.phase2.final_transcript_digest);
    check(proof_hash
              == "6770f8534c21510208a094aac58956b6a9fd3efcb605f246f2450410bf8126a9",
          "proof hash baseline changed: " + proof_hash);
    check(crs_hash
              == "19b24005db27d7fa5609045f52f5d635c26d1f614706b4f9ad21485a3dc2dd04",
          "CRS hash baseline changed: " + crs_hash);
    check(precomputation_hash
              == "84dbed4b1334d8ab579986bd35cc6b429a53bc1cb52cba744f77ae514b4abb83",
          "precomputation hash baseline changed: " + precomputation_hash);
    check(statement_hash
              == "6531cd8b593b0056d855da4bcdfa716069e00d2130bb471def4768a341906112",
          "statement hash baseline changed: " + statement_hash);
    check(transcript_hash
              == "b27657686dcc772bb2da5a2618a347cece1961ab5b41e6df2ede8aa0140f6166",
          "transcript hash baseline changed: " + transcript_hash);
}

void test_wire_round_trips() {
    const auto f = fixture(2);
    DecodeError error = DecodeError::None;
    VmeIbfCRS crs;
    check(deserialize_crs(serialize_crs(f.setup.crs), crs, &error),
          "CRS round trip failed");
    check(crs.digest == f.setup.crs.digest, "CRS round trip digest changed");
    VmeIbfPrecomputation precomputation;
    check(deserialize_precomputation(
              serialize_precomputation(f.setup.crs, f.setup.precomp),
              f.setup.crs, precomputation, &error),
          "precomputation round trip failed");
    check(validate_precomputation(f.setup.crs, precomputation),
          "decoded precomputation audit failed");
    VmeIbfStatement statement;
    check(deserialize_statement(
              serialize_statement(f.setup.crs, f.phase1.statement),
              f.setup.crs, statement, &error),
          "statement round trip failed");
    VmeIbfProof proof;
    const auto proof_wire = serialize_proof(2, f.proof);
    check(deserialize_proof(proof_wire, f.setup.crs, proof, &error),
          "proof round trip failed");
    check(verify_reference(f.setup.crs, f.setup.precomp, statement, proof),
          "decoded objects failed verification");

    auto truncated = proof_wire;
    truncated.pop_back();
    check(!deserialize_proof(truncated, f.setup.crs, proof, &error)
              && error == DecodeError::Truncated,
          "truncated proof was accepted or misclassified");
    auto trailing = proof_wire;
    trailing.push_back(0);
    check(!deserialize_proof(trailing, f.setup.crs, proof, &error)
              && error == DecodeError::TrailingBytes,
          "trailing proof bytes were accepted or misclassified");
    auto wrong_magic = proof_wire;
    wrong_magic[0] ^= 1;
    check(!deserialize_proof(wrong_magic, f.setup.crs, proof, &error)
              && error == DecodeError::WrongMagic,
          "wrong proof magic was accepted or misclassified");
    auto wrong_version = proof_wire;
    wrong_version[9] ^= 1;
    check(!deserialize_proof(wrong_version, f.setup.crs, proof, &error)
              && error == DecodeError::UnsupportedVersion,
          "wrong proof version was accepted or misclassified");
    auto wrong_curve = proof_wire;
    wrong_curve[11] ^= 1;
    check(!deserialize_proof(wrong_curve, f.setup.crs, proof, &error)
              && error == DecodeError::WrongCurve,
          "wrong proof curve was accepted or misclassified");
    auto wrong_dimension = proof_wire;
    wrong_dimension[15] ^= 1;
    check(!deserialize_proof(wrong_dimension, f.setup.crs, proof, &error)
              && error == DecodeError::InvalidDimension,
          "wrong proof dimension was accepted or misclassified");
    const auto unchanged_wire = serialize_proof(2, f.proof);
    proof = f.proof;
    check(!deserialize_proof(truncated, f.setup.crs, proof, &error)
              && serialize_proof(2, proof) == unchanged_wire,
          "failed proof decode modified its output object");

    auto invalid_crs = f.setup.crs;
    invalid_crs.digest[0] ^= 1;
    check_throws([&] { (void)serialize_crs(invalid_crs); },
                 "serializer accepted a CRS with a wrong digest");
    check(!deserialize_precomputation(
              serialize_precomputation(f.setup.crs, f.setup.precomp),
              invalid_crs, precomputation, &error)
              && error == DecodeError::InvalidCrs,
          "precomputation decoder accepted an invalid CRS");
    check(!deserialize_statement(
              serialize_statement(f.setup.crs, f.phase1.statement),
              invalid_crs, statement, &error)
              && error == DecodeError::InvalidCrs,
          "statement decoder accepted an invalid CRS");
    check(!deserialize_proof(proof_wire, invalid_crs, proof, &error)
              && error == DecodeError::InvalidCrs,
          "proof decoder accepted an invalid CRS");
    auto invalid_statement = f.phase1.statement;
    invalid_statement.digest[0] ^= 1;
    check_throws([&] {
        (void)serialize_statement(f.setup.crs, invalid_statement);
    }, "serializer accepted a statement with a wrong digest");
}

void test_tampering() {
    const auto f = fixture(2);
    auto proof = f.proof;
    GT::mul(proof.batch_U[0], proof.batch_U[0], f.setup.precomp.pairing_x[0]);
    check(!verify_reference(f.setup.crs, f.setup.precomp,
                            f.phase1.statement, proof),
          "reference verifier accepted tampered U");
    check(!verify_deferred_combined(f.setup.crs, f.setup.precomp,
                                    f.phase1.statement, proof),
          "combined verifier accepted tampered U");

    proof = f.proof;
    G1::add(proof.R, proof.R, f.setup.crs.G[0]);
    check(!verify_reference(f.setup.crs, f.setup.precomp,
                            f.phase1.statement, proof),
          "reference verifier accepted tampered R");

    auto mutable_fixture = fixture(2);
    const auto token = validate_verification_inputs(
        mutable_fixture.setup.crs, mutable_fixture.setup.precomp,
        mutable_fixture.phase1.statement, mutable_fixture.proof);
    check(token.has_value(), "validation token was not created");
    mutable_fixture.setup.crs.digest[0] ^= 1;
    mutable_fixture.setup.precomp.pairing_x.clear();
    mutable_fixture.phase1.statement.digest[0] ^= 1;
    mutable_fixture.proof.batch_U.clear();
    check(verify_online(*token),
          "validated token changed after source mutation");

    auto bad_precomputation = f.setup.precomp;
    GT::mul(bad_precomputation.pairing_x[0],
            bad_precomputation.pairing_x[0],
            f.setup.precomp.pairing_x[1]);
    check_throws([&] {
        (void)prove_phase2(f.setup.crs, bad_precomputation, f.phase1);
    }, "phase 2 accepted inconsistent precomputation");

    auto bad_phase1 = f.phase1;
    bad_phase1.fresh[1].D0.clear();
    check_throws([&] {
        (void)prove_phase2(f.setup.crs, f.setup.precomp, bad_phase1);
    }, "phase 2 accepted an invalid fresh instance");

    bad_phase1 = f.phase1;
    bad_phase1.dynamic_claims[0].E.clear();
    check_throws([&] {
        (void)prove_phase2(f.setup.crs, f.setup.precomp, bad_phase1);
    }, "phase 2 accepted an invalid phase-1 claim");
}

void test_arithmetic() {
    const auto f = fixture(2);
    std::vector<Fr> scalars(3);
    scalars[0] = 2;
    scalars[1] = 3;
    scalars[2] = 5;
    std::vector<GT> bases{f.setup.precomp.pairing_x[0],
                          f.setup.precomp.pairing_x[1],
                          f.setup.precomp.delta1R[0]};
    check(gt_multiexp_reference(bases, scalars)
              == gt_multiexp_pippenger(bases, scalars),
          "GT multiexponentiation implementations disagree");
    set_gt_multiexp_window_override_for_benchmark(3);
    check(gt_multiexp_reference(bases, scalars)
              == gt_multiexp_pippenger(bases, scalars),
          "GT multiexponentiation window override changed the result");
    set_gt_multiexp_window_override_for_benchmark(0);
    check_throws([] { set_gt_multiexp_window_override_for_benchmark(9); },
                 "invalid GT multiexponentiation window was accepted");
    const auto inverses = batch_invert_nonzero(scalars);
    Fr one;
    one = 1;
    for (std::size_t i = 0; i < scalars.size(); ++i) {
        Fr product;
        Fr::mul(product, scalars[i], inverses[i]);
        check(product == one, "batch inversion is incorrect");
    }
    check(g1_multiexp_reference(f.setup.crs.G, f.setup.statement_input.x)
              == g1_multiexp(f.setup.crs.G, f.setup.statement_input.x),
          "G1 multiexponentiation implementations disagree");
    check(g2_multiexp_reference(f.setup.crs.H, f.setup.statement_input.x)
              == g2_multiexp(f.setup.crs.H, f.setup.statement_input.x),
          "G2 multiexponentiation implementations disagree");
}

void test_size_baselines() {
    constexpr std::size_t proof_sizes[]{2832, 7056, 11280, 15504};
    constexpr std::size_t crs_sizes[]{304, 496, 880, 1648};
    for (std::size_t d = 1; d <= 4; ++d) {
        const auto f = fixture(d);
        check(serialize_proof(d, f.proof).size() == proof_sizes[d - 1],
              "proof size baseline changed");
        check(serialize_crs(f.setup.crs).size() == crs_sizes[d - 1],
              "CRS size baseline changed");
    }
}

} // namespace

int main() {
    try {
        initialize();
        test_determinism_and_protocol();
        test_wire_round_trips();
        test_tampering();
        test_arithmetic();
        test_size_baselines();
        std::cout << "All vmebf tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "vmebf test failure: " << error.what() << '\n';
        return 1;
    }
}

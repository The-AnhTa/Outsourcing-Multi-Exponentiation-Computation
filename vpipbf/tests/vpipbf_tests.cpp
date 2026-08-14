#include "vpip_bf/protocol.hpp"
#include "vpip_bf/serialization.hpp"
#include "vpip_bf/verify_deferred.hpp"
#include "vpip_bf/verify_reference.hpp"
#include "vpip_bf/batch_inversion.hpp"
#include "vpip_bf/gt_multiexp.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using namespace vpip_bf;

void check(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}
template<class Function>
void check_throws(Function&& function, const std::string& message) {
    try { function(); } catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}

std::vector<G2> public_lambda(std::size_t d) {
    const std::size_t n = std::size_t{1} << d;
    G2 base;
    const std::string tag = "VPIPBF/TEST/LAMBDA/d=" + std::to_string(d);
    mcl::bn::hashAndMapToG2(base, tag.data(), tag.size());
    std::vector<G2> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        Fr scalar;
        scalar = static_cast<int>(i + 1);
        result.push_back(g2_pow(base, scalar));
    }
    return result;
}

std::vector<G1> public_x(std::size_t d) {
    const std::size_t n = std::size_t{1} << d;
    G1 base;
    const std::string tag = "VPIPBF/TEST/X/d=" + std::to_string(d);
    mcl::bn::hashAndMapToG1(base, tag.data(), tag.size());
    std::vector<G1> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        Fr scalar;
        scalar = static_cast<int>(i + 3);
        result.push_back(g1_pow(base, scalar));
    }
    return result;
}

struct Fixture {
    SetupResult setup;
    ProveResult proved;
};

Fixture fixture(std::size_t d) {
    DeterministicRng rng("VPIPBF/TEST/SETUP/d=" + std::to_string(d));
    Fixture value;
    value.setup = setup(d, public_lambda(d), public_x(d), rng);
    value.proved = Prove(value.setup.crs, value.setup.precomp,
                         value.setup.statement_input);
    return value;
}

void test_baseline() {
    const auto first = fixture(3);
    const auto second = fixture(3);
    check(first.setup.crs.digest == second.setup.crs.digest,
          "deterministic CRS changed");
    check(first.setup.precomp.digest == second.setup.precomp.digest,
          "deterministic precomputation changed");
    check(first.proved.statement.digest == second.proved.statement.digest,
          "deterministic statement changed");
    check(serialize_proof(3, first.proved.proof)
              == serialize_proof(3, second.proved.proof),
          "deterministic proof changed");
    check(Verify(first.setup.crs, first.setup.precomp,
                 first.proved.statement, first.proved.proof),
          "primary verifier rejected baseline");
    check(verify_reference(first.setup.crs, first.setup.precomp,
                           first.proved.statement, first.proved.proof),
          "reference verifier rejected baseline");
    check(verify_deferred(first.setup.crs, first.setup.precomp,
                          first.proved.statement, first.proved.proof),
          "deferred verifier rejected baseline");
    auto validated = PrevalidateVerificationInputs(
        first.setup.crs, first.setup.precomp,
        first.proved.statement, first.proved.proof);
    check(validated && VerifyOnline(*validated),
          "online verifier rejected baseline");
    const auto proof_wire = serialize_proof(3, first.proved.proof);
    check(VerifySerialized(first.setup.crs, first.setup.precomp,
                           first.proved.statement, proof_wire),
          "serialized verifier rejected baseline");

    check(hex(sha256(serialize_crs(first.setup.crs)))
              == "cd30cf965c496ec316ba2c6bbd327d1957370516048e1e75bbbddd04cb707e8f",
          "CRS golden changed");
    check(hex(sha256(serialize_precomputation(
              first.setup.crs, first.setup.precomp)))
              == "5fc536b0eb19997fabda63f657f51b93fac6202384513885b72326bf1761b65f",
          "precomputation golden changed");
    check(hex(sha256(serialize_statement(
              first.setup.crs, first.proved.statement)))
              == "c00c9132f78c9eb53a06f0b4676d36fe8f4f419dc8d9d0b3b34c4126b88d8028",
          "statement golden changed");
    check(hex(sha256(proof_wire))
              == "dd4daa3dbfdf7f6a4170b86b24d2a96a594a2141b87a02b5e85300d0632173ac",
          "proof golden changed");
    check(hex(first.proved.phase2.final_transcript_digest)
              == "bbafda5609bf8c54391cc0452662d5d3b5325f3b342dc36d40f57ad8232c1594",
          "transcript golden changed");
}

void test_size_matrix() {
    constexpr std::size_t proof_sizes[]{2832, 7056, 11280, 15504};
    constexpr std::size_t crs_sizes[]{272, 464, 848, 1616};
    for (std::size_t d = 1; d <= 4; ++d) {
        const auto value = fixture(d);
        check(serialize_proof(d, value.proved.proof).size()
                  == proof_sizes[d - 1], "proof size changed");
        check(serialize_crs(value.setup.crs).size()
                  == crs_sizes[d - 1], "CRS size changed");
        check(Verify(value.setup.crs, value.setup.precomp,
                     value.proved.statement, value.proved.proof),
              "dimension baseline rejected");
    }
}

void test_object_validation() {
    const auto value = fixture(2);
    check(validate_crs(value.setup.crs), "valid CRS was rejected");
    check(validate_statement_input(value.setup.crs,
                                   value.setup.statement_input),
          "valid statement input was rejected");
    check(audit_precomputation(value.setup.crs, value.setup.precomp),
          "valid precomputation was rejected");
    check(validate_statement_shape(value.setup.crs, value.proved.statement)
              && validate_statement_elements(value.proved.statement)
              && validate_statement_digest(value.setup.crs,
                                             value.proved.statement),
          "valid statement was rejected");
    check(validate_proof_shape(value.setup.crs, value.proved.proof)
              && validate_proof_elements(value.proved.proof),
          "valid proof was rejected");

    auto crs = value.setup.crs;
    crs.digest[0] ^= 1;
    check(!validate_crs(crs), "corrupt CRS digest was accepted");
    auto precomputation = value.setup.precomp;
    GT::mul(precomputation.pairing_x[0],
            precomputation.pairing_x[0],
            value.setup.precomp.pairing_x[1]);
    check(validate_precomputation_elements(precomputation)
              && !audit_precomputation(value.setup.crs, precomputation),
          "inconsistent precomputation was accepted");
    precomputation.digest = compute_precomputation_digest(
        value.setup.crs, precomputation);
    check_throws([&] { (void)serialize_precomputation(
        value.setup.crs, precomputation); },
        "inconsistent precomputation was serialized");
    auto statement = value.proved.statement;
    statement.digest[0] ^= 1;
    check(!validate_statement_digest(value.setup.crs, statement),
          "corrupt statement digest was accepted");
    auto proof = value.proved.proof;
    proof.batch_U.pop_back();
    check(!validate_proof_shape(value.setup.crs, proof),
          "malformed proof was accepted");
    DeterministicRng oversized_rng("VPIPBF/TEST/OVERSIZED");
    check_throws([&] { (void)setup_parameters(
        kMaxSerializedDimension + 1, {}, {}, oversized_rng); },
        "setup accepted a dimension above the protocol limit");
}

void test_wire_rejection() {
    const auto value = fixture(2);
    DecodeError error = DecodeError::None;
    auto wire = serialize_proof(2, value.proved.proof);
    VpipBfProof output = value.proved.proof;
    auto truncated = wire;
    truncated.pop_back();
    check(!deserialize_proof(truncated, value.setup.crs, output, &error)
              && error == DecodeError::Truncated
              && serialize_proof(2, output) == wire,
          "truncated proof was accepted or modified output");
    auto trailing = wire;
    trailing.push_back(0);
    check(!deserialize_proof(trailing, value.setup.crs, output, &error)
              && error == DecodeError::TrailingBytes,
          "trailing proof bytes were accepted");
    auto bad_header = wire;
    bad_header[0] ^= 1;
    check(!deserialize_proof(bad_header, value.setup.crs, output, &error)
              && error == DecodeError::WrongMagic,
          "wrong proof magic was accepted");
    auto wrong_version = wire;
    wrong_version[9] ^= 1;
    check(!deserialize_proof(
              wrong_version, value.setup.crs, output, &error)
              && error == DecodeError::UnsupportedVersion,
          "unsupported proof version was accepted");
    auto wrong_curve = wire;
    wrong_curve[11] ^= 1;
    check(!deserialize_proof(
              wrong_curve, value.setup.crs, output, &error)
              && error == DecodeError::WrongCurve,
          "wrong proof curve was accepted");
    auto wrong_dimension = wire;
    wrong_dimension[15] = 3;
    check(!deserialize_proof(
              wrong_dimension, value.setup.crs, output, &error)
              && error == DecodeError::InvalidDimension,
          "wrong proof dimension was accepted");
    auto bad_crs = value.setup.crs;
    bad_crs.digest[0] ^= 1;
    VpipBfPrecomputation decoded_precomputation = value.setup.precomp;
    const auto pre_wire = serialize_precomputation(
        value.setup.crs, value.setup.precomp);
    check(!deserialize_precomputation(
              pre_wire, bad_crs, decoded_precomputation, &error)
              && error == DecodeError::InvalidCrs
              && decoded_precomputation.pairing_x
                  == value.setup.precomp.pairing_x,
          "precomputation decoder accepted invalid CRS or modified output");
    VpipBfStatement decoded_statement = value.proved.statement;
    const auto statement_wire = serialize_statement(
        value.setup.crs, value.proved.statement);
    check(!deserialize_statement(
              statement_wire, bad_crs, decoded_statement, &error)
              && error == DecodeError::InvalidCrs
              && decoded_statement.digest == value.proved.statement.digest,
          "statement decoder accepted invalid CRS or modified output");
    check(!deserialize_proof(wire, bad_crs, output, &error)
              && error == DecodeError::InvalidCrs,
          "proof decoder accepted invalid CRS");
    check_throws([&] { (void)serialize_crs(bad_crs); },
                 "invalid CRS was serialized");
    auto bad_proof = value.proved.proof;
    bad_proof.batch_U.pop_back();
    check_throws([&] { (void)serialize_proof(2, bad_proof); },
                 "malformed proof was serialized");
}

void test_owning_validation_boundary() {
    auto value = fixture(2);
    auto validated = PrevalidateVerificationInputs(
        value.setup.crs, value.setup.precomp,
        value.proved.statement, value.proved.proof);
    check(validated.has_value(), "validation token was not created");
    value.setup.crs.digest[0] ^= 1;
    value.setup.precomp.pairing_x.clear();
    value.proved.statement.X.clear();
    value.proved.proof.batch_U.clear();
    check(VerifyOnline(*validated),
          "validated token changed after source mutation");
}

void test_phase_boundaries() {
    const auto value = fixture(2);
    auto bad_precomputation = value.setup.precomp;
    GT::mul(bad_precomputation.pairing_x[0],
            bad_precomputation.pairing_x[0],
            value.setup.precomp.pairing_x[1]);
    check_throws([&] { (void)prove_phase1(
        value.setup.crs, bad_precomputation, value.setup.statement_input); },
        "phase 1 accepted inconsistent precomputation");
    auto bad_phase1 = value.proved.phase1;
    bad_phase1.transcript_start[0] ^= 1;
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted wrong transcript start");
    bad_phase1 = value.proved.phase1;
    bad_phase1.fresh[1].Phi.pop_back();
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted malformed fresh instance");
    bad_phase1 = value.proved.phase1;
    bad_phase1.fresh[1].D0.clear();
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted invalid fresh GT target");
    bad_phase1 = value.proved.phase1;
    bad_phase1.dynamic_claims[0].E.clear();
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted invalid dynamic claim");
    bad_phase1 = value.proved.phase1;
    GT::mul(bad_phase1.dynamic_claims[0].E,
            bad_phase1.dynamic_claims[0].E,
            value.setup.precomp.pairing_x[0]);
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted inconsistent dynamic claim");
    bad_phase1 = value.proved.phase1;
    G1::add(bad_phase1.fresh[1].Phi[0],
            bad_phase1.fresh[1].Phi[0], value.setup.crs.G[0]);
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted inconsistent fresh witness");
    bad_phase1 = value.proved.phase1;
    bad_phase1.rho[0] = bad_phase1.rho[1];
    check_throws([&] { (void)prove_phase2(
        value.setup.crs, value.setup.precomp, bad_phase1); },
        "phase 2 accepted mismatched challenges");
}

void test_verifier_equivalence_and_tampering() {
    const auto value = fixture(3);
    const auto reference = verify_reference_diagnostic(
        value.setup.crs, value.setup.precomp,
        value.proved.statement, value.proved.proof);
    const auto combined = verify_deferred_combined_with_trace(
        value.setup.crs, value.setup.precomp,
        value.proved.statement, value.proved.proof);
    const auto validated = PrevalidateVerificationInputs(
        value.setup.crs, value.setup.precomp,
        value.proved.statement, value.proved.proof);
    check(validated.has_value(), "equivalence token was not created");
    const auto online = verify_online_with_trace(*validated);
    const auto deferred = verify_deferred_with_trace(
        value.setup.crs, value.setup.precomp,
        value.proved.statement, value.proved.proof);
    check(reference.accepted && combined.accepted
              && online.accepted
              && reference.rho == combined.rho
              && reference.rho == online.rho
              && reference.beta == combined.beta
              && reference.beta == online.beta
              && reference.alpha == combined.alpha
              && reference.alpha == online.alpha
              && reference.gamma == combined.gamma
              && reference.gamma == online.gamma
              && reference.epsilon == combined.epsilon
              && reference.epsilon == online.epsilon
              && reference.dory_residual
                  == combined.evaluated_dory_residual
              && reference.rexp_residual
                  == combined.evaluated_rexp_residual,
          "verification engines disagree");
    check(deferred.accepted
              && deferred.dory_residual_reference
                  == deferred.dory_residual_pippenger
              && deferred.rexp_residual_reference
                  == deferred.rexp_residual_pippenger
              && combined.evaluated_combined_residual
                  == combined.evaluated_combined_reference,
          "independent deferred residuals disagree");

    auto proof = value.proved.proof;
    GT::mul(proof.batch_U[0], proof.batch_U[0],
            value.setup.precomp.pairing_x[0]);
    check(!Verify(value.setup.crs, value.setup.precomp,
                  value.proved.statement, proof)
              && !verify_reference(value.setup.crs, value.setup.precomp,
                                   value.proved.statement, proof),
          "tampered proof was accepted");

    std::vector<GT> bases{value.setup.precomp.pairing_x[0],
                          value.setup.precomp.pairing_x[1],
                          value.setup.precomp.delta1R[0]};
    std::vector<Fr> scalars(3);
    scalars[0] = 2; scalars[1] = 3; scalars[2] = 5;
    check(gt_multiexp_reference(bases, scalars)
              == gt_multiexp_pippenger(bases, scalars),
          "GT multiexponentiation engines disagree");
    const auto inverses = batch_invert_nonzero(scalars);
    for (std::size_t i = 0; i < scalars.size(); ++i)
        check(inverses[i] == inverse_nonzero(scalars[i]),
              "batch inversion disagrees with reference");
}
} // namespace

int main() {
    try {
        initialize();
        test_baseline();
        test_size_matrix();
        test_object_validation();
        test_wire_rejection();
        test_owning_validation_boundary();
        test_phase_boundaries();
        test_verifier_equivalence_and_tampering();
        std::cout << "All vpipbf baseline tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "vpipbf test failure: " << error.what() << '\n';
        return 1;
    }
}

#include "bp/bp.hpp"
#include "bp/helper_prover.hpp"
#include "bp/helper_verifier.hpp"

#include "hp_protocol_internal.hpp"
#include "hv_internal.hpp"

#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

using namespace bp;

void check(bool condition, std::string_view message) {
  if (!condition) throw std::runtime_error(std::string(message));
}

Scalar scalar(std::uint64_t value) {
  Scalar result;
  result.setStr(std::to_string(value));
  return result;
}

Bytes seed(std::string_view label) {
  return Bytes(label.begin(), label.end());
}

std::vector<Scalar> witness(std::size_t count, std::uint64_t offset) {
  std::vector<Scalar> values(count);
  for (std::size_t i = 0; i < count; ++i) {
    values[i] = scalar(offset + 3 * i + i * i);
  }
  return values;
}

Group add(const Group& left, const Group& right) {
  Group result;
  Group::add(result, left, right);
  return result;
}

void test_base_bp(std::size_t n) {
  const PublicParams parameters = Setup(n, seed("bp-tests/base/v1"));
  const auto x = witness(n, 5);
  const auto y = witness(n, 19);
  const Group statement = commit(parameters, x, y);
  Trace prover_trace;
  const Proof proof = Prove(parameters, statement, x, y, &prover_trace);

  check(prover_trace.challenges.size() == parameters.d,
        "BP prover challenge count is inconsistent");
  for (const bool invariant : prover_trace.round_invariants) {
    check(invariant, "a BP prover round invariant failed");
  }

  Trace verifier_trace;
  check(Verify(parameters, statement, proof, &verifier_trace),
        "valid BP proof was rejected");
  check(prover_trace.challenges == verifier_trace.challenges,
        "BP prover and verifier challenges differ");
  check(statement == commit(parameters, x, y),
        "BP commitment is not deterministic");

  const Bytes encoded = serialize_proof(parameters, proof);
  check(encoded.size() == proof_wire_bytes(n),
        "BP wire-size helper is inconsistent");
  Proof decoded;
  check(deserialize_proof(parameters, encoded, decoded),
        "valid BP proof did not decode");
  check(serialize_proof(parameters, decoded) == encoded,
        "BP proof encoding changed on round trip");
  check(VerifySerialized(parameters, statement, encoded),
        "serialized BP verification failed");

  Proof tampered = proof;
  if (!tampered.rounds.empty()) {
    tampered.rounds[0].A = add(tampered.rounds[0].A, parameters.K);
  } else {
    tampered.x_final = scalar(999);
  }
  check(!Verify(parameters, statement, tampered),
        "tampered BP proof was accepted");
  check(!Verify(parameters, add(statement, parameters.K), proof),
        "tampered BP statement was accepted");

  auto truncated = encoded;
  truncated.pop_back();
  decoded = proof;
  check(!deserialize_proof(parameters, truncated, decoded),
        "truncated BP proof was accepted");
  check(decoded.rounds.empty(),
        "failed BP decode retained stale proof rounds");
  auto trailing = encoded;
  trailing.push_back(0);
  check(!deserialize_proof(parameters, trailing, decoded),
        "BP proof with trailing bytes was accepted");
  auto noncanonical = encoded;
  std::fill(noncanonical.end() - scalar_bytes(), noncanonical.end(), 0xff);
  check(!deserialize_proof(parameters, noncanonical, decoded),
        "non-canonical BP scalar was accepted");
}

void test_helper_prover() {
  constexpr std::size_t n = 2;
  const HpSetupResult setup = setup_hp(n, seed("bp-tests/hp/v1"));
  const auto x = witness(n, 7);
  const auto y = witness(n, 31);
  const Group statement = commit(setup.pp.bp, x, y);
  const HpProof proof =
      prove_hp(setup.pp, setup.helper_precomp, statement, x, y);

  const auto result =
      verify_hp(setup.pp, setup.client_precomp, statement, x, y, proof);
  check(result.has_value(), "valid HP proof was rejected");
  check(Verify(setup.pp.bp, statement, *result),
        "HP output is not a valid BP proof");
  check(bp::hp_internal::verify_hp_reference_for_test(
            setup.pp, setup.client_precomp, statement, x, y, proof),
        "HP reference verifier rejected a valid proof");
  check(!bp::hp_internal::verify_hp_wrong_round_encoding_for_test(
             setup.pp, setup.client_precomp, statement, x, y, proof),
        "HP verifier accepted the wrong round-index encoding");

  const auto trace = bp::hp_internal::reconstruct_hp_trace(
      setup.pp, statement, x, y, proof.rounds);
  check(bp::hp_internal::direct_hp_aggregate_msm_for_test(
            setup.pp, trace.aggregate) == trace.aggregate.X_gamma,
        "HP aggregate does not match direct MSM");

  const Bytes encoded = serialize_hp_proof(setup.pp, proof);
  check(encoded.size() == hp_proof_wire_bytes(n),
        "HP wire-size helper is inconsistent");
  HpProof decoded;
  check(deserialize_hp_proof(setup.pp, encoded, decoded),
        "valid HP proof did not decode");
  check(serialize_hp_proof(setup.pp, decoded) == encoded,
        "HP proof encoding changed on round trip");
  check(verify_hp_serialized(
            setup.pp, setup.client_precomp, statement, x, y, encoded)
            .has_value(),
        "serialized HP verification failed");
  auto truncated_wire = encoded;
  truncated_wire.pop_back();
  decoded = proof;
  check(!deserialize_hp_proof(setup.pp, truncated_wire, decoded),
        "truncated HP proof was accepted");
  check(decoded.rounds.empty() && !decoded.vme_proof,
        "failed HP decode retained stale proof state");
  auto trailing_wire = encoded;
  trailing_wire.push_back(0);
  check(!deserialize_hp_proof(setup.pp, trailing_wire, decoded),
        "HP proof with trailing bytes was accepted");

  HpProof tampered = proof;
  tampered.rounds[0].B = add(tampered.rounds[0].B, setup.pp.bp.K);
  check(!verify_hp(
             setup.pp, setup.client_precomp, statement, x, y, tampered),
        "tampered HP rounds were accepted");
  tampered = proof;
  tampered.vme_proof->theta_final =
      add(tampered.vme_proof->theta_final, setup.pp.bp.K);
  check(!verify_hp(
             setup.pp, setup.client_precomp, statement, x, y, tampered),
        "tampered HP VME proof was accepted");
  check(!bp::hp_internal::verify_hp_reference_for_test(
             setup.pp, setup.client_precomp, statement, x, y, tampered),
        "HP reference verifier accepted a tampered VME proof");
}

void test_helper_verifier() {
  constexpr std::size_t n = 2;
  const HvSetupResult setup =
      setup_helper_verifier(n, seed("bp-tests/hv/v1"));
  const auto x = witness(n, 11);
  const auto y = witness(n, 43);
  HvStatement statement;
  statement.Z = commit(setup.pp.bp, x, y);
  statement.bulletproof =
      Prove(setup.pp.bp, statement.Z, x, y);

  const auto proof = prove_helper_verifier(
      setup.pp, setup.helper_precomp, statement);
  check(proof.has_value(), "HV prover failed for a valid statement");
  check(verify_helper_verifier(
            setup.pp, setup.verifier_precomp, statement, *proof),
        "valid HV proof was rejected");
  check(bp::hv_internal::verify_helper_verifier_reference_for_test(
            setup.pp, setup.verifier_precomp, statement, *proof),
        "HV reference verifier rejected a valid proof");

  bp::hv_internal::PreparedStatement prepared;
  check(bp::hv_internal::prepare_statement_for_test(
            setup.pp, statement, prepared),
        "HV statement preparation failed");
  check(bp::hv_internal::direct_vme_relation_for_test(
            setup.pp, prepared) == prepared.X0,
        "HV prepared relation does not match direct MSM");

  HvPreparedStatementCache cache;
  check(verify_helper_verifier_cached(
            setup.pp, setup.verifier_precomp, statement, *proof, cache),
        "cached HV verification failed on first use");
  check(cache.ready(), "HV cache was not populated");
  check(verify_helper_verifier_cached(
            setup.pp, setup.verifier_precomp, statement, *proof, cache),
        "cached HV verification failed on reuse");
  cache.clear();
  check(!cache.ready(), "HV cache clear did not reset readiness");
  check(verify_helper_verifier_cached(
            setup.pp, setup.verifier_precomp, statement, *proof, cache),
        "cached HV verification failed after clear");

  const Bytes encoded_proof = serialize_hv_proof(setup.pp, *proof);
  const Bytes encoded_statement = serialize_hv_statement(setup.pp, statement);
  check(encoded_proof.size() == hv_proof_wire_bytes(n),
        "HV proof-size helper is inconsistent");
  check(encoded_statement.size() == hv_statement_wire_bytes(n),
        "HV statement-size helper is inconsistent");
  HvProof decoded_proof;
  HvStatement decoded_statement;
  check(deserialize_hv_proof(setup.pp, encoded_proof, decoded_proof),
        "valid HV proof did not decode");
  check(deserialize_hv_statement(
            setup.pp, encoded_statement, decoded_statement),
        "valid HV statement did not decode");
  check(verify_helper_verifier_serialized(
            setup.pp, setup.verifier_precomp, statement, encoded_proof),
        "serialized HV proof verification failed");
  check(verify_helper_verifier_serialized(
            setup.pp, setup.verifier_precomp, encoded_statement,
            encoded_proof),
        "serialized HV statement/proof verification failed");
  auto truncated_proof = encoded_proof;
  truncated_proof.pop_back();
  decoded_proof = *proof;
  check(!deserialize_hv_proof(setup.pp, truncated_proof, decoded_proof),
        "truncated HV proof was accepted");
  check(decoded_proof.vme_proof.dory_folds.empty(),
        "failed HV proof decode retained stale state");
  auto truncated_statement = encoded_statement;
  truncated_statement.pop_back();
  decoded_statement = statement;
  check(!deserialize_hv_statement(
             setup.pp, truncated_statement, decoded_statement),
        "truncated HV statement was accepted");
  check(decoded_statement.bulletproof.rounds.empty(),
        "failed HV statement decode retained stale state");

  HvProof tampered = *proof;
  tampered.vme_proof.theta_final =
      add(tampered.vme_proof.theta_final, setup.pp.bp.K);
  check(!verify_helper_verifier(
             setup.pp, setup.verifier_precomp, statement, tampered),
        "tampered HV proof was accepted");
  HvStatement tampered_statement = statement;
  tampered_statement.Z = add(tampered_statement.Z, setup.pp.bp.K);
  check(!verify_helper_verifier(
             setup.pp, setup.verifier_precomp, tampered_statement, *proof),
        "tampered HV statement was accepted");
  check(!verify_helper_verifier_cached(
             setup.pp, setup.verifier_precomp, tampered_statement, *proof,
             cache),
        "HV cache was reused for a different statement");
}

}  // namespace

int main() {
  try {
    initialize();
    test_base_bp(1);
    test_base_bp(2);
    test_base_bp(4);
    test_helper_prover();
    test_helper_verifier();
    std::cout << "All BP/HP/HV tests passed.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Test failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

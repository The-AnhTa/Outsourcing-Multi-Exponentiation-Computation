#include "blsagg/protocol.hpp"
#include "blsagg/serialization.hpp"
#include "blsagg/transcript.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace blsagg;

Statement make_statement(const PublicParameters& parameters);

void check(bool condition, const char* message) {
  if (!condition) throw std::runtime_error(message);
}

G1 add(const G1& left, const G1& right) {
  G1 out;
  G1::add(out, left, right);
  return out;
}

bool target_equal(const DoryTarget& left, const DoryTarget& right) {
  return left.D0 == right.D0 && left.D1 == right.D1 && left.D2 == right.D2;
}

std::string hex(std::span<const std::uint8_t> bytes) {
  static constexpr char digits[] = "0123456789abcdef";
  std::string out;
  out.reserve(bytes.size() * 2);
  for (const auto byte : bytes) {
    out.push_back(digits[byte >> 4]);
    out.push_back(digits[byte & 0x0f]);
  }
  return out;
}

void test_transcript_golden() {
  const auto setup_result = setup(
      1, AggregationMode::BasicDistinct, "bls-tests/transcript-golden");
  const auto statement = make_statement(setup_result.pp);
  const auto points = hash_messages(setup_result.pp, statement);
  Transcript transcript(setup_result.pp, statement, points);
  const std::array<Bytes, 2> fields{Bytes{0, 1, 2, 3},
                                    Bytes{'g', 'o', 'l', 'd', 'e', 'n'}};
  transcript.absorb("bls-tests/golden-absorb", 7, fields);
  const auto challenge = transcript.challenge_nonzero(
      "bls-tests/golden-challenge", 11);
  const auto actual = hex(serialize(challenge));
  const std::string expected =
      "14045226b1b7a6aece2a67860efbc510090797d90bc0c41ac027f7ca995fb907";
  check(actual == expected, actual.c_str());
}

Statement make_statement(const PublicParameters& parameters) {
  Statement statement;
  std::vector<Fr> secrets;
  secrets.reserve(parameters.k);
  for (std::size_t i = 0; i < parameters.k; ++i) {
    const std::string text = "bls-tests/message/" + std::to_string(i);
    statement.messages.emplace_back(text.begin(), text.end());
    Fr secret;
    secret = static_cast<int>(101 + i);
    secrets.push_back(secret);
    G2 public_key;
    G2::mul(public_key, parameters.H, secret);
    statement.public_keys.push_back(public_key);
  }
  const auto points = hash_messages(parameters, statement);
  statement.sigma_agg.clear();
  for (std::size_t i = 0; i < parameters.k; ++i) {
    G1 signature;
    G1::mul(signature, points[i], secrets[i]);
    G1::add(statement.sigma_agg, statement.sigma_agg, signature);
  }
  return statement;
}

void check_trace_agreement(const VerificationTrace& reference,
                           const VerificationTrace& candidate) {
  check(reference.accepted == candidate.accepted,
        "verifier strategies disagree on acceptance");
  check(reference.g1_rexp_challenges == candidate.g1_rexp_challenges &&
            reference.g2_rexp_challenges == candidate.g2_rexp_challenges &&
            reference.dory_beta == candidate.dory_beta &&
            reference.dory_alpha == candidate.dory_alpha &&
            reference.insert_g1_gamma == candidate.insert_g1_gamma &&
            reference.insert_g2_gamma == candidate.insert_g2_gamma &&
            reference.eta == candidate.eta && reference.zeta == candidate.zeta,
        "verifier strategies derived different challenges");
  check(target_equal(reference.final_dory, candidate.final_dory) &&
            target_equal(reference.final_g1_rexp,
                         candidate.final_g1_rexp) &&
            target_equal(reference.final_g2_rexp,
                         candidate.final_g2_rexp),
        "verifier strategies derived different terminal targets");
}

void test_serialization(const SetupResult& setup_result,
                        const Statement& statement, const Proof& proof) {
  DecodeError error = DecodeError::None;

  const Bytes pp_wire = serialize_public_parameters(setup_result.pp);
  PublicParameters decoded_pp;
  check(deserialize_public_parameters(pp_wire, decoded_pp, &error),
        "public parameters did not decode");
  check(serialize_public_parameters(decoded_pp) == pp_wire,
        "public-parameter round trip changed bytes");

  const Bytes aux_wire =
      serialize_precomputation(setup_result.pp, setup_result.aux);
  Precomputation decoded_aux;
  check(deserialize_precomputation(
            aux_wire, setup_result.pp, decoded_aux, &error),
        "precomputation did not decode");
  check(serialize_precomputation(setup_result.pp, decoded_aux) == aux_wire,
        "precomputation round trip changed bytes");

  const Bytes statement_wire = serialize_statement(setup_result.pp, statement);
  Statement decoded_statement;
  check(deserialize_statement(
            statement_wire, setup_result.pp, decoded_statement, &error),
        "statement did not decode");
  check(serialize_statement(setup_result.pp, decoded_statement) ==
            statement_wire,
        "statement round trip changed bytes");

  const Bytes proof_wire = serialize_proof(setup_result.pp, proof);
  Proof decoded_proof;
  check(deserialize_proof(proof_wire, setup_result.pp, decoded_proof, &error),
        "proof did not decode");
  check(serialize_proof(setup_result.pp, decoded_proof) == proof_wire,
        "proof round trip changed bytes");
  check(proof_wire.size() == proof_wire_bytes(setup_result.pp, proof),
        "proof wire-size helper is inconsistent");
  check(proof_mathematical_payload_bytes(proof) == proof_payload_bytes(proof),
        "proof payload-size helpers disagree");

  auto truncated = proof_wire;
  truncated.pop_back();
  decoded_proof = proof;
  check(!deserialize_proof(truncated, setup_result.pp, decoded_proof, &error) &&
            error == DecodeError::Truncated,
        "truncated proof was not rejected precisely");
  check(decoded_proof.g1_rexp_claims.empty() &&
            decoded_proof.g2_rexp_claims.empty() &&
            decoded_proof.dory_steps.empty() &&
            decoded_proof.insert_g1_u.empty() &&
            decoded_proof.insert_g2_u.empty(),
        "failed proof decode retained stale output");
  auto trailing = proof_wire;
  trailing.push_back(0);
  check(!deserialize_proof(trailing, setup_result.pp, decoded_proof, &error) &&
            error == DecodeError::TrailingBytes,
        "proof trailing bytes were not rejected precisely");
  auto wrong_magic = proof_wire;
  wrong_magic[0] ^= 1;
  check(!deserialize_proof(
             wrong_magic, setup_result.pp, decoded_proof, &error) &&
            error == DecodeError::WrongMagic,
        "wrong proof magic was not rejected precisely");

  auto invalid_count = proof_wire;
  std::size_t first_count_offset = 24;
  for (const auto* value : {&proof.cm_M, &proof.cm_pk, &proof.T})
    first_count_offset += 8 + serialize(*value).size();
  std::fill_n(invalid_count.begin() + first_count_offset, 8, 0xff);
  check(!deserialize_proof(
             invalid_count, setup_result.pp, decoded_proof, &error) &&
            error == DecodeError::InvalidLength,
        "invalid proof count did not preserve its precise decode error");

  auto bad_pp_wire = pp_wire;
  bad_pp_wire.back() ^= 1;
  decoded_pp = setup_result.pp;
  check(!deserialize_public_parameters(bad_pp_wire, decoded_pp, &error) &&
            error == DecodeError::InvalidDigest && decoded_pp.k == 0 &&
            decoded_pp.Gamma.empty() && decoded_pp.Lambda.empty(),
        "invalid parameter digest was not failure-atomic");

  auto bad_statement_wire = statement_wire;
  bad_statement_wire.pop_back();
  decoded_statement = statement;
  check(!deserialize_statement(
             bad_statement_wire, setup_result.pp, decoded_statement, &error) &&
            error == DecodeError::Truncated &&
            decoded_statement.messages.empty() &&
            decoded_statement.public_keys.empty(),
        "failed statement decode retained stale output");
}

void test_protocol(std::size_t d, AggregationMode mode) {
  const std::string mode_name =
      mode == AggregationMode::BasicDistinct ? "basic" : "augmented";
  const SetupResult setup_result =
      setup(d, mode, "bls-tests/" + mode_name + "/d=" + std::to_string(d));
  const Statement statement = make_statement(setup_result.pp);
  const Proof proof = prove(
      setup_result.pp, setup_result.aux, statement);
  const Proof repeated = prove(
      setup_result.pp, setup_result.aux, statement);
  check(serialize_proof(setup_result.pp, proof) ==
            serialize_proof(setup_result.pp, repeated),
        "prover is not deterministic");

  check(direct_bls_verify(setup_result.pp, statement),
        "direct BLS verification failed");
  auto invalid_parameters = setup_result.pp;
  invalid_parameters.digest[0] ^= 1;
  check(!direct_bls_verify(invalid_parameters, statement),
        "direct verification accepted invalid public parameters");
  check(verify_safe(
            setup_result.pp, setup_result.aux, statement, proof),
        "safe proof verification failed");

  const auto context = prepare_verifier_context(
      setup_result.pp, setup_result.aux, statement);
  check(context.has_value(), "verifier context preparation failed");
  const Bytes proof_wire = serialize_proof(setup_result.pp, proof);
  const auto validated = deserialize_and_validate_proof(
      proof_wire, setup_result.pp);
  check(validated.has_value(), "validated proof construction failed");
  check(validated->wire_binding() != Digest{},
        "validated proof did not retain a wire binding");
  check(verify_online(*context, proof),
        "typed online verification failed");
  check(verify_online(*context, *validated),
        "validated online verification failed");
  check(verify_online_sequential_msm(*context, *validated),
        "sequential verifier failed");
  check(verify_online_parallel_msm(*context, *validated),
        "parallel verifier failed");
  check(verify_online_symbolic_gt(*context, *validated),
        "symbolic GT verifier failed");
  check(verify_online_split_g2_msm(*context, *validated),
        "split G2 verifier failed");

  if (mode == AggregationMode::BasicDistinct && d == 1) {
    for (int repetition = 0; repetition < 2; ++repetition) {
      check(verify_online_parallel_msm(*context, *validated),
            "repeated parallel verification was unstable");
      check(verify_online_split_g2_msm(*context, *validated),
            "repeated split-G2 verification was unstable");
    }
  }

  if (mode == AggregationMode::BasicDistinct && d == 1) {
    const auto other_setup = setup(
        d, mode, "bls-tests/wrong-validated-proof-context");
    const auto other_context = prepare_verifier_context(
        other_setup.pp, other_setup.aux, statement);
    check(other_context.has_value() &&
              !verify_online(*other_context, *validated),
          "validated proof binding was accepted by different parameters");
  }

  const VerificationTrace reference =
      verify_online_sequential_msm_diagnostic(*context, *validated);
  check_trace_agreement(
      reference, verify_online_parallel_msm_diagnostic(*context, *validated));
  check_trace_agreement(
      reference,
      verify_online_symbolic_gt_differential_trace(*context, *validated));

  Proof tampered = proof;
  tampered.Phi_final = add(tampered.Phi_final, setup_result.pp.L);
  check(!verify_safe(
             setup_result.pp, setup_result.aux, statement, tampered),
        "tampered proof was accepted");
  check(!verify_online_diagnostic(*context, tampered).accepted,
        "diagnostic verifier accepted a tampered proof");

  test_serialization(setup_result, statement, proof);

  if (mode == AggregationMode::BasicDistinct) {
    Statement duplicate = statement;
    duplicate.messages[1] = duplicate.messages[0];
    check(!direct_bls_verify(setup_result.pp, duplicate),
          "duplicate Basic messages were accepted directly");
    check(!prepare_verifier_context(
               setup_result.pp, setup_result.aux, duplicate),
          "duplicate Basic messages produced a verifier context");
  }

  if (mode == AggregationMode::BasicDistinct && d == 1) {
    check(proof_wire.size() == 5384,
          "d=1 proof wire size changed from baseline");
    check(serialize_public_parameters(setup_result.pp).size() == 488,
          "d=1 CRS wire size changed from baseline");
  }
  if (mode == AggregationMode::BasicDistinct && d == 2) {
    check(proof_wire.size() == 11656,
          "d=2 proof wire size changed from baseline");
    check(serialize_public_parameters(setup_result.pp).size() == 712,
          "d=2 CRS wire size changed from baseline");
  }
}

}  

int main() {
  try {
    initialize();
    test_transcript_golden();
    test_protocol(1, AggregationMode::BasicDistinct);
    test_protocol(2, AggregationMode::BasicDistinct);
    test_protocol(1, AggregationMode::Augmented);
    test_protocol(2, AggregationMode::Augmented);
    std::cout << "All BLS aggregate tests passed.\n";
    return EXIT_SUCCESS;
  } catch (const std::exception& error) {
    std::cerr << "Test failure: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}

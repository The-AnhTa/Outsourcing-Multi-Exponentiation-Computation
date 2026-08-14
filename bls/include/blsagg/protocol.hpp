#pragma once

#include <mcl/bn.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace blsagg {

class ValidatedProof;

using Fr = mcl::bn::Fr;
using G1 = mcl::bn::G1;
using G2 = mcl::bn::G2;
using GT = mcl::bn::Fp12;
using Bytes = std::vector<std::uint8_t>;
using Digest = std::array<std::uint8_t, 32>;

enum class AggregationMode : std::uint8_t { BasicDistinct = 1, Augmented = 2 };

struct RexpClaims { GT E, F, TL, TR; };
struct DoryStep { GT A1L, A1R, A2L, A2R, W1, W2; };
struct DoryTarget { GT D0, D1, D2; };
struct DoryWitness { std::vector<G1> Phi; std::vector<G2> Theta; };
struct DoryInstance { DoryTarget target; DoryWitness witness; };

struct PublicParameters {
  std::size_t d{}, k{};
  AggregationMode mode{AggregationMode::BasicDistinct};
  G2 H;
  std::vector<G1> Gamma;
  std::vector<G2> Lambda;
  G1 L;
  G2 Lprime;
  Digest digest{};
};

struct Precomputation {
  std::vector<std::vector<G1>> gamma_chain;
  std::vector<std::vector<G2>> lambda_chain;
  std::vector<GT> X;
  std::vector<GT> delta1L, delta1R, delta2L, delta2R;
  RexpClaims g1_round0, g2_round0;
  Digest digest{};
};

struct Statement {
  G1 sigma_agg;
  std::vector<Bytes> messages;
  std::vector<G2> public_keys;
};

struct Proof {
  GT cm_M, cm_pk, T;
  std::vector<RexpClaims> g1_rexp_claims;
  std::vector<RexpClaims> g2_rexp_claims;
  G1 R_Gamma;
  G2 R_Lambda;
  GT U1, U2;
  std::vector<DoryStep> dory_steps;
  std::vector<GT> insert_g1_u;
  std::vector<GT> insert_g2_u;
  G1 Phi_final;
  G2 Theta_final;
};

struct SetupResult { PublicParameters pp; Precomputation aux; };

class ValidatedVerifierContext {
 public:
  ValidatedVerifierContext() = delete;
  ValidatedVerifierContext(const ValidatedVerifierContext&) = default;
  ValidatedVerifierContext(ValidatedVerifierContext&&) = default;
  ValidatedVerifierContext& operator=(const ValidatedVerifierContext&) = delete;
  ValidatedVerifierContext& operator=(ValidatedVerifierContext&&) = delete;
  const PublicParameters& parameters() const { return pp_; }
  const Precomputation& precomputation() const { return aux_; }
  const Statement& statement() const { return statement_; }
  std::span<const G1> message_points() const { return message_points_; }
  const Digest& binding() const { return binding_; }
 private:
  ValidatedVerifierContext(PublicParameters, Precomputation, Statement,
                           std::vector<G1>, Digest);
  const PublicParameters pp_;
  const Precomputation aux_;
  const Statement statement_;
  const std::vector<G1> message_points_;
  const Digest binding_;
  friend std::optional<ValidatedVerifierContext> prepare_verifier_context(
      const PublicParameters&, const Precomputation&, const Statement&);
};

struct VerificationTrace {
  bool accepted{};
  G1 Y_M;
  G2 Y_pk;
  DoryTarget final_dory;
  DoryTarget final_g1_rexp;
  DoryTarget final_g2_rexp;
  std::vector<Fr> g1_rexp_challenges;
  std::vector<Fr> g2_rexp_challenges;
  std::vector<Fr> dory_beta;
  std::vector<Fr> dory_alpha;
  std::vector<Fr> insert_g1_gamma;
  std::vector<Fr> insert_g2_gamma;
  Fr eta, zeta;
  std::vector<DoryTarget> g1_rexp_targets;
  std::vector<DoryTarget> g2_rexp_targets;
  std::vector<DoryTarget> dory_fold_targets;
  std::vector<DoryTarget> accumulator_targets;
  double total_online_ms{};
  double proof_validation_ms{};
  double transcript_hashing_ms{};
  double rexp_gt_recurrence_ms{};
  double tensor_reconstruction_ms{};
  double message_g1_msm_ms{};
  double public_key_g2_msm_ms{};
  double public_key_g2_left_msm_ms{};
  double public_key_g2_right_msm_ms{};
  double public_input_msm_wall_ms{};
  double application_batching_ms{};
  double integrated_dory_ms{};
  double gt_multiexponentiation_ms{};
  double final_challenges_ms{};
  double terminal_checks_ms{};
};

void initialize();
SetupResult setup(std::size_t d, AggregationMode mode, std::string_view crs_seed);
Precomputation precompute(const PublicParameters&);
std::vector<G1> hash_messages(const PublicParameters&, const Statement&);
Proof prove(const PublicParameters&, const Precomputation&, const Statement&);
std::optional<ValidatedVerifierContext> prepare_verifier_context(
    const PublicParameters&, const Precomputation&, const Statement&);
bool verify_safe(const PublicParameters&, const Precomputation&, const Statement&, const Proof&);
bool verify_online(const ValidatedVerifierContext&, const Proof&);
VerificationTrace verify_online_diagnostic(const ValidatedVerifierContext&, const Proof&);
bool verify_online(const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_diagnostic(const ValidatedVerifierContext&, const ValidatedProof&);
bool verify_online_sequential_msm(const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_sequential_msm_diagnostic(
    const ValidatedVerifierContext&, const ValidatedProof&);
bool verify_online_parallel_msm(const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_parallel_msm_diagnostic(
    const ValidatedVerifierContext&, const ValidatedProof&);
bool verify_online_symbolic_gt(const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_symbolic_gt_diagnostic(
    const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_symbolic_gt_differential_trace(
    const ValidatedVerifierContext&, const ValidatedProof&);
bool verify_online_split_g2_msm(const ValidatedVerifierContext&, const ValidatedProof&);
VerificationTrace verify_online_split_g2_msm_diagnostic(
    const ValidatedVerifierContext&, const ValidatedProof&);

GT direct_pairing_product(std::span<const G1>, std::span<const G2>);
bool direct_bls_verify(const PublicParameters&, const Statement&);
std::size_t proof_payload_bytes(const Proof&);
std::vector<Fr> tensor_vector(std::span<const Fr>);

}

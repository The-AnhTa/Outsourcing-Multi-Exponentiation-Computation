#pragma once
#include <mcl/bn.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace vpip_bf {
using Fr = mcl::bn::Fr;
using G1 = mcl::bn::G1;
using G2 = mcl::bn::G2;
using GT = mcl::bn::Fp12;
using Digest = std::array<std::uint8_t, 32>;
using Bytes = std::vector<std::uint8_t>;

struct VpipBfCRS { std::size_t d{}, n{}; std::vector<G1> G; std::vector<G2> H; G2 Lprime; Digest digest{}; };
struct VpipBfPrecomputation { std::vector<GT> pairing_x, delta1R, delta2R; Digest digest{}; };
struct VpipBfStatementInput { std::vector<G1> X; Digest digest{}; };
struct VpipBfProverInput { std::vector<G1> X; };
struct VpipBfStatement { std::vector<G1> X; GT C; Digest digest{}; };
struct FreshDoryInstance { GT D0, D1, D2; std::vector<G1> Phi; std::vector<G2> Theta; };
struct RexpClaims { GT E, F, TL, TR; };
struct SetupResult { VpipBfCRS crs; VpipBfPrecomputation precomp; VpipBfStatementInput statement_input; VpipBfProverInput prover_input; };
struct Phase1Result {
  VpipBfStatement statement;
  std::vector<Fr> rho, r;
  std::vector<RexpClaims> dynamic_claims;
  std::vector<FreshDoryInstance> fresh;
  G1 R;
  Digest transcript_start{};
  Digest transcript_after_R{};
};
struct DoryTargetState { GT D0, D1, D2; };
struct DoryWitnessState { std::vector<G1> Phi; std::vector<G2> Theta; };
struct DoryInstanceState { DoryTargetState target; DoryWitnessState witness; };
struct DoryFoldProof { GT D1L, D1R, D2L, D2R, W1, W2; };
struct Phase2ProofData { std::vector<DoryFoldProof> dory_folds; std::vector<GT> batch_U; G1 PhiFinal; G2 ThetaFinal; };
struct Phase2ChallengeTrace { std::vector<Fr> beta, alpha, gamma; Fr epsilon; };
struct Phase2Result {
  Phase2ProofData proof;
  Phase2ChallengeTrace challenges;
  DoryTargetState final_aggregate_target;
  G1 Y;
  Digest final_transcript_digest{};

  DoryInstanceState initial_instance;
  std::vector<DoryInstanceState> folded;
  std::vector<DoryInstanceState> aggregate;
};

using Crs = VpipBfCRS;
using Precomputation = VpipBfPrecomputation;
using Statement = VpipBfStatement;
using Witness = VpipBfProverInput;
}

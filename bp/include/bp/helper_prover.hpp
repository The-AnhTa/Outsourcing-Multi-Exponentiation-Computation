#pragma once

#include "bp/bp.hpp"
#include <optional>

namespace bp {

using G1 = mcl::bn::G1;
using GT = mcl::bn::Fp12;

struct VmeRexpClaims { GT E, F, TL, TR; };
struct VmeDoryFold { GT D1L, D1R, D2L, D2R, W1, W2; };

struct VmeProof {
  std::vector<VmeRexpClaims> rexp_claims;
  G1 R;
  std::vector<VmeDoryFold> dory_folds;
  std::vector<GT> batch_U;
  G1 phi_final;
  Group theta_final;
};

struct VmePublicParams {
  std::size_t dimension{};
  std::size_t log_dimension{};
  std::vector<G1> auxiliary_G;
  std::vector<Group> fixed_P;
  G1 L;
  Group Lprime;
  std::string transcript_domain{"BPVME/HP/VME/G2/v1"};
  Digest digest{};
};

struct VmePrecomputation {
  std::vector<GT> pairing_x;
  std::vector<GT> delta1R;
  std::vector<GT> delta2R;
  GT pairing_LLprime;
  Digest binding_digest{};
};

struct HpPublicParams {
  PublicParams bp;
  std::size_t N{};
  VmePublicParams vme;
  Digest crs_digest{};
};

struct HpHelperPrecomp { VmePrecomputation vme; Digest crs_digest{}; };
struct HpClientPrecomp { VmePrecomputation vme; Digest crs_digest{}; };

struct HpSetupResult {
  HpPublicParams pp;
  HpHelperPrecomp helper_precomp;
  HpClientPrecomp client_precomp;
};

struct HpProof {
  std::vector<RoundMessage> rounds;
  std::optional<VmeProof> vme_proof;
};

struct HpInstance {
  Group Z;
  std::vector<Scalar> x;
  std::vector<Scalar> y;
};

struct HpSetupTimings {
  double setup_ms{};
  double precomputation_ms{};
};

struct HpProveTimings {
  double bp_round_generation_ms{};
  double outer_aggregation_ms{};
  double vme_prove_ms{};
  double total_ms{};
};

struct HpVerifyTimings {
  double parse_and_validation_ms{};
  double bp_transcript_replay_ms{};
  double outer_aggregation_ms{};
  double vme_verify_ms{};
  double vme_transcript_ms{};
  double vme_batch_inversion_ms{};
  double vme_symbolic_recurrence_ms{};
  double vme_terminal_assembly_ms{};
  double vme_gt_multiexp_ms{};
  double vme_multi_pairing_ms{};
  double total_ms{};
};

HpSetupResult setup_hp(std::size_t n, std::span<const std::uint8_t> setup_seed,
                       HpSetupTimings* timings = nullptr);
HpProof prove_hp(const HpPublicParams&, const HpHelperPrecomp&, const Group& Z,
                 std::span<const Scalar> x, std::span<const Scalar> y,
                 HpProveTimings* timings = nullptr);
std::optional<Proof> verify_hp(const HpPublicParams&, const HpClientPrecomp&,
                               const Group& Z, std::span<const Scalar> x,
                               std::span<const Scalar> y, const HpProof&,
                               HpVerifyTimings* timings = nullptr);
std::optional<Proof> verify_hp_serialized(
    const HpPublicParams&, const HpClientPrecomp&, const Group& Z,
    std::span<const Scalar> x, std::span<const Scalar> y,
    std::span<const std::uint8_t> helper_proof,
    HpVerifyTimings* timings = nullptr);

HpInstance generate_hp_instance(const HpPublicParams&);

Bytes serialize_hp_proof(const HpPublicParams&, const HpProof&);
bool deserialize_hp_proof(const HpPublicParams&, std::span<const std::uint8_t>,
                          HpProof&) noexcept;
std::size_t vme_proof_payload_bytes(std::size_t n);
std::size_t hp_proof_payload_bytes(std::size_t n);
std::size_t hp_proof_wire_bytes(std::size_t n);

}

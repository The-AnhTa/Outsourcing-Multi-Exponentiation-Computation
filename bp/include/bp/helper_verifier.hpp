#pragma once

#include "bp/helper_prover.hpp"

namespace bp {

inline constexpr const char* kHvContextDomain = "BPVME/HV/CONTEXT/v1";
inline constexpr const char* kHvVmeDomain = "BPVME/HV/VME/G2/v1";

struct HvPublicParams {
  PublicParams bp;
  std::size_t N{};
  VmePublicParams vme;
  Digest bp_parameter_digest{};
  Digest bp_transcript_parameter_digest{};
  Digest P_digest{};
  Digest crs_digest{};
};

struct HvHelperPrecomputation {
  VmePrecomputation vme;
  Digest crs_digest{};
  Digest bp_transcript_parameter_digest{};
};

struct HvVerifierPrecomputation {
  VmePrecomputation vme;
  Digest crs_digest{};
  Digest bp_parameter_digest{};
  Digest bp_transcript_parameter_digest{};
  Digest P_digest{};
};

struct HvSetupResult {
  HvPublicParams pp;
  HvHelperPrecomputation helper_precomp;
  HvVerifierPrecomputation verifier_precomp;
};

struct HvStatement {
  Group Z;
  Proof bulletproof;
};

struct HvProof { VmeProof vme_proof; };

struct HvPreparedStatementCache {
  bool ready{};
  Digest crs_digest{};
  Digest statement_context{};
  std::vector<Scalar> z_v;
  Group X0;
};

struct HvInstance {
  HvStatement statement;
  std::vector<Scalar> x;
  std::vector<Scalar> y;
};

struct HvSetupTimings {
  double setup_ms{};
  double verifier_precomputation_ms{};
};

struct HvProveTimings {
  double statement_preparation_ms{};
  double relation_msm_ms{};
  double vme_prove_ms{};
  double total_ms{};
};

struct HvVerifyTimings {
  double validation_ms{};
  double bp_transcript_ms{};
  double batch_inversion_ms{};
  double weight_generation_ms{};
  double z_vector_ms{};
  double x0_msm_ms{};
  double context_digest_ms{};
  double vme_transcript_ms{};
  double vme_batch_inversion_ms{};
  double vme_symbolic_recurrence_ms{};
  double vme_terminal_assembly_ms{};
  double vme_gt_multiexp_ms{};
  double vme_multi_pairing_ms{};
  double vme_verify_ms{};
  double total_ms{};
};

HvSetupResult setup_helper_verifier(
    std::size_t n, std::span<const std::uint8_t> setup_seed,
    HvSetupTimings* timings = nullptr);
std::optional<HvProof> prove_helper_verifier(
    const HvPublicParams&, const HvHelperPrecomputation&, const HvStatement&,
    HvProveTimings* timings = nullptr);
bool verify_helper_verifier(
    const HvPublicParams&, const HvVerifierPrecomputation&, const HvStatement&,
    const HvProof&, HvVerifyTimings* timings = nullptr);
bool verify_helper_verifier_cached(
    const HvPublicParams&, const HvVerifierPrecomputation&, const HvStatement&,
    const HvProof&, HvPreparedStatementCache&,
    HvVerifyTimings* timings = nullptr);
bool verify_helper_verifier_serialized(
    const HvPublicParams&, const HvVerifierPrecomputation&, const HvStatement&,
    std::span<const std::uint8_t> proof, HvVerifyTimings* timings = nullptr);
bool verify_helper_verifier_serialized(
    const HvPublicParams&, const HvVerifierPrecomputation&,
    std::span<const std::uint8_t> statement,
    std::span<const std::uint8_t> proof, HvVerifyTimings* timings = nullptr);

HvInstance generate_helper_verifier_instance(const HvPublicParams&);

Bytes serialize_hv_proof(const HvPublicParams&, const HvProof&);
bool deserialize_hv_proof(const HvPublicParams&,
                          std::span<const std::uint8_t>, HvProof&) noexcept;
std::size_t hv_proof_payload_bytes(std::size_t n);
std::size_t hv_proof_wire_bytes(std::size_t n);
Bytes serialize_hv_statement(const HvPublicParams&, const HvStatement&);
bool deserialize_hv_statement(const HvPublicParams&,
                              std::span<const std::uint8_t>,
                              HvStatement&) noexcept;
std::size_t hv_statement_wire_bytes(std::size_t n);

}

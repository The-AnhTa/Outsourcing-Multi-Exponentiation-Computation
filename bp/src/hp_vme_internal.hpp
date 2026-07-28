#pragma once

#include "bp/helper_prover.hpp"

namespace bp::hp_internal {

struct VmeStatement {
  Digest application_context{};
  Group X;
  std::vector<Scalar> z;
  Digest digest{};
};

struct VmeVerificationStats {
  bool accepted{};
  double transcript_ms{};
  double batch_inversion_ms{};
  double recurrence_ms{};
  double terminal_ms{};
  double gt_multiexp_ms{};
  double multi_pairing_ms{};
  double total_ms{};
  std::size_t gt_terms{};
  std::size_t pairing_terms{};
};

VmePublicParams setup_vme(std::span<const Group> fixed_P,
                          std::span<const std::uint8_t> seed,
                          std::string_view transcript_domain =
                              "BPVME/HP/VME/G2/v1");
VmePrecomputation precompute_vme(const VmePublicParams&);
VmeProof prove_vme(const VmePublicParams&, const VmePrecomputation&,
                   const Digest& application_context, const Group& X,
                   std::span<const Scalar> z);
bool verify_vme_reference(const VmePublicParams&, const VmePrecomputation&,
                          const Digest& application_context, const Group& X,
                          std::span<const Scalar> z, const VmeProof&,
                          VmeVerificationStats* = nullptr);
bool verify_vme_optimized(const VmePublicParams&, const VmePrecomputation&,
                          const Digest& application_context, const Group& X,
                          std::span<const Scalar> z, const VmeProof&,
                          VmeVerificationStats* = nullptr);
bool verify_vme_prevalidated_optimized(
    const VmePublicParams&, const VmePrecomputation&,
    const Digest& application_context, const Group& X,
    std::span<const Scalar> z, const VmeProof&,
    VmeVerificationStats* = nullptr);
bool validate_vme_params(const VmePublicParams&) noexcept;
bool validate_vme_precomputation_binding(const VmePublicParams&,
                                         const VmePrecomputation&) noexcept;
bool validate_vme_proof(const VmePublicParams&, const VmeProof&) noexcept;
bool validate_vme_proof_batched(const VmePublicParams&,
                                const VmeProof&) noexcept;

Bytes serialize_g1(const G1&);
Bytes serialize_gt(const GT&);
std::size_t g1_bytes();
std::size_t gt_bytes();
bool deserialize_g1(std::span<const std::uint8_t>, G1&);
bool deserialize_gt(std::span<const std::uint8_t>, GT&);
bool deserialize_gt_canonical_unchecked(std::span<const std::uint8_t>, GT&);

}

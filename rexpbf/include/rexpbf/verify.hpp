#pragma once
#include "rexpbf/prove.hpp"
#include <optional>

namespace rexpbf {
struct VerifyDiagnostics {
    std::size_t gt_terms_before_coalescing{};
    std::size_t gt_terms_after_coalescing{};
    std::size_t gt_zero_scalar_terms{};
    std::size_t gt_nonzero_scalar_terms{};
    std::size_t gt_multiexp_calls{};
    std::size_t terminal_pairing_terms{};

    std::vector<GT> gt_bases;
    std::vector<Fr> gt_scalars;
    std::array<G1, 2> terminal_g1;
    std::array<G2, 2> terminal_g2;
    GT gt_result;
    GT terminal_pairing_result;
};
struct ValidationBreakdown {
    double crs_shape_ms{}, crs_digest_ms{}, crs_group_validation_ms{};
    double precomputation_shape_ms{}, precomputation_digest_binding_ms{}, precomputation_gt_validation_ms{};
    double statement_shape_ms{}, statement_digest_ms{}, statement_g1_validation_ms{}, statement_gt_validation_ms{};
    double proof_shape_ms{}, proof_g1_validation_ms{}, proof_g2_validation_ms{}, proof_gt_validation_ms{};
    double gt_subgroup_validation_ms{}, total_ms{};
    std::size_t g1_elements_checked{}, g2_elements_checked{}, gt_elements_checked{}, gt_subgroup_checks{};
    std::size_t crs_g1_checked{}, crs_g2_checked{}, statement_g1_checked{}, statement_gt_checked{};
    std::size_t precomputation_gt_checked{}, proof_g1_checked{}, proof_g2_checked{}, proof_gt_checked{};
};
struct VerifyCoreBreakdown {
    double transcript_replay_ms{}, transcript_serialization_ms{}, sha256_ms{}, challenge_to_field_ms{};
    double batch_inversion_ms{}, symbolic_initialization_ms{}, symbolic_dory_fold_ms{};
    double symbolic_rexp_fresh_ms{}, symbolic_batch_ms{}, normalization_ms{};
    double terminal_g1_g2_ms{}, gt_multiexp_ms{}, terminal_pairing_ms{}, identity_check_ms{};
    double final_dory_g1_g2_ms{}, final_dory_pairing_ms{}, final_rexp_preparation_ms{};
    double final_rexp_pairing_ms{}, other_profiling_overhead_ms{}, total_ms{};
    std::size_t transcript_entries{}, transcript_bytes_absorbed{}, sha256_calls{}, challenge_derivations{}, rejection_sampling_retries{};
    std::size_t normalization_calls{}, total_terms_before_normalization{}, total_terms_after_normalization{};
    std::size_t duplicate_terms_coalesced{}, zero_terms_removed{};
    std::size_t gt_multiexp_calls{}, terminal_pairing_terms{}, subgroup_validation_calls{};
    std::size_t fresh_rexp_instances{}, dory_fold_steps{}, batching_u_messages{}, gamma_challenges{};
    std::size_t rexp_challenges{}, dory_fold_challenges{}, eta_challenges{};
    std::size_t final_dory_checks{}, final_rexp_checks{}, combined_terminal_checks{};
};
class ValidatedVerificationInputs {
public:
    const CRS& crs() const { return crs_; }
    const Precomputation& precomputation() const { return precomputation_; }
    const Statement& statement() const { return statement_; }
    const Proof& proof() const { return proof_; }
private:
    ValidatedVerificationInputs(const CRS& crs, const Precomputation& precomputation,
                                const Statement& statement, const Proof& proof)
        : crs_(crs), precomputation_(precomputation),
          statement_(statement), proof_(proof) {}
    CRS crs_;
    Precomputation precomputation_;
    Statement statement_;
    Proof proof_;
    friend std::optional<ValidatedVerificationInputs> validate_verification_inputs(
      const CRS&,const Precomputation&,const Statement&,const Proof&,ValidationBreakdown*);
};
std::optional<ValidatedVerificationInputs> validate_verification_inputs(
    const CRS&,const Precomputation&,const Statement&,const Proof&,ValidationBreakdown* = nullptr);
bool verify_prevalidated(const ValidatedVerificationInputs&, VerifyCoreBreakdown* = nullptr);

bool verify_online(const CRS&, const Precomputation&, const Statement&, const Proof&);
bool verify_online_with_breakdown(const CRS&, const Precomputation&, const Statement&, const Proof&,
                                  VerifyCoreBreakdown&);
std::vector<Fr> batch_invert_nonzero(std::span<const Fr> values);
bool verify_reference(const CRS&, const Precomputation&, const Statement&, const Proof&);
bool verify(const CRS&, const Precomputation&, const Statement&, const Proof&);
bool verify_with_diagnostics(const CRS&, const Precomputation&, const Statement&, const Proof&,
                             VerifyDiagnostics&);
}

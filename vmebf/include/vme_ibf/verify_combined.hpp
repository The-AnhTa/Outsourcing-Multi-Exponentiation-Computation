#pragma once

#include "vme_ibf/verify_deferred.hpp"

namespace vme_ibf {

struct CombinedVerificationTrace {
    bool accepted{};
    std::vector<Fr> rho;
    std::vector<Fr> beta;
    std::vector<Fr> alpha;
    std::vector<Fr> gamma;
    Fr epsilon;
    Fr eta;
    Digest final_transcript_digest{};
    std::size_t gt_terms_before_normalize{};
    std::size_t gt_terms_after_normalize{};
    std::size_t pairing_terms_before_normalize{};
    std::size_t pairing_terms_after_normalize{};
    std::size_t gt_multiexp_calls{};
    std::size_t pairing_product_calls{};
    std::size_t intermediate_gt_exponentiations{};
    std::size_t underlying_pairings{};
    std::size_t coalesced_pairing_terms{};
    std::size_t miller_loop_batches{};
    std::size_t miller_loop_terms{};
    std::size_t final_exponentiations{};
    double transcript_ms{};
    double batch_inversion_ms{};
    double recurrence_ms{};
    double terminal_assembly_ms{};
    double combined_normalize_ms{};
    double gt_msm_pairing_ms{};
    double gt_msm_ms{};
    double multi_pairing_ms{};
    double total_ms{};
    GT evaluated_dory_residual;
    GT evaluated_rexp_residual;
    GT evaluated_combined_residual;
    GT evaluated_combined_reference;
};

bool verify_deferred_combined(const VmeIbfCRS&,
                              const VmeIbfPrecomputation&,
                              const VmeIbfStatement&,
                              const VmeIbfProof&);
CombinedVerificationTrace verify_deferred_combined_with_trace(
    const VmeIbfCRS&,
    const VmeIbfPrecomputation&,
    const VmeIbfStatement&,
    const VmeIbfProof&);

}

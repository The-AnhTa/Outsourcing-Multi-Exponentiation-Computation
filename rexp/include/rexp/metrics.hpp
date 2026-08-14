#pragma once

#include <cstddef>
#include <vector>

namespace rexp {

struct RexpProofValidationMetrics {
    std::size_t gt_elements_checked = 0;
    double canonical_decode_ms = 0;
    double gt_subgroup_validation_ms = 0;
};

struct RexpVerifyMetrics {
    std::size_t dory_verifications = 0, gt_multiexponentiations = 0;
    std::size_t dory_terminal_pairings = 0, final_pairings = 0;
    std::size_t actual_gt_bases = 0;
    double transcript_ms = 0, gt_multiexp_ms = 0;
    double dory_pairing_ms = 0, final_pairing_ms = 0, total_ms = 0;
#ifdef REXP_ENABLE_PROFILING
    struct PerDoryProfile {
        std::size_t round = 0, dimension = 0;
        std::size_t raw_gt_terms = 0, normalized_gt_terms = 0;
        std::size_t duplicate_bases_merged = 0;
        std::size_t zero_terms_removed = 0, identity_terms_removed = 0;
        std::size_t gt_multiexp_calls = 0, terminal_pairings = 0;
    };
    std::vector<PerDoryProfile> per_dory;
#endif
};

} 

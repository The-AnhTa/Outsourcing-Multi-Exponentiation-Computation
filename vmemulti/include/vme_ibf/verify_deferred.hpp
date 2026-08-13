#pragma once
#include "vme_ibf/proof.hpp"
#include "vme_ibf/symbolic_gt.hpp"

namespace vme_ibf {

struct DeferredVerificationTrace {
    bool accepted{};
    bool dory_accepted{};
    bool rexp_accepted{};
    std::vector<Fr> rho;
    std::vector<Fr> beta;
    std::vector<Fr> alpha;
    std::vector<Fr> gamma;
    Fr epsilon;
    std::size_t dory_gt_terms_before_normalize{};
    std::size_t dory_gt_terms_after_normalize{};
    std::size_t dory_pairing_terms{};
    std::size_t rexp_gt_terms_before_normalize{};
    std::size_t rexp_gt_terms_after_normalize{};
    std::size_t rexp_pairing_terms{};
    std::size_t gt_multiexp_calls{};
    std::size_t pairing_product_calls{};
    std::size_t ordinary_intermediate_gt_exponentiations{};
    GT dory_residual_reference;
    GT dory_residual_pippenger;
    GT rexp_residual_reference;
    GT rexp_residual_pippenger;
};

bool verify_deferred(const VmeIbfCRS&,
                     const VmeIbfPrecomputation&,
                     const VmeIbfStatement&,
                     const VmeIbfProof&);
DeferredVerificationTrace verify_deferred_with_trace(
    const VmeIbfCRS&,
    const VmeIbfPrecomputation&,
    const VmeIbfStatement&,
    const VmeIbfProof&);

}

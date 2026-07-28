#pragma once
#include "vme_ibf/proof.hpp"
#include "vme_ibf/symbolic_gt.hpp"
namespace vme_ibf {
struct DeferredVerificationTrace {bool accepted{},dory_accepted{},rexp_accepted{};std::vector<Fr>rho,beta,alpha,gamma;Fr epsilon;size_t dory_gt_terms_before_normalize{},dory_gt_terms_after_normalize{},dory_pairing_terms{},rexp_gt_terms_before_normalize{},rexp_gt_terms_after_normalize{},rexp_pairing_terms{},gt_multiexp_calls{},pairing_product_calls{},ordinary_intermediate_gt_exponentiations{};GT dory_residual_reference,dory_residual_pippenger,rexp_residual_reference,rexp_residual_pippenger;};
bool verify_deferred(const VmeIbfCRS&,const VmeIbfPrecomputation&,const VmeIbfStatement&,const VmeIbfProof&);
DeferredVerificationTrace verify_deferred_with_trace(const VmeIbfCRS&,const VmeIbfPrecomputation&,const VmeIbfStatement&,const VmeIbfProof&);
}

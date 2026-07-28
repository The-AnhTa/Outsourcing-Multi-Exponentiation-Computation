#pragma once
#include "vpip_bf/verify_deferred.hpp"
namespace vpip_bf {
struct CombinedVerificationTrace {bool accepted{};std::vector<Fr>rho,beta,alpha,gamma;Fr epsilon,eta;Digest final_transcript_digest{};size_t gt_terms_before_normalize{},gt_terms_after_normalize{},pairing_terms_before_normalize{},pairing_terms_after_normalize{},gt_multiexp_calls{},pairing_product_calls{},intermediate_gt_exponentiations{},underlying_pairings{},coalesced_pairing_terms{},miller_loop_batches{},miller_loop_terms{},final_exponentiations{};double transcript_ms{},batch_inversion_ms{},recurrence_ms{},terminal_assembly_ms{},combined_normalize_ms{},gt_msm_pairing_ms{},gt_msm_ms{},multi_pairing_ms{},total_ms{};GT evaluated_dory_residual,evaluated_rexp_residual,evaluated_combined_residual,evaluated_combined_reference;};
bool verify_deferred_combined(const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&);
CombinedVerificationTrace verify_deferred_combined_with_trace(const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&);
}

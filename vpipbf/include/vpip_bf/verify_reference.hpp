#pragma once
#include "vpip_bf/proof.hpp"

namespace vpip_bf {
struct ReferenceVerificationTrace {
  bool accepted{}, dory_accepted{}, rexp_accepted{};
  std::vector<Fr> rho, beta, alpha, gamma;
  Fr epsilon;
  DoryTargetState final_aggregate;
  GT final_rexp_d1;
  GT dory_residual;
  GT rexp_residual;
};
struct OnlineTimingBreakdown {
  double transcript_init_ms{}, transcript_replay_ms{}, challenge_derivation_ms{};
  double tensor_reconstruction_ms{}, g1_msm_y_ms{}, dory_recurrence_ms{};
  double gt_multiexp_ms{}, terminal_dory_ms{}, terminal_rexp_ms{};
  double other_online_ms{}, online_verify_ms{};
};

bool verify_reference(const VpipBfCRS&, const VpipBfPrecomputation&, const VpipBfStatement&, const VpipBfProof&);
ReferenceVerificationTrace verify_reference_diagnostic(const VpipBfCRS&, const VpipBfPrecomputation&, const VpipBfStatement&, const VpipBfProof&);
bool validate_verification_inputs(const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&);
ReferenceVerificationTrace verify_core_unchecked(const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&,OnlineTimingBreakdown* = nullptr);
void reset_verification_core_call_count_for_testing();
std::size_t verification_core_call_count_for_testing();
}

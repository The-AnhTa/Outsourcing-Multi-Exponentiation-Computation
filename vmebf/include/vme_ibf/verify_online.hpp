#pragma once
#include "vme_ibf/verify_combined.hpp"
namespace vme_ibf {
struct ValidatedVerificationInputs {
  const VmeIbfCRS* crs{};
  const VmeIbfPrecomputation* precomp{};
  const VmeIbfStatement* statement{};
  const VmeIbfProof* proof{};
};
bool prepare_validated_verification_inputs(const VmeIbfCRS&,const VmeIbfPrecomputation&,const VmeIbfStatement&,const VmeIbfProof&,ValidatedVerificationInputs&);
bool verify_online(const ValidatedVerificationInputs&);
CombinedVerificationTrace verify_online_with_trace(const ValidatedVerificationInputs&);
}

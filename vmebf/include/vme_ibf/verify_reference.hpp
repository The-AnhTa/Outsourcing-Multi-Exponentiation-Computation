#pragma once
#include "vme_ibf/proof.hpp"

namespace vme_ibf {
struct ReferenceVerificationTrace {
  bool accepted{}, dory_accepted{}, rexp_accepted{};
  std::vector<Fr> rho, beta, alpha, gamma;
  Fr epsilon;
  DoryTargetState final_aggregate;
  GT final_rexp_d1;
  GT dory_residual;
  GT rexp_residual;
};

bool verify_reference(const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeIbfStatement&, const VmeIbfProof&);
ReferenceVerificationTrace verify_reference_diagnostic(const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeIbfStatement&, const VmeIbfProof&);
}

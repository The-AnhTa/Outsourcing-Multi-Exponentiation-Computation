#pragma once
#include "vme_ibf/verify_combined.hpp"
#include <optional>
namespace vme_ibf {
class ValidatedVerificationInputs {
 public:
  const VmeIbfCRS& crs() const { return crs_; }
  const VmeIbfPrecomputation& precomputation() const { return precomputation_; }
  const VmeIbfStatement& statement() const { return statement_; }
  const VmeIbfProof& proof() const { return proof_; }
 private:
  ValidatedVerificationInputs(const VmeIbfCRS& crs,
                              const VmeIbfPrecomputation& precomputation,
                              const VmeIbfStatement& statement,
                              const VmeIbfProof& proof)
      : crs_(crs), precomputation_(precomputation),
        statement_(statement), proof_(proof) {}
  VmeIbfCRS crs_;
  VmeIbfPrecomputation precomputation_;
  VmeIbfStatement statement_;
  VmeIbfProof proof_;
  friend std::optional<ValidatedVerificationInputs>
  validate_verification_inputs(const VmeIbfCRS&,
                               const VmeIbfPrecomputation&,
                               const VmeIbfStatement&,
                               const VmeIbfProof&);
};
std::optional<ValidatedVerificationInputs> validate_verification_inputs(
    const VmeIbfCRS&, const VmeIbfPrecomputation&,
    const VmeIbfStatement&, const VmeIbfProof&);
bool verify_online(const ValidatedVerificationInputs&);
CombinedVerificationTrace verify_online_with_trace(const ValidatedVerificationInputs&);
}

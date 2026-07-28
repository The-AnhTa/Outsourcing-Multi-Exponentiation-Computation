#pragma once
#include "vme_ibf/phase2.hpp"
#include "vme_ibf/proof.hpp"
#include "vme_ibf/verify_reference.hpp"
#include "vme_ibf/verify_combined.hpp"

namespace vme_ibf {



struct VmeMultiStatement {
  std::vector<std::vector<Fr>> x_instances;
  std::vector<G2> X_instances;
};

struct AggregatedInstance {
  std::vector<Fr> z;
  G2 Y;
};

struct VmeMultiProverTrace {
  Fr aggregate_gamma;
  AggregatedInstance aggregate;
  Phase1Result phase1;
  Phase2Result phase2;
};

struct ValidatedVmeMultiInputs {
  const VmeIbfCRS* crs{};
  const VmeIbfPrecomputation* precomputation{};
  const VmeMultiStatement* statement{};
  const VmeIbfProof* proof{};
  Digest precomputation_binding_digest{};
};

struct VmeMultiOnlineVerificationTrace {
  bool accepted{};
  double transcript_prefix_ms{}, aggregate_gamma_ms{}, aggregate_z_ms{},
      aggregate_Y_ms{}, core_input_construction_ms{}, total_ms{};
  CombinedVerificationTrace core;
};

bool validate_vmemulti_statement(const VmeIbfCRS&, const VmeMultiStatement&);
Digest compute_vmemulti_binding_digest(const VmeIbfCRS&,
                                       const VmeIbfPrecomputation&,
                                       const VmeMultiStatement&);
Fr derive_vmemulti_aggregate_gamma(const VmeIbfCRS&,
                                   const VmeIbfPrecomputation&,
                                   const VmeMultiStatement&,
                                   Digest* transcript_after_gamma = nullptr);
AggregatedInstance aggregate_instances(
    const std::vector<std::vector<Fr>>& x_instances,
    const std::vector<G2>& X_instances, const Fr& gamma);
VmeMultiProverTrace prove_vmemulti_with_trace(const VmeIbfCRS&,
                                              const VmeIbfPrecomputation&,
                                              const VmeMultiStatement&);
VmeIbfProof prove_vmemulti(const VmeIbfCRS&, const VmeIbfPrecomputation&,
                           const VmeMultiStatement&);
ReferenceVerificationTrace verify_vmemulti_diagnostic(
    const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeMultiStatement&,
    const VmeIbfProof&);
CombinedVerificationTrace verify_vmemulti_combined_diagnostic(
    const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeMultiStatement&,
    const VmeIbfProof&);
bool verify_vmemulti(const VmeIbfCRS&, const VmeIbfPrecomputation&,
                     const VmeMultiStatement&, const VmeIbfProof&);
bool prepare_validated_vmemulti_inputs(
    const VmeIbfCRS&, const VmeIbfPrecomputation&, const VmeMultiStatement&,
    const VmeIbfProof&, ValidatedVmeMultiInputs&);
bool verify_vmemulti_online(const ValidatedVmeMultiInputs&);
VmeMultiOnlineVerificationTrace verify_vmemulti_online_with_trace(
    const ValidatedVmeMultiInputs&);

}

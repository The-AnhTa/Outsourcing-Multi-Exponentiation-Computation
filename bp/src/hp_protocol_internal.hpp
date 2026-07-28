#pragma once

#include "bp/helper_prover.hpp"

namespace bp::hp_internal {

struct HpAggregate {
  Group X_gamma;
  std::vector<Scalar> z_gamma;
  Scalar c_gamma;
  Scalar x0;
  Scalar y0;
};

struct HpRoundTrace {
  std::vector<Scalar> alphas;
  Scalar x0;
  Scalar y0;
  Scalar gamma;
  HpAggregate aggregate;
  Digest application_context{};
};

HpRoundTrace reconstruct_hp_trace(const HpPublicParams&, const Group& Z,
                                  std::span<const Scalar> x,
                                  std::span<const Scalar> y,
                                  std::span<const RoundMessage> rounds);
Group direct_hp_aggregate_msm_for_test(const HpPublicParams&,
                                       const HpAggregate&);
bool verify_hp_reference_for_test(const HpPublicParams&, const HpClientPrecomp&,
                                  const Group& Z, std::span<const Scalar> x,
                                  std::span<const Scalar> y,
                                  const HpProof&);
bool verify_hp_wrong_round_encoding_for_test(
    const HpPublicParams&, const HpClientPrecomp&, const Group& Z,
    std::span<const Scalar> x, std::span<const Scalar> y, const HpProof&);

}

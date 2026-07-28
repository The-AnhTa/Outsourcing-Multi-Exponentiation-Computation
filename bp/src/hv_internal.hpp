#pragma once

#include "bp/helper_verifier.hpp"

namespace bp::hv_internal {

struct PreparedStatement {
  std::vector<Scalar> alphas;
  std::vector<Scalar> alpha_inverses;
  std::vector<Scalar> g;
  std::vector<Scalar> h;
  std::vector<Scalar> z_v;
  Group X0;
  Digest context{};
};

bool prepare_statement_for_test(const HvPublicParams&, const HvStatement&,
                                PreparedStatement&);
bool verify_helper_verifier_reference_for_test(
    const HvPublicParams&, const HvVerifierPrecomputation&,
    const HvStatement&, const HvProof&);
Group direct_vme_relation_for_test(const HvPublicParams&,
                                   const PreparedStatement&);

}

#pragma once

#include "bp/bp.hpp"

namespace bp::internal {

bool replay_bp_challenges(const PublicParams&, const Group&, const Proof&,
                          std::vector<Scalar>&) noexcept;
bool replay_bp_challenges_prevalidated(
    const PublicParams&, const Group&, const Proof&,
    std::vector<Scalar>&) noexcept;
bool replay_bp_challenges_prevalidated(
    const PublicParams&, const Digest& public_parameter_digest,
    const Group&, const Proof&, std::vector<Scalar>&) noexcept;
Digest bp_public_parameter_digest(const PublicParams&);

}

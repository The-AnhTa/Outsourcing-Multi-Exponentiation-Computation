#pragma once

#include "rexpbf/prove.hpp"

namespace rexpbf::internal {

ChallengeTrace replay_challenges_prevalidated(
    const CRS& crs,
    const Statement& statement,
    const Proof& proof,
    TranscriptMetrics* metrics = nullptr);

} // namespace rexpbf::internal

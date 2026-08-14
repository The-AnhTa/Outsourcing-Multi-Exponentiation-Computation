#pragma once

#include "rexp/dory.hpp"

namespace rexp::internal {

Digest dory_initial(const DoryCRS& crs, const DoryStatement& statement);
Digest dory_batch_initial(
    const DoryCRS& crs, const std::vector<DoryStatement>& statements);
Digest dory_absorb_beta(
    const Digest& current, std::size_t round, std::size_t dimension,
    const DoryRound& message);
Digest dory_absorb_alpha(
    const Digest& beta_digest, std::size_t round, std::size_t dimension,
    const DoryRound& message);
Digest dory_absorb_final(
    const Digest& current, const G1& phi, const G2& theta);
Digest dory_absorb_merge(
    const Digest& current, std::size_t instance, const GT& cross_term);
Digest dory_enter_batch(const Digest& current);

struct DoryTranscriptReplay {
    DoryChallenges challenges;
    Digest end{};
};

DoryTranscriptReplay replay_dory_transcript(
    const DoryCRS& crs, const DoryProof& proof, Digest transcript);

} 

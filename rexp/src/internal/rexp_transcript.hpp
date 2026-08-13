#pragma once

#include "rexp/proof.hpp"

namespace rexp::internal {

Digest rexp_crs_digest(const RawRexpCRS& crs);
Digest rexp_statement_digest(
    const PreparedPublicParameters& params,
    const std::vector<G1>& points,
    const GT& d1,
    const GT& e,
    const GT& f,
    const GT& tl,
    const GT& tr);
Digest rexp_initial_transcript(
    const PreparedPublicParameters& params,
    const PreparedStatement& statement);
Digest rexp_absorb_round(
    const Digest& transcript,
    std::size_t round,
    std::size_t dimension,
    const RexpRoundMessage& message);
Digest rexp_enter_dory(
    const Digest& transcript,
    std::size_t round,
    std::size_t half_dimension,
    const DoryStatement& statement);
Digest rexp_leave_dory(const Digest& transcript, std::size_t round);
Digest rexp_absorb_final(const Digest& transcript, const G1& result);
RexpRoundMessage rexp_initial_round_message(const PreparedStatement& statement);

} // namespace rexp::internal

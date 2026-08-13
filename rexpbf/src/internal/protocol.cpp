#include "protocol.hpp"

#include <cstdint>

namespace rexpbf::internal {

void initialize_protocol_transcript(
    Transcript& transcript, const CRS& crs, const Statement& statement) {
    transcript.append_bytes("crs-digest", crs.digest);
    transcript.append_bytes("statement-digest", statement.digest);
    transcript.append_u64("d", crs.d);
    transcript.append_u64("n", crs.n);
}

void append_round_metadata(
    Transcript& transcript,
    std::string_view phase,
    std::size_t level,
    std::size_t outer_round,
    std::size_t current_dimension,
    std::size_t next_dimension) {
    transcript.append_bytes(
        "phase",
        {reinterpret_cast<const std::uint8_t*>(phase.data()), phase.size()});
    transcript.append_u64("level", level);
    transcript.append_u64("outer-round", outer_round);
    transcript.append_u64("current-dimension", current_dimension);
    transcript.append_u64("next-dimension", next_dimension);
}

Fr derive_initial_rho(Transcript& transcript, const CRS& crs) {
    append_round_metadata(
        transcript, "REXP-BF-G1-INITIAL-REXP-V1", 1, 0, crs.n, crs.n / 2);
    return transcript.challenge_nonzero_fr("REXP-BF-G1-RHO-V1", 0);
}

} // namespace rexpbf::internal

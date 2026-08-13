#pragma once

#include "rexpbf/transcript.hpp"
#include "rexpbf/types.hpp"

#include <cstddef>
#include <string_view>

namespace rexpbf::internal {

inline constexpr std::string_view transcript_domain = "REXP-BF-G1-FS-v1";

void initialize_protocol_transcript(
    Transcript& transcript, const CRS& crs, const Statement& statement);

void append_round_metadata(
    Transcript& transcript,
    std::string_view phase,
    std::size_t level,
    std::size_t outer_round,
    std::size_t current_dimension,
    std::size_t next_dimension);

Fr derive_initial_rho(Transcript& transcript, const CRS& crs);

} // namespace rexpbf::internal

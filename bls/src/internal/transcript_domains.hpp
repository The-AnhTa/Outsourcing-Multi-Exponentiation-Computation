#pragma once

#include "blsagg/protocol.hpp"

#include <string_view>

namespace blsagg::internal::transcript_domain {

inline constexpr std::string_view protocol = "bls-agg-bf-nonzk-v1";
inline constexpr std::string_view curve = "BN254/mcl-v3.00/e:G1xG2->GT";
inline constexpr std::string_view absorb_operation = "absorb";
inline constexpr std::string_view challenge_operation = "challenge-nonzero";

// Largest acceptable SHA-256 digest for unbiased reduction into the BN254
// scalar field. Values below this limit are rejection-sampled.
inline constexpr Digest scalar_rejection_limit = {
    0xde, 0xd4, 0x5b, 0x0d, 0x80, 0x00, 0x00, 0x0a,
    0x5d, 0x39, 0xd1, 0x00, 0x00, 0x00, 0x00, 0x2f,
    0xfd, 0xbd, 0x00, 0x00, 0x00, 0x00, 0x00, 0x63,
    0xc6, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x4e};

}  // namespace blsagg::internal::transcript_domain

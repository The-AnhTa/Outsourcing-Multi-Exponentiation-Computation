#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace pippenger::internal {

inline constexpr std::size_t kCanonicalScalarBytes = 32;
inline constexpr std::size_t kMaximumWindowWidth = 16;
using CanonicalScalar =
    std::array<std::uint8_t, kCanonicalScalarBytes>;

}  // namespace pippenger::internal

namespace pinkas::internal {

inline constexpr std::string_view kProtocolVersion = "PinkasFS128BN-v2";
inline constexpr std::string_view kGroupIdentifier = "BN254/G2/mcl-v3.00";
inline constexpr std::size_t kChallengeBitLength = 128;
inline constexpr std::size_t kChallengeBytes = kChallengeBitLength / 8;
inline constexpr std::size_t kDigestBytes = 32;
using Digest = std::array<std::uint8_t, kDigestBytes>;
using CanonicalScalar = pippenger::internal::CanonicalScalar;

}  // namespace pinkas::internal

#pragma once
#include "vpip_bf/setup.hpp"
#include "vpip_bf/transcript.hpp"
namespace vpip_bf::internal {
RexpClaims rexp_claim(std::size_t, const VpipBfPrecomputation&,
                      std::span<const RexpClaims>);
void absorb_rexp_claim(Transcript&, std::size_t, std::size_t,
                       const RexpClaims&);
}

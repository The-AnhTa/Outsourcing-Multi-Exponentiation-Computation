#include "protocol.hpp"
#include <array>
#include <stdexcept>
namespace vpip_bf::internal {
RexpClaims rexp_claim(std::size_t round, const VpipBfPrecomputation& p,
                      std::span<const RexpClaims> dynamic) {
    if (round == 0) return {p.delta1R.at(0), p.delta2R.at(0),
                            p.pairing_x.at(1), p.delta1R.at(0)};
    if (round - 1 >= dynamic.size())
        throw std::invalid_argument("missing dynamic REXP claim");
    return dynamic[round - 1];
}
void absorb_rexp_claim(Transcript& transcript, std::size_t round,
                       std::size_t dimension, const RexpClaims& claim) {
    std::array<Bytes, 6> fields{encode_u64_be(round),
        encode_u64_be(dimension), serialize(claim.E), serialize(claim.F),
        serialize(claim.TL), serialize(claim.TR)};
    transcript.absorb("vpipbf/rexp-message/v1", fields);
}
} // namespace vpip_bf::internal

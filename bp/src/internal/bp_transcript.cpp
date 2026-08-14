#include "internal/bp_transcript.hpp"

#include "internal/protocol_utils.hpp"

#include <algorithm>
#include <array>

namespace bp::internal {
namespace {

constexpr std::string_view kProtocol = "BP-IPA/BN254-G2/v1";
constexpr std::string_view kAlphaDomain = "BP-IPA-ALPHA-v1";
constexpr std::string_view kRoundDomain = "BP-IPA-ROUND-v1";
constexpr Digest kRejectionLimit = {
    0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,
    0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,
    0xfd,0xbd,0x00,0x00,0x00,0x00,0x00,0x63,
    0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};

}  

Digest bp_transcript_parameter_digest(const PublicParams& pp) {
  Bytes input;
  frame(input, "BP-IPA-PUBLIC-PARAMS-v1");
  frame(input, pp.group_identifier);
  frame(input, pp.scalar_modulus);
  frame(input, pp.hash_suite_identifier);
  frame(input, pp.transcript_domain);
  frame(input, u64be(pp.n));
  frame(input, u64be(pp.d));
  for (const auto& point : pp.G) frame(input, serialize(point));
  for (const auto& point : pp.H) frame(input, serialize(point));
  frame(input, serialize(pp.K));
  return sha256(input);
}

BpTranscript::BpTranscript(const PublicParams& pp, const Group& statement,
                           const Digest* cached_parameter_digest) {
  if (pp.transcript_domain == kHpTranscriptDomain) {
    hp_.emplace(pp, statement);
    return;
  }
  Bytes input;
  frame(input, kProtocol);
  frame(input, pp.transcript_domain);
  frame(input, pp.group_identifier);
  frame(input, pp.scalar_modulus);
  frame(input, pp.hash_suite_identifier);
  frame(input, cached_parameter_digest ? *cached_parameter_digest
                                       : bp_transcript_parameter_digest(pp));
  frame(input, serialize(statement));
  state_ = sha256(input);
}

Scalar BpTranscript::round(std::size_t round_index, const Group& left,
                           const Group& right) {
  if (hp_) return hp_->round(round_index, left, right);
  Bytes message;
  frame(message, kRoundDomain);
  frame(message, state_);
  frame(message, u64be(round_index));
  frame(message, serialize(left));
  frame(message, serialize(right));
  const Digest round_state = sha256(message);
  const Scalar challenge = hash_to_nonzero(round_state);
  Bytes next;
  frame(next, kProtocol);
  frame(next, round_state);
  frame(next, serialize(challenge));
  state_ = sha256(next);
  return challenge;
}

Scalar BpTranscript::hash_to_nonzero(const Digest& input) {
  for (std::uint64_t counter = 0;; ++counter) {
    Bytes candidate;
    frame(candidate, kAlphaDomain);
    frame(candidate, input);
    frame(candidate, u64be(counter));
    const Digest hash = sha256(candidate);
    if (!std::lexicographical_compare(hash.begin(), hash.end(),
                                      kRejectionLimit.begin(),
                                      kRejectionLimit.end()))
      continue;
    Scalar out;
    out.setBigEndianMod(hash.data(), hash.size());
    if (!out.isZero()) return out;
  }
}

Group derive_bp_generator(std::string_view domain,
                          std::span<const std::uint8_t> seed,
                          std::uint64_t index, bool indexed) {
  for (std::uint64_t counter = 0;; ++counter) {
    Bytes input;
    frame(input, kProtocol);
    frame(input, domain);
    frame(input, seed);
    if (indexed) frame(input, u64be(index));
    frame(input, u64be(counter));
    Group out;
    mcl::bn::hashAndMapToG2(out, input.data(), input.size());
    if (valid_group(out, true)) return out;
  }
}

}  

#include "hp_transcript_internal.hpp"

#include <algorithm>
#include <stdexcept>

namespace bp::hp_internal {
namespace {

constexpr Digest kRejectionLimit = {
    0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,
    0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,
    0xfd,0xbd,0x00,0x00,0x00,0x00,0x00,0x63,
    0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};

Bytes u64be(std::uint64_t value) {
  Bytes out;
  for (int shift = 56; shift >= 0; shift -= 8)
    out.push_back(static_cast<std::uint8_t>(value >> shift));
  return out;
}

void raw(Bytes& out, std::span<const std::uint8_t> field) {
  out.insert(out.end(), field.begin(), field.end());
}

void frame(Bytes& out, std::span<const std::uint8_t> field) {
  raw(out, u64be(field.size()));
  raw(out, field);
}

void frame(Bytes& out, std::string_view field) {
  frame(out, {reinterpret_cast<const std::uint8_t*>(field.data()), field.size()});
}

}

HpBpTranscript::HpBpTranscript(const PublicParams& pp, const Group& Z) {
  if (pp.transcript_domain != kHpTranscriptDomain || !pp.transcript_crs_digest)
    throw std::invalid_argument("not HP Bulletproof parameters");
  Bytes initial;
  frame(initial, kHpTranscriptDomain);
  frame(initial, *pp.transcript_crs_digest);
  frame(initial, u64be(pp.n));
  frame(initial, serialize(Z));
  state_ = sha256(initial);
}

void HpBpTranscript::absorb(std::string_view label,
                            std::span<const Bytes> fields) {
  Bytes input;
  frame(input, "BPVME/HP/TRANSCRIPT/ABSORB/v1");
  frame(input, state_);
  frame(input, label);
  for (const auto& field : fields) frame(input, field);
  state_ = sha256(input);
}

void HpBpTranscript::absorb(std::string_view label, const Bytes& field) {
  absorb(label, std::span<const Bytes>(&field, 1));
}

Scalar HpBpTranscript::challenge_nonzero(std::string_view label) {
  for (std::uint64_t counter = 0;; ++counter) {
    Bytes input;
    frame(input, "BPVME/HP/TRANSCRIPT/CHALLENGE/v1");
    frame(input, state_);
    frame(input, label);
    frame(input, u64be(counter));
    const Digest hash = sha256(input);
    if (!std::lexicographical_compare(
            hash.begin(), hash.end(), kRejectionLimit.begin(), kRejectionLimit.end()))
      continue;
    Scalar result;
    result.setBigEndianMod(hash.data(), hash.size());
    if (!result.isZero()) return result;
  }
}

Scalar HpBpTranscript::round(std::size_t k, const Group& A, const Group& B) {
  return round_encoded(u64be(k), A, B);
}

Scalar HpBpTranscript::round_encoded(const Bytes& index, const Group& A,
                                     const Group& B) {
  std::array<Bytes, 3> fields{index, serialize(A), serialize(B)};
  absorb("BPVME/HP/BP-ROUND/v1", fields);
  Scalar alpha = challenge_nonzero("bp-alpha");
  absorb("bp-alpha-value", serialize(alpha));
  return alpha;
}

Scalar HpBpTranscript::round_with_index_bytes_for_test(
    std::span<const std::uint8_t> index, const Group& A, const Group& B) {
  return round_encoded(Bytes(index.begin(), index.end()), A, B);
}

Scalar HpBpTranscript::outer_gamma(const Scalar& x0, const Scalar& y0) {
  std::array<Bytes, 2> fields{serialize(x0), serialize(y0)};
  absorb("BPVME/HP/OUTER-BATCH/v1", fields);
  return challenge_nonzero("hp-gamma");
}

}

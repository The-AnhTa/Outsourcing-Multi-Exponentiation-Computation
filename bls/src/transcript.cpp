#include "blsagg/transcript.hpp"
#include "internal/crypto.hpp"
#include "internal/transcript_domains.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <chrono>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace blsagg {

void initialize() {
  static std::once_flag once;
  std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}
Bytes encode_u64(std::uint64_t x) {
  Bytes out;
  for (int s = 56; s >= 0; s -= 8) out.push_back(static_cast<std::uint8_t>(x >> s));
  return out;
}
void append_frame(Bytes& out, std::span<const std::uint8_t> x) {
  internal::append_raw(out, encode_u64(x.size()));
  internal::append_raw(out, x);
}
void append_frame(Bytes& out, std::string_view x) {
  append_frame(out, {reinterpret_cast<const std::uint8_t*>(x.data()), x.size()});
}
Digest sha256(std::span<const std::uint8_t> in) {
  if (in.size() > std::numeric_limits<std::uint32_t>::max()) throw std::length_error("hash input");
  Digest out{};
  if (mcl::fp::sha256(out.data(), static_cast<std::uint32_t>(out.size()), in.data(),
                      static_cast<std::uint32_t>(in.size())) != out.size())
    throw std::runtime_error("sha256 failed");
  return out;
}
Bytes serialize(const Fr& x) { return internal::encode(x); }
Bytes serialize(const G1& x) { return internal::encode(x); }
Bytes serialize(const G2& x) { return internal::encode(x); }
Bytes serialize(const GT& x) { return internal::encode(x); }

Transcript::Transcript(const PublicParameters& pp, const Statement& s,
                       std::span<const G1> points) {
  const auto start=std::chrono::steady_clock::now();
  Bytes in;
  append_frame(in, internal::transcript_domain::protocol);
  append_frame(in, internal::transcript_domain::curve);
  append_frame(in, pp.mode == AggregationMode::BasicDistinct ? "basic-distinct" : "augmented");
  append_frame(in, encode_u64(pp.k));
  append_frame(in, encode_u64(pp.d));
  append_frame(in, pp.digest);
  append_frame(in, serialize(pp.H));
  append_frame(in, serialize(s.sigma_agg));
  append_frame(in, encode_u64(s.messages.size()));
  for (std::size_t i = 0; i < s.messages.size(); ++i) {
    append_frame(in, encode_u64(i));
    append_frame(in, s.messages[i]);
  }
  append_frame(in, encode_u64(s.public_keys.size()));
  for (std::size_t i = 0; i < s.public_keys.size(); ++i) {
    append_frame(in, encode_u64(i));
    append_frame(in, serialize(s.public_keys[i]));
  }
  append_frame(in, encode_u64(points.size()));
  for (std::size_t i = 0; i < points.size(); ++i) {
    append_frame(in, encode_u64(i));
    append_frame(in, serialize(points[i]));
  }
  state_ = sha256(in);
  timing_ms_+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
}
void Transcript::absorb(std::string_view label, std::uint64_t index,
                        std::span<const Bytes> fields) {
  const auto start=std::chrono::steady_clock::now();
  Bytes in;
  append_frame(in, internal::transcript_domain::protocol);
  append_frame(in, internal::transcript_domain::absorb_operation);
  append_frame(in, label);
  append_frame(in, encode_u64(index));
  append_frame(in, state_);
  append_frame(in, encode_u64(fields.size()));
  for (const auto& f : fields) append_frame(in, f);
  state_ = sha256(in);
  timing_ms_+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
}
void Transcript::absorb(std::string_view label, std::uint64_t index, const Bytes& field) {
  absorb(label, index, std::span<const Bytes>(&field, 1));
}
Fr Transcript::challenge_nonzero(std::string_view label, std::uint64_t index) {
  for (std::uint64_t counter = 0;; ++counter) {
    const auto start=std::chrono::steady_clock::now();
    Bytes in;
    append_frame(in, internal::transcript_domain::protocol);
    append_frame(in, internal::transcript_domain::challenge_operation);
    append_frame(in, label);
    append_frame(in, encode_u64(index));
    append_frame(in, state_);
    append_frame(in, encode_u64(counter));
    const auto h = sha256(in);
    if (!std::lexicographical_compare(
            h.begin(), h.end(),
            internal::transcript_domain::scalar_rejection_limit.begin(),
            internal::transcript_domain::scalar_rejection_limit.end())) continue;
    Fr out;
    out.setBigEndianMod(h.data(), h.size());
    if (out.isZero()) continue;
    timing_ms_+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    absorb(std::string(label) + "/value", index, serialize(out));
    return out;
  }
}
}

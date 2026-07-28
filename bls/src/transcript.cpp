#include "blsagg/transcript.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <chrono>
#include <limits>
#include <mutex>
#include <stdexcept>

namespace blsagg {
namespace {
constexpr std::string_view kProtocol = "bls-agg-bf-nonzk-v1";
constexpr std::string_view kCurve = "BN254/mcl-v3.00/e:G1xG2->GT";

constexpr Digest kLimit = {0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,
  0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,0xfd,0xbd,0x00,0x00,
  0x00,0x00,0x00,0x63,0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};

template<class T> Bytes encode(const T& x) {
  Bytes out(2048);
  const auto n = x.serialize(out.data(), out.size());
  if (!n) throw std::runtime_error("mcl serialization failed");
  out.resize(n);
  return out;
}
void raw(Bytes& out, std::span<const std::uint8_t> x) {
  out.insert(out.end(), x.begin(), x.end());
}
}

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
  raw(out, encode_u64(x.size()));
  raw(out, x);
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
Bytes serialize(const Fr& x) { return encode(x); }
Bytes serialize(const G1& x) { return encode(x); }
Bytes serialize(const G2& x) { return encode(x); }
Bytes serialize(const GT& x) { return encode(x); }

Transcript::Transcript(const PublicParameters& pp, const Statement& s,
                       std::span<const G1> points) {
  const auto start=std::chrono::steady_clock::now();
  Bytes in;
  append_frame(in, kProtocol);
  append_frame(in, kCurve);
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
  append_frame(in, kProtocol);
  append_frame(in, "absorb");
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
    append_frame(in, kProtocol);
    append_frame(in, "challenge-nonzero");
    append_frame(in, label);
    append_frame(in, encode_u64(index));
    append_frame(in, state_);
    append_frame(in, encode_u64(counter));
    const auto h = sha256(in);
    if (!std::lexicographical_compare(h.begin(), h.end(), kLimit.begin(), kLimit.end())) continue;
    Fr out;
    out.setBigEndianMod(h.data(), h.size());
    if (out.isZero()) continue;
    timing_ms_+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
    absorb(std::string(label) + "/value", index, serialize(out));
    return out;
  }
}
}

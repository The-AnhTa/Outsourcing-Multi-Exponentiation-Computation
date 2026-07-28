#pragma once
#include "blsagg/protocol.hpp"

namespace blsagg {

class Transcript {
 public:
  Transcript(const PublicParameters&, const Statement&, std::span<const G1> message_points);
  void absorb(std::string_view label, std::uint64_t index, std::span<const Bytes> fields);
  void absorb(std::string_view label, std::uint64_t index, const Bytes& field);
  Fr challenge_nonzero(std::string_view label, std::uint64_t index);
  Digest digest() const { return state_; }
  double timing_ms() const { return timing_ms_; }
 private:
  Digest state_{};
  double timing_ms_{};
};

Bytes encode_u64(std::uint64_t);
Bytes serialize(const Fr&);
Bytes serialize(const G1&);
Bytes serialize(const G2&);
Bytes serialize(const GT&);
Digest sha256(std::span<const std::uint8_t>);
void append_frame(Bytes&, std::span<const std::uint8_t>);
void append_frame(Bytes&, std::string_view);

}

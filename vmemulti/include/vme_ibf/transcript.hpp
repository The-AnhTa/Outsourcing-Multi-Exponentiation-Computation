#pragma once
#include "vme_ibf/group_utils.hpp"
#include <string_view>

namespace vme_ibf {
class Transcript {
 public:
  explicit Transcript(std::span<const std::uint8_t> statement_digest);
  static Transcript resume(const Digest& current_state);
  void absorb(std::string_view label, std::span<const Bytes> fields);
  void absorb(std::string_view label, const Bytes& field);
  Fr challenge_nonzero(std::string_view label, std::uint64_t logical_index);
  Digest digest() const { return state_; }
 private:
  struct ResumeTag {};
  Transcript(const Digest& state, ResumeTag) : state_(state) {}
  Digest state_{};
};
}

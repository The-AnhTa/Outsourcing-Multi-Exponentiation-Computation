#pragma once

#include "bp/bp.hpp"
#include "hp_transcript_internal.hpp"

#include <optional>
#include <string_view>

namespace bp::internal {

Digest bp_transcript_parameter_digest(const PublicParams& pp);

class BpTranscript {
 public:
  BpTranscript(const PublicParams& pp, const Group& statement,
               const Digest* cached_parameter_digest = nullptr);
  Scalar round(std::size_t round_index, const Group& left,
               const Group& right);

 private:
  static Scalar hash_to_nonzero(const Digest& input);
  Digest state_{};
  std::optional<hp_internal::HpBpTranscript> hp_;
};

Group derive_bp_generator(std::string_view domain,
                          std::span<const std::uint8_t> seed,
                          std::uint64_t index, bool indexed);

}  

#pragma once

#include "bp/bp.hpp"

namespace bp::hp_internal {

class HpBpTranscript {
 public:
  HpBpTranscript(const PublicParams& pp, const Group& Z);
  Scalar round(std::size_t k, const Group& A, const Group& B);
  Scalar round_with_index_bytes_for_test(std::span<const std::uint8_t> index,
                                         const Group& A, const Group& B);
  Scalar outer_gamma(const Scalar& x0, const Scalar& y0);
  const Digest& state() const { return state_; }

 private:
  void absorb(std::string_view label, std::span<const Bytes> fields);
  void absorb(std::string_view label, const Bytes& field);
  Scalar challenge_nonzero(std::string_view label);
  Scalar round_encoded(const Bytes& index, const Group& A, const Group& B);
  Digest state_{};
};

}

#include "pippenger/pippenger.hpp"

#include <mutex>

namespace pippenger {

void initialize() {
  static std::once_flag once;
  std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}

std::size_t scalar_bit_length() {
  initialize();
  return Scalar::getBitSize();
}

}  // namespace pippenger

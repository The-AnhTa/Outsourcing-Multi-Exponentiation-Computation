#include "pinkas/pinkas.hpp"

#include "internal/constants.hpp"
#include <string>

namespace pinkas {
namespace {

using internal::kChallengeBitLength;
using internal::kGroupIdentifier;
using internal::kProtocolVersion;

}

PublicParameters setup(std::size_t msm_window_width) {
  pippenger::initialize();

  static_cast<void>(pippenger::number_of_windows(msm_window_width));
  return {
      std::string(kGroupIdentifier),
      std::string(kProtocolVersion),
      Scalar::getBitSize(),
      kChallengeBitLength,
      msm_window_width};
}

}  

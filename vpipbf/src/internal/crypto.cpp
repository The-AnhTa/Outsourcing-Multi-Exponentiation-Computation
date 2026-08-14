#include "crypto.hpp"

#include "vpip_bf/group_utils.hpp"

#include <limits>
#include <stdexcept>

namespace vpip_bf::internal {

std::size_t dimension_size(std::size_t d) {
    if (d < 1 || d > 20 || d >= std::numeric_limits<std::size_t>::digits)
        throw std::invalid_argument("invalid dimension");
    return std::size_t{1} << d;
}

bool canonical_fr(const Fr& value) {
    try {
        const Bytes encoded = serialize(value);
        Fr decoded;
        return decoded.deserialize(encoded.data(), encoded.size())
                == encoded.size()
            && decoded == value && serialize(decoded) == encoded;
    } catch (...) {
        return false;
    }
}

bool valid_gt(const GT& value) {
    try {
        const Bytes encoded = serialize(value);
        GT decoded;
        return decoded.deserialize(encoded.data(), encoded.size())
                == encoded.size()
            && decoded == value && serialize(decoded) == encoded
            && mcl::bn::isValidGT(value);
    } catch (...) {
        return false;
    }
}

} // namespace vpip_bf::internal

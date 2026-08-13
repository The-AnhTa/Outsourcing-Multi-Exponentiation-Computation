#include "crypto.hpp"

#include "vme_ibf/group_utils.hpp"

#include <mcl/gmp_util.hpp>

#include <limits>
#include <stdexcept>

namespace vme_ibf::internal {

std::size_t dimension_size(std::size_t d) {
    if (d < 1 || d >= std::numeric_limits<std::size_t>::digits)
        throw std::invalid_argument("invalid dimension");
    return std::size_t{1} << d;
}

Fr fr_one() {
    Fr value;
    value = 1;
    return value;
}

Fr fr_minus_one() {
    Fr value;
    value = -1;
    return value;
}

Fr fr_multiply(const Fr& left, const Fr& right) {
    Fr result;
    Fr::mul(result, left, right);
    return result;
}

bool valid_gt(const GT& value) {
    try {
        const Bytes encoded = serialize(value);
        GT decoded;
        if (decoded.deserialize(encoded.data(), encoded.size()) != encoded.size()
            || decoded != value || serialize(decoded) != encoded)
            return false;
        mpz_class order;
        mcl::gmp::setStr(order, Fr::getModulo(), 10);
        GT powered;
        GT::pow(powered, value, order);
        GT one;
        one.setOne();
        return powered == one;
    } catch (...) {
        return false;
    }
}

double milliseconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}

} // namespace vme_ibf::internal

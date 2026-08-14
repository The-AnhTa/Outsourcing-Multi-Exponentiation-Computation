#include "crypto.hpp"

#include <limits>
#include <stdexcept>
#include <vector>

namespace rexpbf::internal {

void validate_dimension(std::size_t d) {
    if (d == 0) throw std::invalid_argument("d must be at least one");
    if (d >= std::numeric_limits<std::size_t>::digits) {
        throw std::invalid_argument("d overflows 1 << d");
    }
    const std::size_t n = std::size_t{1} << d;
    if (n > std::vector<G1>().max_size()
        || n > std::vector<G2>().max_size()) {
        throw std::invalid_argument("d exceeds vector maximum size");
    }
    if (n > (std::numeric_limits<std::uint64_t>::max() - 3) / 4) {
        throw std::invalid_argument("d overflows pairing count");
    }
}

Fr fr_one() { Fr result; result = 1; return result; }
Fr fr_multiply(const Fr& left, const Fr& right) {
    Fr result; Fr::mul(result, left, right); return result;
}
Fr fr_inverse_nonzero(const Fr& value) {
    if (value.isZero()) throw std::invalid_argument("cannot invert zero");
    Fr result; Fr::inv(result, value); return result;
}
G1 g1_multiply(const G1& point, const Fr& scalar) {
    G1 result; G1::mul(result, point, scalar); return result;
}
G2 g2_multiply(const G2& point, const Fr& scalar) {
    G2 result; G2::mul(result, point, scalar); return result;
}
G1 g1_add(const G1& left, const G1& right) {
    G1 result; G1::add(result, left, right); return result;
}
G2 g2_add(const G2& left, const G2& right) {
    G2 result; G2::add(result, left, right); return result;
}
GT gt_power(const GT& value, const Fr& scalar) {
    GT result; GT::pow(result, value, scalar); return result;
}
GT gt_multiply(const GT& left, const GT& right) {
    GT result; GT::mul(result, left, right); return result;
}
GT gt_product(std::initializer_list<GT> values) {
    GT result; result.setOne();
    for (const GT& value : values) result = gt_multiply(result, value);
    return result;
}

} 

#pragma once

#include "rexpbf/types.hpp"

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <string>

namespace rexpbf::internal {

void validate_dimension(std::size_t d);

template<class Point>
void require_point(const Point& point, bool nonzero, const char* name) {
    if (!point.isValid() || !point.isValidOrder()
        || (nonzero && point.isZero())) {
        throw std::invalid_argument(std::string("invalid ") + name);
    }
}

Fr fr_one();
Fr fr_multiply(const Fr& left, const Fr& right);
Fr fr_inverse_nonzero(const Fr& value);
G1 g1_multiply(const G1& point, const Fr& scalar);
G2 g2_multiply(const G2& point, const Fr& scalar);
G1 g1_add(const G1& left, const G1& right);
G2 g2_add(const G2& left, const G2& right);
GT gt_power(const GT& value, const Fr& scalar);
GT gt_multiply(const GT& left, const GT& right);
GT gt_product(std::initializer_list<GT> values);

} // namespace rexpbf::internal

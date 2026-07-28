#include "rexpbf/pairing.hpp"
#include <mutex>
#include <stdexcept>

namespace rexpbf {
void initialize_bn254() {
    static std::once_flag once;
    std::call_once(once, [] { mcl::bn::initPairing(mcl::BN254); });
}
GT pairing_product(std::span<const G1> left, std::span<const G2> right) {
    initialize_bn254();
    if (left.size() != right.size()) throw std::invalid_argument("pairing-product length mismatch");
    if (left.empty()) { GT one; one.setOne(); return one; }
    GT miller;
    mcl::bn::millerLoopVec(miller, left.data(), right.data(), left.size(), true);
    GT result;
    mcl::bn::finalExp(result, miller);
    return result;
}
}

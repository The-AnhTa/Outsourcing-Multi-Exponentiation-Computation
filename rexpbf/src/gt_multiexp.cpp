#include "rexpbf/gt_multiexp.hpp"
#include <stdexcept>
#include <vector>

namespace rexpbf {
GT gt_multiexp_pippenger(std::span<const GT> bases, std::span<const Fr> scalars) {
    if (bases.size()!=scalars.size()) throw std::invalid_argument("GT multiexp length mismatch");
    GT out; out.setOne();
    if (!bases.empty()) GT::powVec(out,bases.data(),scalars.data(),bases.size());
    return out;
}
GT gt_multiexp_naive(std::span<const GT> bases, std::span<const Fr> scalars) {
    if (bases.size()!=scalars.size()) throw std::invalid_argument("GT multiexp length mismatch");
    GT out;out.setOne();
    for(std::size_t i=0;i<bases.size();++i){GT p;GT::pow(p,bases[i],scalars[i]);GT::mul(out,out,p);}
    return out;
}
}

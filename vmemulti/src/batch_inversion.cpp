#include "vme_ibf/batch_inversion.hpp"
#include <stdexcept>
namespace vme_ibf {
std::vector<Fr> batch_invert_nonzero(std::span<const Fr>v){if(v.empty())return {};std::vector<Fr>prefix(v.size()),out(v.size());Fr acc;acc=1;for(size_t i=0;i<v.size();++i){if(v[i].isZero())throw std::invalid_argument("batch inversion contains zero");prefix[i]=acc;Fr::mul(acc,acc,v[i]);}Fr inv;Fr::inv(inv,acc);for(size_t i=v.size();i-->0;){Fr::mul(out[i],inv,prefix[i]);Fr::mul(inv,inv,v[i]);}return out;}
}

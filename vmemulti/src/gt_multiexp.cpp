#include "vme_ibf/gt_multiexp.hpp"
#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <string>
namespace vme_ibf { namespace {
GT one(){GT x;x.setOne();return x;} GT mul(const GT&a,const GT&b){GT x;GT::mul(x,a,b);return x;}
std::atomic_size_t window_override{};
unsigned nibble(char c){if(c>='0'&&c<='9')return c-'0';if(c>='a'&&c<='f')return c-'a'+10;if(c>='A'&&c<='F')return c-'A'+10;throw std::runtime_error("invalid scalar hex");}
unsigned digit(const std::string&h,size_t off,size_t width){unsigned d=0;for(size_t b=0;b<width;++b){size_t a=off+b,n=a/4;if(n<h.size())d|=((nibble(h[h.size()-1-n])>>(a%4))&1U)<<b;}return d;}
size_t window(size_t n){const auto override_value=window_override.load(std::memory_order_relaxed);if(override_value)return override_value;if(n<=4)return 2;if(n<=16)return 3;if(n<=160)return 4;if(n<=256)return 5;return 6;}
}}
namespace vme_ibf {
void set_gt_multiexp_window_override_for_benchmark(size_t w){if(w>8)throw std::invalid_argument("invalid GT MSM window override");window_override.store(w,std::memory_order_relaxed);}
GT gt_multiexp_reference(std::span<const GT>b,std::span<const Fr>e){if(b.size()!=e.size())throw std::invalid_argument("GT multiexp length mismatch");GT out=one();for(size_t i=0;i<b.size();++i){GT p;GT::pow(p,b[i],e[i]);out=mul(out,p);}return out;}
GT gt_multiexp_pippenger(std::span<const GT>b,std::span<const Fr>e){if(b.size()!=e.size())throw std::invalid_argument("GT multiexp length mismatch");if(b.empty())return one();size_t w=window(b.size()),bits=0;std::vector<std::string>hex;for(auto&s:e){hex.push_back(s.getStr(16));bits=std::max(bits,hex.back().size()*4);}size_t windows=(bits+w-1)/w,count=size_t{1}<<w;std::vector<GT>buckets(count);GT acc=one();for(size_t win=windows;win-->0;){for(size_t k=0;k<w;++k)acc=mul(acc,acc);for(auto&x:buckets)x=one();for(size_t i=0;i<b.size();++i){unsigned d=digit(hex[i],win*w,w);if(d)buckets[d]=mul(buckets[d],b[i]);}GT run=one();for(size_t k=count;k-->1;){run=mul(run,buckets[k]);acc=mul(acc,run);}}return acc;}
}

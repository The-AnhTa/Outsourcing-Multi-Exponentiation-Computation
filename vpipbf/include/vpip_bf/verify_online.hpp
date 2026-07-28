#pragma once
#include "vpip_bf/verify_combined.hpp"
#include <optional>
namespace vpip_bf {
struct OnlineTimingBreakdown;
class ValidatedVerificationInputs {
 public:
  const VpipBfCRS& crs() const { return *crs_; }
  const VpipBfPrecomputation& precomputation() const { return *precomp_; }
  const VpipBfStatement& statement() const { return *statement_; }
  const VpipBfProof& proof() const { return *proof_; }
 private:
  ValidatedVerificationInputs(const VpipBfCRS&c,const VpipBfPrecomputation&p,
      const VpipBfStatement&s,const VpipBfProof&v)
      :crs_(&c),precomp_(&p),statement_(&s),proof_(&v){}
  const VpipBfCRS* crs_;
  const VpipBfPrecomputation* precomp_;
  const VpipBfStatement* statement_;
  const VpipBfProof* proof_;
  friend std::optional<ValidatedVerificationInputs> PrevalidateVerificationInputs(
      const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&);
};
std::optional<ValidatedVerificationInputs> PrevalidateVerificationInputs(
    const VpipBfCRS&,const VpipBfPrecomputation&,const VpipBfStatement&,const VpipBfProof&);
bool verify_online(const ValidatedVerificationInputs&);
CombinedVerificationTrace verify_online_with_trace(const ValidatedVerificationInputs&);
bool verify_online_with_timing(const ValidatedVerificationInputs&,OnlineTimingBreakdown&);
}

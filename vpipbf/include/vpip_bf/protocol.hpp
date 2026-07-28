#pragma once
#include "vpip_bf/phase2.hpp"
#include "vpip_bf/proof.hpp"
#include "vpip_bf/verify_online.hpp"
#include <span>

namespace vpip_bf {
struct ProveResult { Statement statement; Proof proof; Phase1Result phase1; Phase2Result phase2; };
struct VerificationIoCallCounts {
  std::size_t full_proof_serialization_calls{};
  std::size_t full_proof_parse_calls{};
};
inline SetupResult Setup(std::size_t d,const std::vector<G2>&lambda,const std::vector<G1>&X,ScalarRng&rng){return setup_parameters(d,lambda,X,rng);}
inline Precomputation Precompute(const Crs& crs){return precompute(crs);}
ProveResult Prove(const Crs&,const Precomputation&,const VpipBfStatementInput&);
std::optional<ValidatedVerificationInputs> PrevalidateVerificationInputs(
    const Crs&,const Precomputation&,const Statement&,const Proof&);
bool Verify(const Crs&,const Precomputation&,const Statement&,const Proof&);
bool VerifySerialized(const Crs&,const Precomputation&,const Statement&,
    std::span<const std::uint8_t> proof_bytes);
bool VerifyWithIoCountsForTesting(const Crs&,const Precomputation&,const Statement&,
    const Proof&,VerificationIoCallCounts&);
bool VerifySerializedWithIoCountsForTesting(const Crs&,const Precomputation&,
    const Statement&,std::span<const std::uint8_t> proof_bytes,VerificationIoCallCounts&);
bool VerifyPrevalidated(const ValidatedVerificationInputs&);
bool VerifyOnline(const ValidatedVerificationInputs&);
bool VerifyOnline(const ValidatedVerificationInputs&,OnlineTimingBreakdown&);
}

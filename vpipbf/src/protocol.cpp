#include "vpip_bf/protocol.hpp"
#include "vpip_bf/serialization.hpp"
namespace vpip_bf {
namespace {
bool verify_parsed(const Crs&c,const Precomputation&p,const Statement&s,const Proof&v,
    VerificationIoCallCounts*){
  auto inputs=PrevalidateVerificationInputs(c,p,s,v);
  return inputs&&verify_online(*inputs);
}
bool verify_serialized(const Crs&c,const Precomputation&p,const Statement&s,
    std::span<const std::uint8_t>proof_bytes,VerificationIoCallCounts*counts){
  Proof proof;DecodeError error;if(counts)++counts->full_proof_parse_calls;
  if(!deserialize_proof(proof_bytes,c,proof,&error))return false;
  return verify_parsed(c,p,s,proof,counts);
}
}
ProveResult Prove(const Crs&c,const Precomputation&p,const VpipBfStatementInput&i){ProveResult z;z.phase1=prove_phase1(c,p,i);z.phase2=prove_phase2(c,p,z.phase1);z.statement=z.phase1.statement;z.proof=assemble_public_proof(z.phase1,z.phase2);return z;}
bool Verify(const Crs&c,const Precomputation&p,const Statement&s,const Proof&v){try{return verify_parsed(c,p,s,v,nullptr);}catch(...){return false;}}
bool VerifySerialized(const Crs&c,const Precomputation&p,const Statement&s,std::span<const std::uint8_t>b){try{return verify_serialized(c,p,s,b,nullptr);}catch(...){return false;}}
bool VerifyWithIoCountsForTesting(const Crs&c,const Precomputation&p,const Statement&s,const Proof&v,VerificationIoCallCounts&counts){counts={};try{return verify_parsed(c,p,s,v,&counts);}catch(...){return false;}}
bool VerifySerializedWithIoCountsForTesting(const Crs&c,const Precomputation&p,const Statement&s,std::span<const std::uint8_t>b,VerificationIoCallCounts&counts){counts={};try{return verify_serialized(c,p,s,b,&counts);}catch(...){return false;}}
bool VerifyPrevalidated(const ValidatedVerificationInputs&i){return verify_online(i);}
bool VerifyOnline(const ValidatedVerificationInputs&i){return verify_online(i);}
bool VerifyOnline(const ValidatedVerificationInputs&i,OnlineTimingBreakdown&t){return verify_online_with_timing(i,t);}
}

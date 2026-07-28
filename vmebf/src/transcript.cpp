#include "vme_ibf/transcript.hpp"
#include <algorithm>
#include <stdexcept>

namespace vme_ibf {


static constexpr Digest kRejectionLimit = {0xde,0xd4,0x5b,0x0d,0x80,0x00,0x00,0x0a,0x5d,0x39,0xd1,0x00,0x00,0x00,0x00,0x2f,0xfd,0xbd,0x00,0x00,0x00,0x00,0x00,0x63,0xc6,0x00,0x00,0x00,0x00,0x00,0x00,0x4e};
Transcript::Transcript(std::span<const std::uint8_t> d){Bytes b;append_frame(b,"VME.BF.G2/TRANSCRIPT/V2");append_frame(b,d);state_=sha256(b);}
Transcript Transcript::resume(const Digest& state){return Transcript(state,ResumeTag{});}
void Transcript::absorb(std::string_view label,std::span<const Bytes> fields){Bytes b;append_frame(b,"VME.BF.G2/ABSORB/V2");append_frame(b,label);append_frame(b,state_);for(const auto&f:fields)append_frame(b,f);state_=sha256(b);}
void Transcript::absorb(std::string_view label,const Bytes& field){absorb(label,std::span<const Bytes>(&field,1));}
Fr Transcript::challenge_nonzero(std::string_view label,std::uint64_t index){
  for(std::uint64_t counter=0;;++counter){Bytes b;append_frame(b,"VME.BF.G2/CHALLENGE/V2");append_frame(b,label);append_frame(b,state_);auto i=encode_u64_be(index),c=encode_u64_be(counter);b.insert(b.end(),i.begin(),i.end());b.insert(b.end(),c.begin(),c.end());Digest h=sha256(b);if(!std::lexicographical_compare(h.begin(),h.end(),kRejectionLimit.begin(),kRejectionLimit.end()))continue;Fr out;out.setBigEndianMod(h.data(),h.size());if(out.isZero())continue;Bytes encoded=serialize(out);std::string value_label(label);value_label+="/value";absorb(value_label,encoded);return out;}
}
}

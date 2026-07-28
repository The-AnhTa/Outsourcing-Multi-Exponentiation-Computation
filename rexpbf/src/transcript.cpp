#include "rexpbf/transcript.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <chrono>

namespace rexpbf {
namespace {
std::span<const std::uint8_t> bytes(std::string_view s) {
    return {reinterpret_cast<const std::uint8_t*>(s.data()), s.size()};
}
bool acceptable(const Digest32& x) {

    static constexpr Digest32 limit = {
        0xf1,0xf5,0x88,0x3e,0x65,0xf8,0x20,0xd0,0x99,0x91,0x5c,0x90,0x87,0x86,0xb9,0xd1,
        0xc9,0x03,0x89,0x6a,0x60,0x9f,0x32,0xd6,0x53,0x69,0xcb,0xe3,0xb0,0x00,0x00,0x05};
    return std::lexicographical_compare(x.begin(), x.end(), limit.begin(), limit.end());
}
}
Transcript::Transcript(std::string_view domain, TranscriptMetrics* metrics) : metrics_(metrics) {
    if (metrics_) ++metrics_->transcript_initializations;
    append_entry("REXP-BF-G1-TRANSCRIPT-INIT-V1", "protocol-domain", bytes(domain));
}
void Transcript::append_entry(std::string_view domain, std::string_view label,
                              std::span<const std::uint8_t> value) {
    std::chrono::steady_clock::time_point begin;
    if(metrics_)begin=std::chrono::steady_clock::now();std::vector<std::uint8_t> entry;
    entry.reserve(domain.size()+label.size()+value.size()+24);
    append_frame(entry, bytes(domain)); append_frame(entry, bytes(label)); append_frame(entry, value);
    std::chrono::steady_clock::time_point framed;
    if(metrics_)framed=std::chrono::steady_clock::now();hash_.update(entry.data(),entry.size());
    if(metrics_){auto updated=std::chrono::steady_clock::now();++metrics_->transcript_entries;metrics_->bytes_absorbed+=entry.size();metrics_->serialization_ms+=std::chrono::duration<double,std::milli>(framed-begin).count();metrics_->sha256_ms+=std::chrono::duration<double,std::milli>(updated-framed).count();}
}
void Transcript::append_u64(std::string_view label, std::uint64_t value) {
    const auto b = encode_u64_be(value); append_entry("REXP-BF-G1-TRANSCRIPT-U64-V1", label, b);
}
void Transcript::append_bytes(std::string_view label, std::span<const std::uint8_t> value) {
    append_entry("REXP-BF-G1-TRANSCRIPT-BYTES-V1", label, value);
}
void Transcript::append_g1(std::string_view label, const G1& p) {
    std::chrono::steady_clock::time_point t;if(metrics_)t=std::chrono::steady_clock::now();
    const auto b = serialize_g1(p); append_entry("REXP-BF-G1-TRANSCRIPT-G1-V1", label, b);
    if(metrics_){metrics_->serialization_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();++metrics_->g1_serializations;}
}
void Transcript::append_g2(std::string_view label, const G2& p) {
    std::chrono::steady_clock::time_point t;if(metrics_)t=std::chrono::steady_clock::now();
    const auto b = serialize_g2(p); append_entry("REXP-BF-G1-TRANSCRIPT-G2-V1", label, b);
    if(metrics_){metrics_->serialization_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();++metrics_->g2_serializations;}
}
void Transcript::append_gt(std::string_view label, const GT& p) {
    std::chrono::steady_clock::time_point t;if(metrics_)t=std::chrono::steady_clock::now();
    const auto b = serialize_gt(p); append_entry("REXP-BF-G1-TRANSCRIPT-GT-V1", label, b);
    if(metrics_){metrics_->serialization_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();++metrics_->gt_serializations;}
}
Digest32 Transcript::digest() const {std::chrono::steady_clock::time_point t;if(metrics_)t=std::chrono::steady_clock::now();cybozu::Sha256 copy=hash_;Digest32 d{};copy.digest(d.data(),d.size(),nullptr,0);if(metrics_){metrics_->sha256_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t).count();++metrics_->sha256_calls;}return d;}
Fr Transcript::challenge_nonzero_fr(std::string_view label, std::uint64_t index) {
    if(metrics_)++metrics_->challenge_derivations;
    for (std::uint64_t counter = 0;; ++counter) {
        std::vector<std::uint8_t> input;
        append_frame(input, bytes("REXP-BF-G1-CHALLENGE-NONZERO-FR-V1"));
        append_frame(input, bytes(label));
        append_frame(input, digest());
        const auto i = encode_u64_be(index), c = encode_u64_be(counter);
        input.insert(input.end(), i.begin(), i.end()); input.insert(input.end(), c.begin(), c.end());
        std::chrono::steady_clock::time_point ht;if(metrics_)ht=std::chrono::steady_clock::now();const auto candidate = sha256(input);if(metrics_){metrics_->sha256_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-ht).count();++metrics_->sha256_calls;}
        if (acceptable(candidate)) {
            std::chrono::steady_clock::time_point ft;if(metrics_)ft=std::chrono::steady_clock::now();Fr value; value.setBigEndianMod(candidate.data(), candidate.size());if(metrics_)metrics_->challenge_to_field_ms+=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-ft).count();
            if (!value.isZero()) {
                const auto encoded = serialize_fr(value);
                append_entry("REXP-BF-G1-CHALLENGE-RESULT-V1", label, encoded);
                return value;
            }
        }
        if(metrics_)++metrics_->rejection_sampling_retries;
        if (counter == std::numeric_limits<std::uint64_t>::max())
            throw std::runtime_error("Fiat-Shamir rejection counter exhausted");
    }
}
}

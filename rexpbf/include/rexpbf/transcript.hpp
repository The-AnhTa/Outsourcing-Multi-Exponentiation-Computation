#pragma once
#include "rexpbf/serialization.hpp"
#include <cybozu/sha2.hpp>
#include <string_view>

namespace rexpbf {
struct TranscriptMetrics {
    double serialization_ms{}, sha256_ms{}, challenge_to_field_ms{};
    std::size_t g1_serializations{}, g2_serializations{}, gt_serializations{};
    std::size_t transcript_initializations{}, transcript_entries{}, bytes_absorbed{}, sha256_calls{},
        challenge_derivations{}, rejection_sampling_retries{};
};
class Transcript {
public:
    explicit Transcript(std::string_view protocol_domain, TranscriptMetrics* metrics = nullptr);
    void append_u64(std::string_view label, std::uint64_t value);
    void append_bytes(std::string_view label, std::span<const std::uint8_t> bytes);
    void append_g1(std::string_view label, const G1& point);
    void append_g2(std::string_view label, const G2& point);
    void append_gt(std::string_view label, const GT& value);
    Fr challenge_nonzero_fr(std::string_view label, std::uint64_t logical_index);
    Digest32 digest() const;
private:
    cybozu::Sha256 hash_;
    TranscriptMetrics* metrics_{};
    void append_entry(std::string_view domain, std::string_view label,
                      std::span<const std::uint8_t> value);
};
}

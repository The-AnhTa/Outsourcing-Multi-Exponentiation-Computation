#pragma once

#include <mcl/bn.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <optional>
#include <vector>

namespace bp {

using Scalar = mcl::bn::Fr;
using Group = mcl::bn::G2;
using Bytes = std::vector<std::uint8_t>;
using Digest = std::array<std::uint8_t, 32>;

inline constexpr const char* kGroupIdentifier = "BN254/G2/mcl-v3.00";
inline constexpr const char* kScalarModulus =
    "0x2523648240000001ba344d8000000007ff9f800000000010a10000000000000d";
inline constexpr const char* kHashSuiteIdentifier = "SHA-256+mcl-hashAndMapToG2";
inline constexpr const char* kTranscriptDomain = "BP-IPA-FS-v1";
inline constexpr const char* kHpTranscriptDomain = "BPVME/HP/BP/v1";

struct PublicParams {
  std::size_t n{};
  std::size_t d{};
  std::vector<Group> G;
  std::vector<Group> H;
  Group K;
  std::string group_identifier{kGroupIdentifier};
  std::string scalar_modulus{kScalarModulus};
  std::string hash_suite_identifier{kHashSuiteIdentifier};
  std::string transcript_domain{kTranscriptDomain};
  std::optional<Digest> transcript_crs_digest;
};

struct RoundMessage {
  Group A;
  Group B;
};

struct Proof {
  std::vector<RoundMessage> rounds;
  Scalar x_final;
  Scalar y_final;
};

struct Trace {
  std::vector<Scalar> challenges;
  std::vector<bool> round_invariants;
};

void initialize();
PublicParams Setup(std::size_t n, std::span<const std::uint8_t> public_seed);
Proof Prove(const PublicParams& pp, const Group& Z,
            std::span<const Scalar> x, std::span<const Scalar> y,
            Trace* trace = nullptr);
bool Verify(const PublicParams& pp, const Group& Z, const Proof& proof,
            Trace* trace = nullptr) noexcept;

bool validate_public_params(const PublicParams& pp) noexcept;
Group commit(const PublicParams& pp, std::span<const Scalar> x,
             std::span<const Scalar> y);

Bytes serialize_proof(const PublicParams& pp, const Proof& proof);
bool deserialize_proof(const PublicParams& pp, std::span<const std::uint8_t> bytes,
                       Proof& proof) noexcept;
bool VerifySerialized(const PublicParams& pp, const Group& Z,
                      std::span<const std::uint8_t> bytes) noexcept;

Bytes serialize(const Scalar& value);
Bytes serialize(const Group& value);
Digest sha256(std::span<const std::uint8_t> input);
std::size_t scalar_bytes();
std::size_t group_bytes();
std::size_t proof_payload_bytes(std::size_t n);
std::size_t proof_wire_bytes(std::size_t n);

}

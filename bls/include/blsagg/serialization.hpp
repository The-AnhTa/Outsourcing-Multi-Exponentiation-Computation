#pragma once
#include "blsagg/protocol.hpp"

namespace blsagg {

enum class DecodeError {
  None,
  Truncated,
  TrailingBytes,
  WrongMagic,
  WrongVersion,
  WrongMode,
  InvalidDimension,
  IntegerOverflow,
  InvalidLength,
  InvalidG1,
  InvalidG2,
  InvalidGT,
  NonCanonical,
  IdentityNotAllowed,
  InvalidDigest,
  InvalidShape
};

class ValidatedProof {
 public:
  ValidatedProof() = delete;
  ValidatedProof(const ValidatedProof&) = default;
  ValidatedProof(ValidatedProof&&) = default;
  ValidatedProof& operator=(const ValidatedProof&) = delete;
  ValidatedProof& operator=(ValidatedProof&&) = delete;
  const Proof& proof() const { return proof_; }
  const Digest& parameter_digest() const { return parameter_digest_; }
  const Digest& wire_binding() const { return wire_binding_; }
 private:
  ValidatedProof(Proof, Digest, Digest);
  const Proof proof_;
  const Digest parameter_digest_;
  const Digest wire_binding_;
  friend std::optional<ValidatedProof> deserialize_and_validate_proof(
      std::span<const std::uint8_t>, const PublicParameters&, DecodeError*);
};

Bytes serialize_public_parameters(const PublicParameters&);
bool deserialize_public_parameters(std::span<const std::uint8_t>,
                                   PublicParameters&, DecodeError* = nullptr);

Bytes serialize_precomputation(const PublicParameters&, const Precomputation&);
bool deserialize_precomputation(std::span<const std::uint8_t>,
                                const PublicParameters&, Precomputation&,
                                DecodeError* = nullptr);

Bytes serialize_statement(const PublicParameters&, const Statement&);
bool deserialize_statement(std::span<const std::uint8_t>,
                           const PublicParameters&, Statement&,
                           DecodeError* = nullptr);

Bytes serialize_proof(const PublicParameters&, const Proof&);
bool deserialize_proof(std::span<const std::uint8_t>,
                       const PublicParameters&, Proof&,
                       DecodeError* = nullptr);
std::optional<ValidatedProof> deserialize_and_validate_proof(
    std::span<const std::uint8_t>, const PublicParameters&,
    DecodeError* = nullptr);

std::size_t proof_mathematical_payload_bytes(const Proof&);
std::size_t proof_wire_bytes(const PublicParameters&, const Proof&);

}

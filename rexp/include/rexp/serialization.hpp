#pragma once

#include "rexp/metrics.hpp"
#include "rexp/proof.hpp"

#include <cstdint>
#include <vector>

namespace rexp {

std::vector<std::uint8_t> SerializeRexpCRS(const RawRexpCRS&);
std::vector<std::uint8_t> SerializeRexpPrecomputation(
    const PreparedPublicParameters&);
std::vector<std::uint8_t> SerializeRexpStatement(const RawRexpStatement&);
std::vector<std::uint8_t> SerializePreparedStatement(const PreparedStatement&);
std::vector<std::uint8_t> SerializeRexpProof(const RexpProof&, std::size_t d);
std::vector<std::uint8_t> SerializeRexpCRSWire(const RawRexpCRS&);
std::vector<std::uint8_t> SerializeRexpStatementWire(
    const RawRexpStatement&, std::size_t n);
std::vector<std::uint8_t> SerializeRexpProofWire(
    const RexpProof&, std::size_t d);
RawRexpCRS DeserializeRexpCRSWire(const std::vector<std::uint8_t>&);
RawRexpStatement DeserializeRexpStatementWire(
    const std::vector<std::uint8_t>&, std::size_t expected_n);
RexpProof DeserializeRexpProofWire(
    const std::vector<std::uint8_t>&, std::size_t expected_d,
    RexpProofValidationMetrics* metrics = nullptr);
ValidatedRexpProof DeserializeValidatedRexpProofWire(
    const std::vector<std::uint8_t>&, std::size_t expected_d,
    RexpProofValidationMetrics* metrics = nullptr);

} // namespace rexp

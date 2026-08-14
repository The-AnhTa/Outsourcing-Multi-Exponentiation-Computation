#include "rexp/dory.hpp"

#include "internal/crypto.hpp"

#include <stdexcept>

namespace rexp {

std::vector<std::uint8_t> SerializeProof(
    const DoryProof& proof, std::size_t d) {
    if (proof.rounds.size() != d) {
        throw std::invalid_argument("proof round count differs from d");
    }
    internal::Bytes out;
    internal::frame(out, "DORY-PROOF-WIRE-BN254-V1");
    internal::frame(out, "BN254");
    internal::append_u64(out, d);
    internal::append_u64(out, proof.rounds.size());
    for (const DoryRound& round : proof.rounds) {
        internal::frame_element(out, round.D1L);
        internal::frame_element(out, round.D1R);
        internal::frame_element(out, round.D2L);
        internal::frame_element(out, round.D2R);
        internal::frame_element(out, round.W1);
        internal::frame_element(out, round.W2);
    }
    internal::frame_element(out, proof.PhiFinal);
    internal::frame_element(out, proof.ThetaFinal);
    return out;
}

std::vector<std::uint8_t> SerializeBatchProof(
    const DoryBatchProof& proof, std::size_t d, std::size_t batch_size) {
    if (batch_size == 0 || proof.batchCrossTerms.size() != batch_size - 1) {
        throw std::invalid_argument("batch proof cross-term count mismatch");
    }
    internal::Bytes out;
    internal::frame(out, "DORY-BATCH-PROOF-WIRE-BN254-V1");
    internal::frame(out, "BN254");
    internal::append_u64(out, d);
    internal::append_u64(out, batch_size);
    internal::append_u64(out, proof.batchCrossTerms.size());
    for (const GT& cross : proof.batchCrossTerms) {
        internal::frame_element(out, cross);
    }
    const internal::Bytes ordinary = SerializeProof(proof.doryProof, d);
    internal::frame(out, ordinary.data(), ordinary.size());
    return out;
}

ProofSizes MeasureProofSizes(const DoryProof& proof, std::size_t d) {
    if (proof.rounds.size() != d) {
        throw std::invalid_argument("proof round count differs from d");
    }
    std::size_t payload = 0;
    for (const DoryRound& round : proof.rounds) {
        payload += internal::encode(round.D1L).size();
        payload += internal::encode(round.D1R).size();
        payload += internal::encode(round.D2L).size();
        payload += internal::encode(round.D2R).size();
        payload += internal::encode(round.W1).size();
        payload += internal::encode(round.W2).size();
    }
    payload += internal::encode(proof.PhiFinal).size();
    payload += internal::encode(proof.ThetaFinal).size();
    return {payload, SerializeProof(proof, d).size()};
}

ProofSizes MeasureProofSizes(
    const DoryBatchProof& proof, std::size_t d, std::size_t batch_size) {
    if (batch_size == 0 || proof.batchCrossTerms.size() != batch_size - 1) {
        throw std::invalid_argument("batch proof cross-term count mismatch");
    }
    std::size_t payload =
        MeasureProofSizes(proof.doryProof, d).mathematical_payload_bytes;
    for (const GT& cross : proof.batchCrossTerms) {
        payload += internal::encode(cross).size();
    }
    return {payload, SerializeBatchProof(proof, d, batch_size).size()};
}

} 

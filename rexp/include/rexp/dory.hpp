#pragma once

#include "rexp/dory_setup.hpp"

#include <cstddef>
#include <string_view>
#include <vector>

namespace rexp {

struct DoryRound {
    GT D1L;
    GT D1R;
    GT D2L;
    GT D2R;
    GT W1;
    GT W2;
};

struct DoryProof {
    std::vector<DoryRound> rounds;
    G1 PhiFinal;
    G2 ThetaFinal;
};

struct DoryBatchProof {
    std::vector<GT> batchCrossTerms;
    DoryProof doryProof;
};

struct DoryChallenges {
    std::vector<Fr> beta;
    std::vector<Fr> alpha;
    Fr epsilon;
};

struct DoryBatchChallenges {
    std::vector<Fr> gamma;
    DoryChallenges dory;
};


Fr ChallengeNonzeroFr(
    const Digest& transcript_digest,
    std::string_view label,
    std::size_t index);

DoryProof Prove(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryWitness& witness);


DoryProof ProveEmbedded(
    const DoryCRS& crs,
    const DoryWitness& witness,
    const Digest& transcript_in,
    Digest* transcript_end = nullptr);

DoryBatchProof ProveBatch(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const std::vector<DoryWitness>& witnesses);


DoryChallenges DeriveChallenges(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryProof& proof);

std::vector<Fr> DeriveBatchGammas(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof);

DoryBatchChallenges DeriveBatchChallenges(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof);


std::vector<std::uint8_t> SerializeProof(
    const DoryProof& proof,
    std::size_t d);
std::vector<std::uint8_t> SerializeBatchProof(
    const DoryBatchProof& proof,
    std::size_t d,
    std::size_t batch_size);

struct ProofSizes {
    std::size_t mathematical_payload_bytes = 0;
    std::size_t wire_bytes = 0;
};



ProofSizes MeasureProofSizes(const DoryProof& proof, std::size_t d);
ProofSizes MeasureProofSizes(
    const DoryBatchProof& proof,
    std::size_t d,
    std::size_t batch_size);

}

#include "dory_transcript.hpp"

#include "crypto.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace rexp {
namespace internal {

Digest dory_initial(const DoryCRS& crs, const DoryStatement& statement) {
    Bytes input;
    frame(input, "DORY-NI-NONZK-V1");
    frame_digest(input, crs.digest);
    append_u64(input, crs.d);
    append_u64(input, crs.n);
    frame_element(input, statement.D0);
    frame_element(input, statement.D1);
    frame_element(input, statement.D2);
    return sha256(input);
}

Digest dory_batch_initial(
    const DoryCRS& crs, const std::vector<DoryStatement>& statements) {
    Bytes input;
    frame(input, "DORY-BATCH-NI-NONZK-V1");
    frame_digest(input, crs.digest);
    append_u64(input, crs.d);
    append_u64(input, crs.n);
    append_u64(input, statements.size());
    for (std::size_t j = 0; j < statements.size(); ++j) {
        append_u64(input, j);
        frame_element(input, statements[j].D0);
        frame_element(input, statements[j].D1);
        frame_element(input, statements[j].D2);
    }
    return sha256(input);
}

Digest dory_absorb_beta(
    const Digest& current, std::size_t round, std::size_t dimension,
    const DoryRound& message) {
    Bytes input;
    frame(input, "DORY-FOLD-BETA-MESSAGE-V1");
    frame_digest(input, current);
    append_u64(input, round);
    append_u64(input, dimension);
    frame_element(input, message.D1L);
    frame_element(input, message.D1R);
    frame_element(input, message.D2L);
    frame_element(input, message.D2R);
    return sha256(input);
}

Digest dory_absorb_alpha(
    const Digest& beta_digest, std::size_t round, std::size_t dimension,
    const DoryRound& message) {
    Bytes input;
    frame(input, "DORY-FOLD-ALPHA-MESSAGE-V1");
    frame_digest(input, beta_digest);
    append_u64(input, round);
    append_u64(input, dimension);
    frame_element(input, message.W1);
    frame_element(input, message.W2);
    return sha256(input);
}

Digest dory_absorb_final(
    const Digest& current, const G1& phi, const G2& theta) {
    Bytes input;
    frame(input, "DORY-FINAL-CHECK-V1");
    frame_digest(input, current);
    frame_element(input, phi);
    frame_element(input, theta);
    return sha256(input);
}

Digest dory_absorb_merge(
    const Digest& current, std::size_t instance, const GT& cross_term) {
    Bytes input;
    frame(input, "DORY-BATCH-MERGE-V1");
    frame_digest(input, current);
    append_u64(input, instance);
    frame_element(input, cross_term);
    return sha256(input);
}

Digest dory_enter_batch(const Digest& current) {
    Bytes input;
    frame(input, "DORY-BATCH-AGGREGATE-SUBPROTOCOL-V1");
    frame_digest(input, current);
    return sha256(input);
}

DoryTranscriptReplay replay_dory_transcript(
    const DoryCRS& crs, const DoryProof& proof, Digest transcript) {
    if (proof.rounds.size() != crs.d) {
        throw std::invalid_argument("Dory proof round count differs from CRS");
    }
    DoryTranscriptReplay out;
    out.challenges.beta.reserve(crs.d);
    out.challenges.alpha.reserve(crs.d);
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t dimension = crs.n >> k;
        const Digest beta_digest =
            dory_absorb_beta(transcript, k, dimension, proof.rounds[k]);
        out.challenges.beta.push_back(
            ChallengeNonzeroFr(beta_digest, "DORY-BETA-V1", k));
        transcript =
            dory_absorb_alpha(beta_digest, k, dimension, proof.rounds[k]);
        out.challenges.alpha.push_back(
            ChallengeNonzeroFr(transcript, "DORY-ALPHA-V1", k));
    }
    out.end = dory_absorb_final(
        transcript, proof.PhiFinal, proof.ThetaFinal);
    out.challenges.epsilon =
        ChallengeNonzeroFr(out.end, "DORY-EPSILON-V1", crs.d);
    return out;
}

namespace {

bool digest_less_than_limit(const Digest& value) {
    static constexpr Digest limit = {
        0xf1,0xf5,0x88,0x3e,0x65,0xf8,0x20,0xd0,
        0x99,0x91,0x5c,0x90,0x87,0x86,0xb9,0xd1,
        0xc9,0x03,0x89,0x6a,0x60,0x9f,0x32,0xd6,
        0x53,0x69,0xcb,0xe3,0xb0,0x00,0x00,0x05};
    return std::lexicographical_compare(
        value.begin(), value.end(), limit.begin(), limit.end());
}

} // namespace
} // namespace internal

Fr ChallengeNonzeroFr(
    const Digest& transcript_digest,
    std::string_view label,
    std::size_t index) {
    initialize();
    for (std::uint64_t counter = 0;; ++counter) {
        internal::Bytes input;
        internal::frame(input, "FS-NONZERO-FR-V1");
        internal::frame(input, label);
        internal::frame_digest(input, transcript_digest);
        internal::append_u64(input, index);
        internal::append_u64(input, counter);
        const Digest candidate = internal::sha256(input);
        if (!internal::digest_less_than_limit(candidate)) continue;
        Fr value;
        value.setBigEndianMod(candidate.data(), candidate.size());
        if (!value.isZero()) return value;
        if (counter == std::numeric_limits<std::uint64_t>::max()) {
            throw std::runtime_error("Fiat-Shamir rejection counter exhausted");
        }
    }
}

} // namespace rexp

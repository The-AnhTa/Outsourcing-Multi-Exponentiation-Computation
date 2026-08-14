#include "rexp/dory.hpp"

#include "internal/crypto.hpp"
#include "internal/dory_transcript.hpp"

#include <stdexcept>
#include <utility>

namespace rexp {
namespace {

using internal::dory_absorb_alpha;
using internal::dory_absorb_beta;
using internal::dory_absorb_final;
using internal::dory_absorb_merge;
using internal::dory_batch_initial;
using internal::dory_enter_batch;
using internal::dory_initial;
using internal::fr_mul;
using internal::g1_add;
using internal::g1_mul;
using internal::g2_add;
using internal::g2_mul;
using internal::gt_mul;
using internal::gt_pow;
using internal::inverse;
using internal::pairing_product;
using internal::replay_dory_transcript;

void validate_witness(const DoryCRS& crs, const DoryWitness& witness) {
    if (!ValidateCRS(crs)) throw std::invalid_argument("invalid CRS");
    if (witness.Phi.size() != crs.n || witness.Theta.size() != crs.n) {
        throw std::invalid_argument("witness vector length differs from CRS");
    }
    for (const G1& point : witness.Phi) {
        if (!point.isValid() || !point.isValidOrder()) {
            throw std::invalid_argument("invalid Phi point");
        }
    }
    for (const G2& point : witness.Theta) {
        if (!point.isValid() || !point.isValidOrder()) {
            throw std::invalid_argument("invalid Theta point");
        }
    }
}

struct FoldResult {
    DoryProof proof;
    DoryChallenges challenges;
    Digest transcript_end{};
};

FoldResult prove_from_transcript(
    const DoryCRS& crs,
    const DoryWitness& witness,
    Digest transcript) {
    std::vector<G1> phi = witness.Phi;
    std::vector<G2> theta = witness.Theta;
    FoldResult result;
    result.proof.rounds.reserve(crs.d);
    result.challenges.beta.reserve(crs.d);
    result.challenges.alpha.reserve(crs.d);
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        const std::size_t h = m / 2;
        DoryRound round;
        round.D1L = pairing_product(phi, 0, crs.Lambda, 0, h);
        round.D1R = pairing_product(phi, h, crs.Lambda, 0, h);
        round.D2L = pairing_product(crs.Gamma, 0, theta, 0, h);
        round.D2R = pairing_product(crs.Gamma, 0, theta, h, h);

        const Digest beta_digest = dory_absorb_beta(transcript, k, m, round);
        const Fr beta = ChallengeNonzeroFr(beta_digest, "DORY-BETA-V1", k);
        const Fr beta_inv = inverse(beta);
        std::vector<G1> phi_circle(m);
        std::vector<G2> theta_circle(m);
        for (std::size_t i = 0; i < m; ++i) {
            phi_circle[i] = g1_add(phi[i], g1_mul(crs.Gamma[i], beta));
            theta_circle[i] =
                g2_add(theta[i], g2_mul(crs.Lambda[i], beta_inv));
        }
        round.W1 = pairing_product(phi_circle, 0, theta_circle, h, h);
        round.W2 = pairing_product(phi_circle, h, theta_circle, 0, h);

        const Digest alpha_digest =
            dory_absorb_alpha(beta_digest, k, m, round);
        const Fr alpha =
            ChallengeNonzeroFr(alpha_digest, "DORY-ALPHA-V1", k);
        const Fr alpha_inv = inverse(alpha);
        std::vector<G1> phi_next(h);
        std::vector<G2> theta_next(h);
        for (std::size_t i = 0; i < h; ++i) {
            phi_next[i] =
                g1_add(g1_mul(phi_circle[i], alpha), phi_circle[h + i]);
            theta_next[i] =
                g2_add(g2_mul(theta_circle[i], alpha_inv), theta_circle[h + i]);
        }
        result.proof.rounds.push_back(std::move(round));
        result.challenges.beta.push_back(beta);
        result.challenges.alpha.push_back(alpha);
        phi = std::move(phi_next);
        theta = std::move(theta_next);
        transcript = alpha_digest;
    }
    result.proof.PhiFinal = phi.front();
    result.proof.ThetaFinal = theta.front();
    result.transcript_end = dory_absorb_final(
        transcript, result.proof.PhiFinal, result.proof.ThetaFinal);
    result.challenges.epsilon = ChallengeNonzeroFr(
        result.transcript_end, "DORY-EPSILON-V1", crs.d);
    return result;
}

} 

DoryProof Prove(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryWitness& witness) {
    validate_witness(crs, witness);
    return prove_from_transcript(
        crs, witness, dory_initial(crs, statement)).proof;
}

DoryProof ProveEmbedded(
    const DoryCRS& crs,
    const DoryWitness& witness,
    const Digest& transcript_in,
    Digest* transcript_end) {
    validate_witness(crs, witness);
    FoldResult result = prove_from_transcript(crs, witness, transcript_in);
    if (transcript_end) *transcript_end = result.transcript_end;
    return result.proof;
}

DoryBatchProof ProveBatch(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const std::vector<DoryWitness>& witnesses) {
    if (!ValidateCRS(crs)) throw std::invalid_argument("invalid CRS");
    if (statements.empty() || statements.size() != witnesses.size()) {
        throw std::invalid_argument("batch statements/witnesses mismatch or empty");
    }
    for (const DoryWitness& witness : witnesses) validate_witness(crs, witness);

    Digest transcript = dory_batch_initial(crs, statements);
    DoryStatement statement_acc = statements.front();
    DoryWitness witness_acc = witnesses.front();
    DoryBatchProof result;
    result.batchCrossTerms.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        const GT left = pairing_product(
            witness_acc.Phi, 0, witnesses[j].Theta, 0, crs.n);
        const GT right = pairing_product(
            witnesses[j].Phi, 0, witness_acc.Theta, 0, crs.n);
        const GT cross = gt_mul(left, right);
        result.batchCrossTerms.push_back(cross);
        transcript = dory_absorb_merge(transcript, j, cross);
        const Fr gamma =
            ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j);
        const Fr gamma_sq = fr_mul(gamma, gamma);
        statement_acc.D0 = gt_mul(
            gt_mul(gt_pow(statement_acc.D0, gamma_sq), gt_pow(cross, gamma)),
            statements[j].D0);
        statement_acc.D1 =
            gt_mul(gt_pow(statement_acc.D1, gamma), statements[j].D1);
        statement_acc.D2 =
            gt_mul(gt_pow(statement_acc.D2, gamma), statements[j].D2);
        for (std::size_t i = 0; i < crs.n; ++i) {
            witness_acc.Phi[i] = g1_add(
                g1_mul(witness_acc.Phi[i], gamma), witnesses[j].Phi[i]);
            witness_acc.Theta[i] = g2_add(
                g2_mul(witness_acc.Theta[i], gamma), witnesses[j].Theta[i]);
        }
    }
    transcript = dory_enter_batch(transcript);
    result.doryProof = prove_from_transcript(crs, witness_acc, transcript).proof;
    return result;
}

DoryChallenges DeriveChallenges(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryProof& proof) {
    if (!ValidateCRS(crs) || proof.rounds.size() != crs.d) {
        throw std::invalid_argument("invalid CRS or proof round count");
    }
    return replay_dory_transcript(
        crs, proof, dory_initial(crs, statement)).challenges;
}

std::vector<Fr> DeriveBatchGammas(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof) {
    if (!ValidateCRS(crs) || statements.empty()
        || proof.batchCrossTerms.size() != statements.size() - 1) {
        throw std::invalid_argument("invalid batch proof shape");
    }
    Digest transcript = dory_batch_initial(crs, statements);
    std::vector<Fr> out;
    out.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        transcript = dory_absorb_merge(
            transcript, j, proof.batchCrossTerms[j - 1]);
        out.push_back(
            ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j));
    }
    return out;
}

DoryBatchChallenges DeriveBatchChallenges(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof) {
    if (!ValidateCRS(crs) || statements.empty()
        || proof.batchCrossTerms.size() != statements.size() - 1
        || proof.doryProof.rounds.size() != crs.d) {
        throw std::invalid_argument("invalid batch proof shape");
    }
    Digest transcript = dory_batch_initial(crs, statements);
    DoryBatchChallenges out;
    out.gamma.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        transcript = dory_absorb_merge(
            transcript, j, proof.batchCrossTerms[j - 1]);
        out.gamma.push_back(
            ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j));
    }
    transcript = dory_enter_batch(transcript);
    out.dory = replay_dory_transcript(
        crs, proof.doryProof, transcript).challenges;
    return out;
}

} 

#include "rexp/dory.hpp"
#include "rexp/dory_verify.hpp"

#include <mcl/fp.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <chrono>
#include <utility>

namespace rexp {
namespace {

using Bytes = std::vector<std::uint8_t>;

void append(Bytes& out, const void* data, std::size_t size) {
    const auto* first = static_cast<const std::uint8_t*>(data);
    out.insert(out.end(), first, first + size);
}

void append(Bytes& out, std::string_view text) {
    append(out, text.data(), text.size());
}

void append_u64(Bytes& out, std::uint64_t value) {
    for (int shift = 56; shift >= 0; shift -= 8) {
        out.push_back(static_cast<std::uint8_t>(value >> shift));
    }
}

void frame(Bytes& out, const void* data, std::size_t size) {
    append_u64(out, static_cast<std::uint64_t>(size));
    append(out, data, size);
}

void frame(Bytes& out, std::string_view text) {
    frame(out, text.data(), text.size());
}

template<class T>
Bytes encode(const T& value) {
    Bytes out(2048);
    const std::size_t written = value.serialize(out.data(), out.size());
    if (written == 0) throw std::runtime_error("mcl serialization failed");
    out.resize(written);
    return out;
}

template<class T>
void frame_element(Bytes& out, const T& value) {
    const Bytes bytes = encode(value);
    frame(out, bytes.data(), bytes.size());
}

Digest hash(const Bytes& input) {
    Digest out{};
    const std::uint32_t written = mcl::fp::sha256(
        out.data(), static_cast<std::uint32_t>(out.size()),
        input.data(), static_cast<std::uint32_t>(input.size()));
    if (written != out.size()) throw std::runtime_error("SHA-256 failed");
    return out;
}

void frame_digest(Bytes& out, const Digest& digest) {
    frame(out, digest.data(), digest.size());
}

GT pair_product(
    const std::vector<G1>& left,
    std::size_t left_offset,
    const std::vector<G2>& right,
    std::size_t right_offset,
    std::size_t count) {
    if (left_offset + count > left.size() || right_offset + count > right.size()) {
        throw std::invalid_argument("pair-product slice is out of range");
    }
    std::vector<G1> a(left.begin() + left_offset, left.begin() + left_offset + count);
    std::vector<G2> b(right.begin() + right_offset, right.begin() + right_offset + count);
    GT miller;
    mcl::bn::millerLoopVec(miller, a.data(), b.data(), count, true);
    GT out;
    mcl::bn::finalExp(out, miller);
    return out;
}

Fr inverse(const Fr& value) {
    if (value.isZero()) throw std::invalid_argument("cannot invert zero challenge");
    Fr out;
    Fr::inv(out, value);
    return out;
}

Fr fr_mul(const Fr& a, const Fr& b) {
    Fr out;
    Fr::mul(out, a, b);
    return out;
}

G1 g1_mul(const G1& point, const Fr& scalar) {
    G1 out;
    G1::mul(out, point, scalar);
    return out;
}

G2 g2_mul(const G2& point, const Fr& scalar) {
    G2 out;
    G2::mul(out, point, scalar);
    return out;
}

G1 g1_add(const G1& a, const G1& b) {
    G1 out;
    G1::add(out, a, b);
    return out;
}

G2 g2_add(const G2& a, const G2& b) {
    G2 out;
    G2::add(out, a, b);
    return out;
}

GT gt_mul(const GT& a, const GT& b) {
    GT out;
    GT::mul(out, a, b);
    return out;
}

GT gt_pow(const GT& value, const Fr& scalar) {
    GT out;
    GT::pow(out, value, scalar);
    return out;
}

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

Digest ordinary_initial(
    const DoryCRS& crs,
    const DoryStatement& statement) {
    Bytes input;
    frame(input, "DORY-NI-NONZK-V1");
    frame_digest(input, crs.digest);
    append_u64(input, crs.d);
    append_u64(input, crs.n);
    frame_element(input, statement.D0);
    frame_element(input, statement.D1);
    frame_element(input, statement.D2);
    return hash(input);
}

Digest batch_initial(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements) {
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
    return hash(input);
}

Digest absorb_beta(
    const Digest& current,
    std::size_t k,
    std::size_t m,
    const DoryRound& round) {
    Bytes input;
    frame(input, "DORY-FOLD-BETA-MESSAGE-V1");
    frame_digest(input, current);
    append_u64(input, k);
    append_u64(input, m);
    frame_element(input, round.D1L);
    frame_element(input, round.D1R);
    frame_element(input, round.D2L);
    frame_element(input, round.D2R);
    return hash(input);
}

Digest absorb_alpha(
    const Digest& beta_digest,
    std::size_t k,
    std::size_t m,
    const DoryRound& round) {
    Bytes input;
    frame(input, "DORY-FOLD-ALPHA-MESSAGE-V1");
    frame_digest(input, beta_digest);
    append_u64(input, k);
    append_u64(input, m);
    frame_element(input, round.W1);
    frame_element(input, round.W2);
    return hash(input);
}

Digest absorb_final(
    const Digest& current,
    const G1& phi,
    const G2& theta) {
    Bytes input;
    frame(input, "DORY-FINAL-CHECK-V1");
    frame_digest(input, current);
    frame_element(input, phi);
    frame_element(input, theta);
    return hash(input);
}

Digest absorb_merge(const Digest& current, std::size_t j, const GT& cross) {
    Bytes input;
    frame(input, "DORY-BATCH-MERGE-V1");
    frame_digest(input, current);
    append_u64(input, j);
    frame_element(input, cross);
    return hash(input);
}

Digest enter_batch_dory(const Digest& current) {
    Bytes input;
    frame(input, "DORY-BATCH-AGGREGATE-SUBPROTOCOL-V1");
    frame_digest(input, current);
    return hash(input);
}

struct FoldResult {
    DoryProof proof;
    DoryChallenges challenges;
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
        round.D1L = pair_product(phi, 0, crs.Lambda, 0, h);
        round.D1R = pair_product(phi, h, crs.Lambda, 0, h);
        round.D2L = pair_product(crs.Gamma, 0, theta, 0, h);
        round.D2R = pair_product(crs.Gamma, 0, theta, h, h);

        const Digest beta_digest = absorb_beta(transcript, k, m, round);
        const Fr beta = ChallengeNonzeroFr(beta_digest, "DORY-BETA-V1", k);
        const Fr beta_inv = inverse(beta);

        std::vector<G1> phi_circle(m);
        std::vector<G2> theta_circle(m);
        for (std::size_t i = 0; i < m; ++i) {
            phi_circle[i] = g1_add(phi[i], g1_mul(crs.Gamma[i], beta));
            theta_circle[i] = g2_add(theta[i], g2_mul(crs.Lambda[i], beta_inv));
        }
        round.W1 = pair_product(phi_circle, 0, theta_circle, h, h);
        round.W2 = pair_product(phi_circle, h, theta_circle, 0, h);

        const Digest alpha_digest = absorb_alpha(beta_digest, k, m, round);
        const Fr alpha = ChallengeNonzeroFr(alpha_digest, "DORY-ALPHA-V1", k);
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
    const Digest final_digest =
        absorb_final(transcript, result.proof.PhiFinal, result.proof.ThetaFinal);
    result.challenges.epsilon =
        ChallengeNonzeroFr(final_digest, "DORY-EPSILON-V1", crs.d);
    return result;
}

bool digest_less_than_limit(const Digest& value) {

    static constexpr Digest limit = {
        0xf1,0xf5,0x88,0x3e,0x65,0xf8,0x20,0xd0,
        0x99,0x91,0x5c,0x90,0x87,0x86,0xb9,0xd1,
        0xc9,0x03,0x89,0x6a,0x60,0x9f,0x32,0xd6,
        0x53,0x69,0xcb,0xe3,0xb0,0x00,0x00,0x05};
    return std::lexicographical_compare(
        value.begin(), value.end(), limit.begin(), limit.end());
}

using Clock = std::chrono::steady_clock;

double elapsed_ms(Clock::time_point begin, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct ExpressionEntry {
    GT base;
    Fr coefficient;
};

class TargetExpression {
public:
    explicit TargetExpression(VerifyMetrics* metrics) : metrics_(metrics) {}

    void atom(const GT& base, const Fr& coefficient) {
        if (metrics_) ++metrics_->symbolic_atom_insertions;
        add(base, coefficient);
    }

    void merge(const TargetExpression& other, const Fr& scale) {
        for (const auto& item : other.terms_) {
            add(item.base, fr_mul(item.coefficient, scale));
        }
    }

    void merge(const TargetExpression& other) {
        Fr one;
        one = 1;
        merge(other, one);
    }

    const std::vector<ExpressionEntry>& terms() const { return terms_; }

private:
    void add(const GT& base, const Fr& coefficient) {
        if (coefficient.isZero()) {
            if (metrics_) ++metrics_->zero_coefficients_removed;
            return;
        }
        if (base.isOne()) {
            if (metrics_) ++metrics_->identity_bases_removed;
            return;
        }
        auto found = std::find_if(
            terms_.begin(), terms_.end(),
            [&](const ExpressionEntry& entry) { return entry.base == base; });
        if (found == terms_.end()) {
            terms_.push_back(ExpressionEntry{base, coefficient});
            return;
        }
        if (metrics_) ++metrics_->coalesced_duplicate_bases;
        Fr combined;
        Fr::add(combined, found->coefficient, coefficient);
        if (combined.isZero()) {
            terms_.erase(found);
            if (metrics_) ++metrics_->zero_coefficients_removed;
        } else {
            found->coefficient = combined;
        }
    }

    VerifyMetrics* metrics_;
    std::vector<ExpressionEntry> terms_;
};

std::vector<Fr> batch_invert(const std::vector<Fr>& values) {
    if (values.empty()) return {};
    std::vector<Fr> prefix(values.size());
    Fr accumulator;
    accumulator = 1;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (values[i].isZero()) throw std::invalid_argument("zero challenge");
        Fr::mul(accumulator, accumulator, values[i]);
        prefix[i] = accumulator;
    }
    Fr inverse_accumulator;
    Fr::inv(inverse_accumulator, accumulator);
    std::vector<Fr> inverses(values.size());
    for (std::size_t i = values.size(); i-- > 0;) {
        Fr before;
        if (i == 0) before = 1;
        else before = prefix[i - 1];
        Fr::mul(inverses[i], inverse_accumulator, before);
        Fr::mul(inverse_accumulator, inverse_accumulator, values[i]);
    }
    return inverses;
}

struct Replay {
    DoryChallenges challenges;
    Digest end;
};

Replay replay_from(
    const DoryCRS& crs,
    const DoryProof& proof,
    Digest transcript) {
    Replay out;
    out.challenges.beta.reserve(crs.d);
    out.challenges.alpha.reserve(crs.d);
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        const Digest beta_digest = absorb_beta(transcript, k, m, proof.rounds[k]);
        out.challenges.beta.push_back(
            ChallengeNonzeroFr(beta_digest, "DORY-BETA-V1", k));
        transcript = absorb_alpha(beta_digest, k, m, proof.rounds[k]);
        out.challenges.alpha.push_back(
            ChallengeNonzeroFr(transcript, "DORY-ALPHA-V1", k));
    }
    const Digest final_digest =
        absorb_final(transcript, proof.PhiFinal, proof.ThetaFinal);
    out.challenges.epsilon =
        ChallengeNonzeroFr(final_digest, "DORY-EPSILON-V1", crs.d);
    out.end = final_digest;
    return out;
}

bool valid_proof_shape(const DoryCRS& crs, const DoryProof& proof) {
    if (proof.rounds.size() != crs.d
        || !proof.PhiFinal.isValid() || !proof.PhiFinal.isValidOrder()
        || !proof.ThetaFinal.isValid() || !proof.ThetaFinal.isValidOrder()) {
        return false;
    }
    return true;
}

GT evaluate_expression(
    const TargetExpression& expression,
    VerifyMetrics* metrics,
    VerifyAuditTrace* trace) {
    std::vector<GT> bases;
    std::vector<Fr> coefficients;
    bases.reserve(expression.terms().size());
    coefficients.reserve(expression.terms().size());
    Fr one;
    one = 1;
    for (const auto& item : expression.terms()) {
        bases.push_back(item.base);
        coefficients.push_back(item.coefficient);
        if (metrics && item.coefficient == one) {
            ++metrics->coefficient_one_bases;
        }
    }
    if (metrics) {
        metrics->actual_gt_bases = bases.size();
        metrics->nonzero_gt_coefficients = coefficients.size();
    }
    if (trace) {
        trace->gt_bases = bases;
        trace->gt_exponents = coefficients;
    }
    GT result;
    if (bases.empty()) result.setOne();
    else {
        GT::powVec(result, bases.data(), coefficients.data(), bases.size());
        if (metrics) ++metrics->gt_multiexp_calls;
    }
    return result;
}

TargetExpression deferred_dory_expression(
    const DoryPrecomputation& precomp,
    const DoryProof& proof,
    const DoryChallenges& challenges,
    const std::vector<Fr>& inverses,
    const Fr& epsilon_inv,
    TargetExpression e0,
    TargetExpression e1,
    TargetExpression e2,
    VerifyMetrics* metrics) {
    Fr one;
    one = 1;
    for (std::size_t k = 0; k < proof.rounds.size(); ++k) {
        const DoryRound& round = proof.rounds[k];
        TargetExpression next0(metrics);
        next0.merge(e0);
        next0.atom(precomp.X[k], one);
        next0.merge(e1, inverses[k]);
        next0.merge(e2, challenges.beta[k]);
        next0.atom(round.W1, challenges.alpha[k]);
        next0.atom(round.W2, inverses[challenges.beta.size() + k]);

        TargetExpression next1(metrics);
        next1.atom(round.D1L, challenges.alpha[k]);
        next1.atom(round.D1R, one);
        next1.atom(
            precomp.X[k + 1], fr_mul(challenges.alpha[k], challenges.beta[k]));
        next1.atom(precomp.Delta1R[k], challenges.beta[k]);

        TargetExpression next2(metrics);
        next2.atom(round.D2L, inverses[challenges.beta.size() + k]);
        next2.atom(round.D2R, one);
        next2.atom(
            precomp.X[k + 1], fr_mul(
                inverses[challenges.beta.size() + k], inverses[k]));
        next2.atom(precomp.Delta2R[k], inverses[k]);
        e0 = std::move(next0);
        e1 = std::move(next1);
        e2 = std::move(next2);
    }
    TargetExpression lhs(metrics);
    lhs.merge(e0);
    lhs.merge(e1, epsilon_inv);
    lhs.merge(e2, challenges.epsilon);
    lhs.atom(precomp.X.back(), one);
    return lhs;
}

}

Fr ChallengeNonzeroFr(
    const Digest& transcript_digest,
    std::string_view label,
    std::size_t index) {
    initialize();
    for (std::uint64_t counter = 0;; ++counter) {
        Bytes input;
        frame(input, "FS-NONZERO-FR-V1");
        frame(input, label);
        frame_digest(input, transcript_digest);
        append_u64(input, index);
        append_u64(input, counter);
        const Digest candidate = hash(input);
        if (!digest_less_than_limit(candidate)) continue;
        Fr value;
        value.setBigEndianMod(candidate.data(), candidate.size());
        if (!value.isZero()) return value;
        if (counter == std::numeric_limits<std::uint64_t>::max()) {
            throw std::runtime_error("Fiat-Shamir rejection counter exhausted");
        }
    }
}

DoryProof Prove(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryWitness& witness) {
    validate_witness(crs, witness);
    return prove_from_transcript(
        crs, witness, ordinary_initial(crs, statement)).proof;
}

DoryProof ProveEmbedded(
    const DoryCRS& crs,
    const DoryWitness& witness,
    const Digest& transcript_in,
    Digest* transcript_end) {
    validate_witness(crs, witness);
    FoldResult result = prove_from_transcript(crs, witness, transcript_in);
    if (transcript_end) {
        *transcript_end = absorb_final(
            [&] {
                Digest t = transcript_in;
                for (std::size_t k = 0; k < crs.d; ++k) {
                    const std::size_t m = crs.n >> k;
                    const Digest b = absorb_beta(t, k, m, result.proof.rounds[k]);
                    t = absorb_alpha(b, k, m, result.proof.rounds[k]);
                }
                return t;
            }(),
            result.proof.PhiFinal,
            result.proof.ThetaFinal);
    }
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

    Digest transcript = batch_initial(crs, statements);
    DoryStatement statement_acc = statements.front();
    DoryWitness witness_acc = witnesses.front();
    DoryBatchProof result;
    result.batchCrossTerms.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        const GT left = pair_product(
            witness_acc.Phi, 0, witnesses[j].Theta, 0, crs.n);
        const GT right = pair_product(
            witnesses[j].Phi, 0, witness_acc.Theta, 0, crs.n);
        const GT cross = gt_mul(left, right);
        result.batchCrossTerms.push_back(cross);
        transcript = absorb_merge(transcript, j, cross);
        const Fr gamma =
            ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j);
        const Fr gamma_sq = fr_mul(gamma, gamma);
        statement_acc.D0 = gt_mul(
            gt_mul(gt_pow(statement_acc.D0, gamma_sq), gt_pow(cross, gamma)),
            statements[j].D0);
        statement_acc.D1 = gt_mul(
            gt_pow(statement_acc.D1, gamma), statements[j].D1);
        statement_acc.D2 = gt_mul(
            gt_pow(statement_acc.D2, gamma), statements[j].D2);
        for (std::size_t i = 0; i < crs.n; ++i) {
            witness_acc.Phi[i] = g1_add(
                g1_mul(witness_acc.Phi[i], gamma), witnesses[j].Phi[i]);
            witness_acc.Theta[i] = g2_add(
                g2_mul(witness_acc.Theta[i], gamma), witnesses[j].Theta[i]);
        }
    }
    transcript = enter_batch_dory(transcript);
    result.doryProof =
        prove_from_transcript(crs, witness_acc, transcript).proof;
    return result;
}

DoryChallenges DeriveChallenges(
    const DoryCRS& crs,
    const DoryStatement& statement,
    const DoryProof& proof) {
    if (!ValidateCRS(crs) || proof.rounds.size() != crs.d) {
        throw std::invalid_argument("invalid CRS or proof round count");
    }
    Digest transcript = ordinary_initial(crs, statement);
    DoryChallenges out;
    out.beta.reserve(crs.d);
    out.alpha.reserve(crs.d);
    for (std::size_t k = 0; k < crs.d; ++k) {
        const std::size_t m = crs.n >> k;
        const Digest beta_digest = absorb_beta(transcript, k, m, proof.rounds[k]);
        out.beta.push_back(
            ChallengeNonzeroFr(beta_digest, "DORY-BETA-V1", k));
        transcript = absorb_alpha(beta_digest, k, m, proof.rounds[k]);
        out.alpha.push_back(
            ChallengeNonzeroFr(transcript, "DORY-ALPHA-V1", k));
    }
    const Digest final_digest =
        absorb_final(transcript, proof.PhiFinal, proof.ThetaFinal);
    out.epsilon =
        ChallengeNonzeroFr(final_digest, "DORY-EPSILON-V1", crs.d);
    return out;
}

std::vector<Fr> DeriveBatchGammas(
    const DoryCRS& crs,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof) {
    if (!ValidateCRS(crs) || statements.empty()
        || proof.batchCrossTerms.size() != statements.size() - 1) {
        throw std::invalid_argument("invalid batch proof shape");
    }
    Digest transcript = batch_initial(crs, statements);
    std::vector<Fr> out;
    out.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        transcript = absorb_merge(transcript, j, proof.batchCrossTerms[j - 1]);
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
    Digest transcript = batch_initial(crs, statements);
    DoryBatchChallenges out;
    out.gamma.reserve(statements.size() - 1);
    for (std::size_t j = 1; j < statements.size(); ++j) {
        transcript = absorb_merge(transcript, j, proof.batchCrossTerms[j - 1]);
        out.gamma.push_back(
            ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j));
    }
    transcript = enter_batch_dory(transcript);
    out.dory = replay_from(crs, proof.doryProof, transcript).challenges;
    return out;
}

std::vector<std::uint8_t> SerializeProof(
    const DoryProof& proof,
    std::size_t d) {
    if (proof.rounds.size() != d) {
        throw std::invalid_argument("proof round count differs from d");
    }
    Bytes out;
    frame(out, "DORY-PROOF-WIRE-BN254-V1");
    frame(out, "BN254");
    append_u64(out, d);
    append_u64(out, proof.rounds.size());
    for (const DoryRound& round : proof.rounds) {
        frame_element(out, round.D1L);
        frame_element(out, round.D1R);
        frame_element(out, round.D2L);
        frame_element(out, round.D2R);
        frame_element(out, round.W1);
        frame_element(out, round.W2);
    }
    frame_element(out, proof.PhiFinal);
    frame_element(out, proof.ThetaFinal);
    return out;
}

std::vector<std::uint8_t> SerializeBatchProof(
    const DoryBatchProof& proof,
    std::size_t d,
    std::size_t batch_size) {
    if (batch_size == 0 || proof.batchCrossTerms.size() != batch_size - 1) {
        throw std::invalid_argument("batch proof cross-term count mismatch");
    }
    Bytes out;
    frame(out, "DORY-BATCH-PROOF-WIRE-BN254-V1");
    frame(out, "BN254");
    append_u64(out, d);
    append_u64(out, batch_size);
    append_u64(out, proof.batchCrossTerms.size());
    for (const GT& cross : proof.batchCrossTerms) frame_element(out, cross);
    const Bytes ordinary = SerializeProof(proof.doryProof, d);
    frame(out, ordinary.data(), ordinary.size());
    return out;
}

ProofSizes MeasureProofSizes(const DoryProof& proof, std::size_t d) {
    if (proof.rounds.size() != d) {
        throw std::invalid_argument("proof round count differs from d");
    }
    std::size_t payload = 0;
    for (const DoryRound& round : proof.rounds) {
        payload += encode(round.D1L).size();
        payload += encode(round.D1R).size();
        payload += encode(round.D2L).size();
        payload += encode(round.D2R).size();
        payload += encode(round.W1).size();
        payload += encode(round.W2).size();
    }
    payload += encode(proof.PhiFinal).size();
    payload += encode(proof.ThetaFinal).size();
    return {payload, SerializeProof(proof, d).size()};
}

ProofSizes MeasureProofSizes(
    const DoryBatchProof& proof,
    std::size_t d,
    std::size_t batch_size) {
    if (batch_size == 0 || proof.batchCrossTerms.size() != batch_size - 1) {
        throw std::invalid_argument("batch proof cross-term count mismatch");
    }
    std::size_t payload =
        MeasureProofSizes(proof.doryProof, d).mathematical_payload_bytes;
    for (const GT& cross : proof.batchCrossTerms) {
        payload += encode(cross).size();
    }
    return {payload, SerializeBatchProof(proof, d, batch_size).size()};
}

bool verify_impl(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics* metrics,
    bool validate_parameters,
    VerifyAuditTrace* trace = nullptr,
    const Digest* embedded_transcript = nullptr,
    Digest* transcript_end = nullptr) {
    const Clock::time_point total_start =
        metrics ? Clock::now() : Clock::time_point{};
    if (metrics) *metrics = VerifyMetrics{};
    try {
        if ((validate_parameters && !ValidatePrecomputation(crs, precomp))
            || !valid_proof_shape(crs, proof)) return false;

        const Clock::time_point transcript_start =
            metrics ? Clock::now() : Clock::time_point{};
        const Replay replay = replay_from(
            crs, proof,
            embedded_transcript ? *embedded_transcript
                                : ordinary_initial(crs, statement));
        if (transcript_end) *transcript_end = replay.end;
        const Clock::time_point transcript_stop =
            metrics ? Clock::now() : Clock::time_point{};

        std::vector<Fr> all(2 * crs.d + 1);
        std::copy(replay.challenges.beta.begin(), replay.challenges.beta.end(),
                  all.begin());
        std::copy(replay.challenges.alpha.begin(), replay.challenges.alpha.end(),
                  all.begin() + crs.d);
        all.back() = replay.challenges.epsilon;
        const Clock::time_point inversion_start =
            metrics ? Clock::now() : Clock::time_point{};
        const std::vector<Fr> inverses = batch_invert(all);
        const Clock::time_point inversion_stop =
            metrics ? Clock::now() : Clock::time_point{};
        const Fr epsilon_inv = inverses.back();

        const Clock::time_point symbolic_start =
            metrics ? Clock::now() : Clock::time_point{};
        Fr one;
        one = 1;
        TargetExpression e0(metrics), e1(metrics), e2(metrics);
        e0.atom(statement.D0, one);
        e1.atom(statement.D1, one);
        e2.atom(statement.D2, one);
        const TargetExpression lhs_expression = deferred_dory_expression(
            precomp, proof, replay.challenges, inverses,
            epsilon_inv, std::move(e0), std::move(e1), std::move(e2), metrics);
        const Clock::time_point symbolic_stop =
            metrics ? Clock::now() : Clock::time_point{};

        const G1 left = g1_add(
            proof.PhiFinal, g1_mul(crs.Gamma.front(), replay.challenges.epsilon));
        const G2 right = g2_add(
            proof.ThetaFinal, g2_mul(crs.Lambda.front(), epsilon_inv));
        if (trace) {
            trace->terminal_left_g1 = left;
            trace->terminal_right_g2 = right;
        }
        const Clock::time_point pairing_start =
            metrics ? Clock::now() : Clock::time_point{};
        GT rhs;
        mcl::bn::pairing(rhs, left, right);
        if (trace) trace->terminal_pairing_rhs = rhs;
        const Clock::time_point pairing_stop =
            metrics ? Clock::now() : Clock::time_point{};
        if (metrics) metrics->terminal_pairings = 1;

        const Clock::time_point multiexp_start =
            metrics ? Clock::now() : Clock::time_point{};
        const GT lhs = evaluate_expression(lhs_expression, metrics, trace);
        if (trace) trace->evaluated_lhs = lhs;
        const Clock::time_point multiexp_stop =
            metrics ? Clock::now() : Clock::time_point{};
        if (metrics) {
            metrics->transcript_ms =
                elapsed_ms(transcript_start, transcript_stop);
            metrics->batch_inversion_ms =
                elapsed_ms(inversion_start, inversion_stop);
            metrics->symbolic_replay_ms =
                elapsed_ms(symbolic_start, symbolic_stop);
            metrics->terminal_pairing_ms =
                elapsed_ms(pairing_start, pairing_stop);
            metrics->gt_multiexp_ms =
                elapsed_ms(multiexp_start, multiexp_stop);
            metrics->total_ms = elapsed_ms(total_start, Clock::now());
        }
        return lhs == rhs;
    } catch (...) {
        if (metrics) metrics->total_ms = elapsed_ms(total_start, Clock::now());
        return false;
    }
}

bool verify_batch_impl(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics* metrics,
    bool validate_parameters,
    VerifyAuditTrace* trace = nullptr) {
    const Clock::time_point total_start = Clock::now();
    if (metrics) *metrics = VerifyMetrics{};
    try {
        if ((validate_parameters && !ValidatePrecomputation(crs, precomp))
            || statements.empty()
            || proof.batchCrossTerms.size() != statements.size() - 1
            || !valid_proof_shape(crs, proof.doryProof)) return false;

        const Clock::time_point transcript_start = Clock::now();
        Digest transcript = batch_initial(crs, statements);
        std::vector<Fr> gammas;
        gammas.reserve(statements.size() - 1);
        for (std::size_t j = 1; j < statements.size(); ++j) {
            transcript = absorb_merge(transcript, j, proof.batchCrossTerms[j - 1]);
            gammas.push_back(
                ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j));
        }
        transcript = enter_batch_dory(transcript);
        const Replay replay = replay_from(crs, proof.doryProof, transcript);
        const Clock::time_point transcript_stop = Clock::now();

        std::vector<Fr> all(2 * crs.d + 1);
        std::copy(replay.challenges.beta.begin(), replay.challenges.beta.end(),
                  all.begin());
        std::copy(replay.challenges.alpha.begin(), replay.challenges.alpha.end(),
                  all.begin() + crs.d);
        all.back() = replay.challenges.epsilon;
        const Clock::time_point inversion_start = Clock::now();
        const std::vector<Fr> inverses = batch_invert(all);
        const Clock::time_point inversion_stop = Clock::now();
        const Fr epsilon_inv = inverses.back();

        const Clock::time_point symbolic_start = Clock::now();
        Fr one;
        one = 1;
        TargetExpression a0(metrics), a1(metrics), a2(metrics);
        a0.atom(statements.front().D0, one);
        a1.atom(statements.front().D1, one);
        a2.atom(statements.front().D2, one);
        for (std::size_t j = 1; j < statements.size(); ++j) {
            const Fr gamma_sq = fr_mul(gammas[j - 1], gammas[j - 1]);
            TargetExpression next0(metrics), next1(metrics), next2(metrics);
            next0.merge(a0, gamma_sq);
            next0.atom(proof.batchCrossTerms[j - 1], gammas[j - 1]);
            next0.atom(statements[j].D0, one);
            next1.merge(a1, gammas[j - 1]);
            next1.atom(statements[j].D1, one);
            next2.merge(a2, gammas[j - 1]);
            next2.atom(statements[j].D2, one);
            a0 = std::move(next0);
            a1 = std::move(next1);
            a2 = std::move(next2);
        }
        const TargetExpression lhs_expression = deferred_dory_expression(
            precomp, proof.doryProof, replay.challenges, inverses,
            epsilon_inv, std::move(a0), std::move(a1), std::move(a2), metrics);
        const Clock::time_point symbolic_stop = Clock::now();

        const G1 left = g1_add(
            proof.doryProof.PhiFinal,
            g1_mul(crs.Gamma.front(), replay.challenges.epsilon));
        const G2 right = g2_add(
            proof.doryProof.ThetaFinal,
            g2_mul(crs.Lambda.front(), epsilon_inv));
        if (trace) {
            trace->terminal_left_g1 = left;
            trace->terminal_right_g2 = right;
        }
        const Clock::time_point pairing_start = Clock::now();
        GT rhs;
        mcl::bn::pairing(rhs, left, right);
        if (trace) trace->terminal_pairing_rhs = rhs;
        const Clock::time_point pairing_stop = Clock::now();
        if (metrics) metrics->terminal_pairings = 1;
        const Clock::time_point multiexp_start = Clock::now();
        const GT lhs = evaluate_expression(lhs_expression, metrics, trace);
        if (trace) trace->evaluated_lhs = lhs;
        const Clock::time_point multiexp_stop = Clock::now();
        if (metrics) {
            metrics->transcript_ms =
                elapsed_ms(transcript_start, transcript_stop);
            metrics->batch_inversion_ms =
                elapsed_ms(inversion_start, inversion_stop);
            metrics->symbolic_replay_ms =
                elapsed_ms(symbolic_start, symbolic_stop);
            metrics->terminal_pairing_ms =
                elapsed_ms(pairing_start, pairing_stop);
            metrics->gt_multiexp_ms =
                elapsed_ms(multiexp_start, multiexp_stop);
            metrics->total_ms = elapsed_ms(total_start, Clock::now());
        }
        return lhs == rhs;
    } catch (...) {
        if (metrics) metrics->total_ms = elapsed_ms(total_start, Clock::now());
        return false;
    }
}

std::optional<PreparedVerifier> PrepareVerifier(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp) {


    if (!ValidateCRS(crs) || !ValidatePrecomputation(crs, precomp)) {
        return std::nullopt;
    }
    return PreparedVerifier(crs, precomp);
}

bool Verify(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics* metrics) {
    return verify_impl(crs, precomp, statement, proof, metrics, true);
}

bool VerifyBatch(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics* metrics) {
    return verify_batch_impl(
        crs, precomp, statements, proof, metrics, true);
}

bool VerifyPrepared(
    const PreparedVerifier& prepared,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics* metrics) {
    return verify_impl(
        prepared.crs(), prepared.precomp(), statement, proof, metrics, false);
}

bool VerifyEmbeddedDeferred(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    const Digest& transcript_in,
    Digest* transcript_end,
    VerifyMetrics* metrics) {
    return verify_impl(
        crs, precomp, statement, proof, metrics, false, nullptr,
        &transcript_in, transcript_end);
}

bool VerifyEmbeddedReference(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp,
    const DoryStatement& statement,
    const DoryProof& proof,
    const Digest& transcript_in,
    Digest* transcript_end) {
    try {
        if (!ValidateCRS(crs) || !ValidatePrecomputation(crs, precomp)
            || !valid_proof_shape(crs, proof)) return false;
        const Replay replay = replay_from(crs, proof, transcript_in);
        if (transcript_end) *transcript_end = replay.end;
        GT d0 = statement.D0;
        GT d1 = statement.D1;
        GT d2 = statement.D2;
        for (std::size_t k = 0; k < crs.d; ++k) {
            const Fr beta_inv = inverse(replay.challenges.beta[k]);
            const Fr alpha_inv = inverse(replay.challenges.alpha[k]);
            const DoryRound& r = proof.rounds[k];
            GT next0 = gt_mul(d0, precomp.X[k]);
            next0 = gt_mul(next0, gt_pow(d1, beta_inv));
            next0 = gt_mul(next0, gt_pow(d2, replay.challenges.beta[k]));
            next0 = gt_mul(next0, gt_pow(r.W1, replay.challenges.alpha[k]));
            next0 = gt_mul(next0, gt_pow(r.W2, alpha_inv));
            GT next1 = gt_mul(
                gt_pow(r.D1L, replay.challenges.alpha[k]), r.D1R);
            next1 = gt_mul(next1, gt_pow(
                precomp.X[k + 1],
                fr_mul(replay.challenges.alpha[k],
                       replay.challenges.beta[k])));
            next1 = gt_mul(
                next1, gt_pow(precomp.Delta1R[k],
                              replay.challenges.beta[k]));
            GT next2 = gt_mul(gt_pow(r.D2L, alpha_inv), r.D2R);
            next2 = gt_mul(next2, gt_pow(
                precomp.X[k + 1], fr_mul(alpha_inv, beta_inv)));
            next2 = gt_mul(
                next2, gt_pow(precomp.Delta2R[k], beta_inv));
            d0 = std::move(next0);
            d1 = std::move(next1);
            d2 = std::move(next2);
        }
        const Fr epsilon_inv = inverse(replay.challenges.epsilon);
        GT lhs = gt_mul(d0, gt_pow(d1, epsilon_inv));
        lhs = gt_mul(lhs, gt_pow(d2, replay.challenges.epsilon));
        lhs = gt_mul(lhs, precomp.X.back());
        const G1 left = g1_add(
            proof.PhiFinal,
            g1_mul(crs.Gamma.front(), replay.challenges.epsilon));
        const G2 right = g2_add(
            proof.ThetaFinal,
            g2_mul(crs.Lambda.front(), epsilon_inv));
        GT rhs;
        mcl::bn::pairing(rhs, left, right);
        return lhs == rhs;
    } catch (...) {
        return false;
    }
}

bool VerifyBatchPrepared(
    const PreparedVerifier& prepared,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics* metrics) {
    return verify_batch_impl(
        prepared.crs(), prepared.precomp(), statements, proof, metrics, false);
}

bool VerifyPreparedAudit(
    const PreparedVerifier& prepared,
    const DoryStatement& statement,
    const DoryProof& proof,
    VerifyMetrics& metrics,
    VerifyAuditTrace& trace) {
    trace = VerifyAuditTrace{};
    return verify_impl(
        prepared.crs(), prepared.precomp(), statement, proof, &metrics, false,
        &trace);
}

bool VerifyBatchPreparedAudit(
    const PreparedVerifier& prepared,
    const std::vector<DoryStatement>& statements,
    const DoryBatchProof& proof,
    VerifyMetrics& metrics,
    VerifyAuditTrace& trace) {
    trace = VerifyAuditTrace{};
    return verify_batch_impl(
        prepared.crs(), prepared.precomp(), statements, proof, &metrics, false,
        &trace);
}

}

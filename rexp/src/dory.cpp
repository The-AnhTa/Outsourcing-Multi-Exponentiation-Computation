#include "rexp/dory.hpp"
#include "rexp/dory_verify.hpp"
#include "internal/crypto.hpp"
#include "internal/dory_transcript.hpp"

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

using internal::Bytes;
using internal::append;
using internal::append_u64;
using internal::encode;
using internal::frame;
using internal::frame_digest;
using internal::frame_element;
using internal::fr_mul;
using internal::g1_add;
using internal::g1_mul;
using internal::g2_add;
using internal::g2_mul;
using internal::gt_mul;
using internal::gt_pow;
using internal::inverse;
using internal::pairing_product;
using internal::sha256;
using internal::DoryTranscriptReplay;
using internal::dory_absorb_alpha;
using internal::dory_absorb_beta;
using internal::dory_absorb_final;
using internal::dory_absorb_merge;
using internal::dory_batch_initial;
using internal::dory_enter_batch;
using internal::dory_initial;
using internal::replay_dory_transcript;

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
        const DoryTranscriptReplay replay = replay_dory_transcript(
            crs, proof,
            embedded_transcript ? *embedded_transcript
                                : dory_initial(crs, statement));
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
        Digest transcript = dory_batch_initial(crs, statements);
        std::vector<Fr> gammas;
        gammas.reserve(statements.size() - 1);
        for (std::size_t j = 1; j < statements.size(); ++j) {
            transcript = dory_absorb_merge(
                transcript, j, proof.batchCrossTerms[j - 1]);
            gammas.push_back(
                ChallengeNonzeroFr(transcript, "DORY-BATCH-GAMMA-V1", j));
        }
        transcript = dory_enter_batch(transcript);
        const DoryTranscriptReplay replay =
            replay_dory_transcript(crs, proof.doryProof, transcript);
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
        const DoryTranscriptReplay replay =
            replay_dory_transcript(crs, proof, transcript_in);
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

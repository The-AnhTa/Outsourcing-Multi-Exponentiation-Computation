#include "rexpbf/verify.hpp"

#include "rexpbf/serialization.hpp"

#include <chrono>
#include <limits>

namespace rexpbf {
namespace {
using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}
}

std::optional<ValidatedVerificationInputs> validate_verification_inputs(
    const CRS& crs,
    const Precomputation& precomputation,
    const Statement& statement,
    const Proof& proof,
    ValidationBreakdown* output) {
    ValidationBreakdown local;
    auto& breakdown = output ? *output : local;
    breakdown = {};
    const auto total_start = Clock::now();
    auto phase = total_start;

    const auto reject = [&]() -> std::optional<ValidatedVerificationInputs> {
        breakdown.total_ms = milliseconds(total_start, Clock::now());
        return std::nullopt;
    };

    try {
        const bool crs_shape = crs.d >= 1
            && crs.d < std::numeric_limits<std::size_t>::digits
            && crs.n == (std::size_t{1} << crs.d)
            && crs.gamma.size() == crs.n && crs.lambda.size() == crs.n;
        breakdown.crs_shape_ms = milliseconds(phase, Clock::now());
        if (!crs_shape) return reject();

        phase = Clock::now();
        const bool crs_digest = crs.digest == compute_crs_digest(crs);
        breakdown.crs_digest_ms = milliseconds(phase, Clock::now());
        if (!crs_digest) return reject();

        phase = Clock::now();
        for (const auto& point : crs.gamma) {
            ++breakdown.g1_elements_checked;
            ++breakdown.crs_g1_checked;
            if (!point.isValid() || !point.isValidOrder() || point.isZero())
                return reject();
        }
        for (const auto& point : crs.lambda) {
            ++breakdown.g2_elements_checked;
            ++breakdown.crs_g2_checked;
            if (!point.isValid() || !point.isValidOrder() || point.isZero())
                return reject();
        }
        breakdown.crs_group_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        const bool precomputation_shape = precomputation.x.size() == crs.d + 1
            && precomputation.delta1_right.size() == crs.d
            && precomputation.delta2_right.size() == crs.d
            && precomputation.pairing_terms == 4 * std::uint64_t(crs.n) - 3;
        breakdown.precomputation_shape_ms = milliseconds(phase, Clock::now());
        if (!precomputation_shape) return reject();

        phase = Clock::now();
        const bool precomputation_binding = precomputation.crs_digest == crs.digest;
        breakdown.precomputation_digest_binding_ms = milliseconds(phase, Clock::now());
        if (!precomputation_binding) return reject();

        const auto check_gt = [&](const GT& value) {
            const auto start = Clock::now();
            const bool valid = mcl::bn::isValidGT(value);
            breakdown.gt_subgroup_validation_ms += milliseconds(start, Clock::now());
            ++breakdown.gt_elements_checked;
            ++breakdown.gt_subgroup_checks;
            return valid;
        };
        phase = Clock::now();
        for (const auto& value : precomputation.x) {
            ++breakdown.precomputation_gt_checked;
            if (!check_gt(value)) return reject();
        }
        for (const auto& value : precomputation.delta1_right) {
            ++breakdown.precomputation_gt_checked;
            if (!check_gt(value)) return reject();
        }
        for (const auto& value : precomputation.delta2_right) {
            ++breakdown.precomputation_gt_checked;
            if (!check_gt(value)) return reject();
        }
        breakdown.precomputation_gt_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        const bool statement_shape = statement.crs_digest == crs.digest
            && statement.h.size() == crs.n
            && statement.pairing_terms == 3 * std::uint64_t(crs.n);
        breakdown.statement_shape_ms = milliseconds(phase, Clock::now());
        if (!statement_shape) return reject();

        phase = Clock::now();
        const bool statement_digest = statement.digest == compute_statement_digest(crs, statement);
        breakdown.statement_digest_ms = milliseconds(phase, Clock::now());
        if (!statement_digest) return reject();

        phase = Clock::now();
        for (const auto& point : statement.h) {
            ++breakdown.g1_elements_checked;
            ++breakdown.statement_g1_checked;
            if (!point.isValid() || !point.isValidOrder() || point.isZero())
                return reject();
        }
        breakdown.statement_g1_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        for (const GT* value : {&statement.d1_initial, &statement.e0, &statement.f0,
                                &statement.t_left0, &statement.t_right0}) {
            ++breakdown.statement_gt_checked;
            if (!check_gt(*value)) return reject();
        }
        breakdown.statement_gt_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        const bool proof_shape = proof.steps.size() == crs.d - 1;
        breakdown.proof_shape_ms = milliseconds(phase, Clock::now());
        if (!proof_shape) return reject();

        phase = Clock::now();
        for (const G1* point : {&proof.phi_final, &proof.r_final}) {
            ++breakdown.g1_elements_checked;
            ++breakdown.proof_g1_checked;
            if (!point->isValid() || !point->isValidOrder() || point->isZero())
                return reject();
        }
        breakdown.proof_g1_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        ++breakdown.g2_elements_checked;
        ++breakdown.proof_g2_checked;
        if (!proof.theta_final.isValid() || !proof.theta_final.isValidOrder()
            || proof.theta_final.isZero()) return reject();
        breakdown.proof_g2_validation_ms = milliseconds(phase, Clock::now());

        phase = Clock::now();
        for (const auto& step : proof.steps) {
            for (const GT* value : {&step.dory_fold.d1_left, &step.dory_fold.d1_right,
                    &step.dory_fold.d2_left, &step.dory_fold.d2_right,
                    &step.dory_fold.w1, &step.dory_fold.w2,
                    &step.rexp_round.e, &step.rexp_round.f,
                    &step.rexp_round.t_left, &step.rexp_round.t_right, &step.u}) {
                ++breakdown.proof_gt_checked;
                if (!check_gt(*value)) return reject();
            }
        }
        breakdown.proof_gt_validation_ms = milliseconds(phase, Clock::now());
        breakdown.total_ms = milliseconds(total_start, Clock::now());
        return ValidatedVerificationInputs(crs, precomputation, statement, proof);
    } catch (...) {
        return reject();
    }
}

} 

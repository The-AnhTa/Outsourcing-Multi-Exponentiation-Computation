#include "symbolic.hpp"

#include "crypto.hpp"

#include <chrono>
#include <map>

namespace rexpbf::internal {
namespace {
using Clock = std::chrono::steady_clock;
double milliseconds(Clock::time_point start, Clock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - start).count();
}
}

GTExpression GTExpression::atom(const GT& base) {
    GTExpression expression;
    expression.multiply_atom(base, fr_one());
    return expression;
}

void GTExpression::multiply(const GTExpression& other) {
    terms_.insert(terms_.end(), other.terms_.begin(), other.terms_.end());
}

void GTExpression::multiply_atom(const GT& base, const Fr& scalar) {
    terms_.push_back({&base, scalar});
}

void GTExpression::scale(const Fr& scalar) {
    for (auto& term : terms_) Fr::mul(term.scalar, term.scalar, scalar);
}

std::size_t GTExpression::normalize() {
    std::map<const GT*, Fr> coalesced;
    for (const auto& term : terms_) {
        const auto [position, inserted] = coalesced.emplace(term.base, term.scalar);
        if (!inserted) Fr::add(position->second, position->second, term.scalar);
    }
    std::size_t zeros = 0;
    terms_.clear();
    for (const auto& [base, scalar] : coalesced) {
        if (scalar.isZero()) ++zeros;
        else terms_.push_back({base, scalar});
    }
    return zeros;
}
SymbolicExpressions build_symbolic_expressions(
    const Precomputation& p,
    const Statement& s,
    const Proof& proof,
    const ChallengeTrace& ch,
    std::span<const Fr> rho_inverse,
    std::span<const Fr> beta_inverse,
    std::span<const Fr> alpha_inverse,
    SymbolicMetrics* metrics) {
    const auto initialization_start = Clock::now();
    auto a0 = GTExpression::atom(s.d1_initial);
    a0.multiply_atom(s.e0, ch.rho[0]);
    a0.multiply_atom(s.f0, rho_inverse[0]);
    auto a1 = GTExpression::atom(s.t_left0);
    a1.multiply_atom(s.t_right0, ch.rho[0]);
    auto outer = a1;
    auto a2 = GTExpression::atom(p.x[1]);
    a2.multiply_atom(p.delta2_right[0], rho_inverse[0]);
    if (metrics)
        metrics->initialization_ms = milliseconds(initialization_start, Clock::now());

    for (std::size_t i = 0; i < proof.steps.size(); ++i) {
        const auto dory_start = Clock::now();
        const std::size_t t = i + 2;
        const std::size_t level = i + 1;
        const auto& step = proof.steps[i];

        auto folded0 = a0;
        folded0.multiply_atom(p.x[level], fr_one());
        auto expression = a1;
        expression.scale(beta_inverse[i]);
        folded0.multiply(expression);
        expression = a2;
        expression.scale(ch.beta[i]);
        folded0.multiply(expression);
        folded0.multiply_atom(step.dory_fold.w1, ch.alpha[i]);
        folded0.multiply_atom(step.dory_fold.w2, alpha_inverse[i]);

        GTExpression folded1;
        folded1.multiply_atom(step.dory_fold.d1_left, ch.alpha[i]);
        folded1.multiply_atom(step.dory_fold.d1_right, fr_one());
        folded1.multiply_atom(p.x[t], fr_multiply(ch.alpha[i], ch.beta[i]));
        folded1.multiply_atom(p.delta1_right[level], ch.beta[i]);

        GTExpression folded2;
        folded2.multiply_atom(step.dory_fold.d2_left, alpha_inverse[i]);
        folded2.multiply_atom(step.dory_fold.d2_right, fr_one());
        folded2.multiply_atom(p.x[t], fr_multiply(alpha_inverse[i], beta_inverse[i]));
        folded2.multiply_atom(p.delta2_right[level], beta_inverse[i]);

        const auto fresh_start = Clock::now();
        if (metrics) metrics->dory_fold_ms += milliseconds(dory_start, fresh_start);

        auto fresh0 = outer;
        fresh0.multiply_atom(step.rexp_round.e, ch.rho[t - 1]);
        fresh0.multiply_atom(step.rexp_round.f, rho_inverse[t - 1]);
        GTExpression fresh1;
        fresh1.multiply_atom(step.rexp_round.t_left, fr_one());
        fresh1.multiply_atom(step.rexp_round.t_right, ch.rho[t - 1]);
        GTExpression fresh2;
        fresh2.multiply_atom(p.x[t], fr_one());
        fresh2.multiply_atom(p.delta2_right[level], rho_inverse[t - 1]);

        const auto batch_start = Clock::now();
        if (metrics) metrics->fresh_rexp_ms += milliseconds(fresh_start, batch_start);
        const Fr gamma_squared = fr_multiply(ch.gamma[i], ch.gamma[i]);
        folded0.scale(gamma_squared);
        folded0.multiply_atom(step.u, ch.gamma[i]);
        folded0.multiply(fresh0);
        folded1.scale(ch.gamma[i]);
        folded1.multiply(fresh1);
        folded2.scale(ch.gamma[i]);
        folded2.multiply(fresh2);

        const auto normalization_start = Clock::now();
        const std::size_t terms_before = folded0.terms().size()
            + folded1.terms().size() + folded2.terms().size()
            + fresh1.terms().size();
        const std::size_t zeros = folded0.normalize() + folded1.normalize()
            + folded2.normalize() + fresh1.normalize();
        const auto normalization_end = Clock::now();
        if (metrics) {
            metrics->batch_ms += milliseconds(batch_start, normalization_start);
            metrics->normalization_ms += milliseconds(normalization_start, normalization_end);
            metrics->normalization_calls += 4;
            metrics->terms_before_normalization += terms_before;
            metrics->terms_after_normalization += folded0.terms().size()
                + folded1.terms().size() + folded2.terms().size()
                + fresh1.terms().size();
            metrics->zero_terms_removed += zeros;
        }

        a0 = std::move(folded0);
        a1 = std::move(folded1);
        a2 = std::move(folded2);
        outer = std::move(fresh1);
    }
    return {std::move(a0), std::move(a1), std::move(a2), std::move(outer)};
}

} // namespace rexpbf::internal

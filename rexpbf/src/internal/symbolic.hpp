#pragma once

#include "rexpbf/verify.hpp"

#include <span>

namespace rexpbf::internal {

struct GTTerm {
    const GT* base;
    Fr scalar;
};

class GTExpression {
public:
    static GTExpression atom(const GT& base);
    void multiply(const GTExpression& other);
    void multiply_atom(const GT& base, const Fr& scalar);
    void scale(const Fr& scalar);
    std::size_t normalize();
    const std::vector<GTTerm>& terms() const { return terms_; }
private:
    std::vector<GTTerm> terms_;
};

struct SymbolicExpressions {
    GTExpression accumulated;
    GTExpression left;
    GTExpression right;
    GTExpression outer;
};

struct SymbolicMetrics {
    double initialization_ms{};
    double dory_fold_ms{};
    double fresh_rexp_ms{};
    double batch_ms{};
    double normalization_ms{};
    std::size_t normalization_calls{};
    std::size_t terms_before_normalization{};
    std::size_t terms_after_normalization{};
    std::size_t zero_terms_removed{};
};

SymbolicExpressions build_symbolic_expressions(
    const Precomputation& precomputation,
    const Statement& statement,
    const Proof& proof,
    const ChallengeTrace& challenges,
    std::span<const Fr> rho_inverses,
    std::span<const Fr> beta_inverses,
    std::span<const Fr> alpha_inverses,
    SymbolicMetrics* metrics);

} 

#include "rexpbf/verify.hpp"

#include "rexpbf/setup.hpp"
#include "internal/crypto.hpp"

namespace rexpbf {
namespace {

bool valid_proof(const CRS& crs, const Proof& proof) {
    if (proof.steps.size() != crs.d - 1
        || !proof.phi_final.isValid() || !proof.phi_final.isValidOrder()
        || !proof.theta_final.isValid() || !proof.theta_final.isValidOrder()
        || !proof.r_final.isValid() || !proof.r_final.isValidOrder()
        || proof.phi_final.isZero() || proof.theta_final.isZero()
        || proof.r_final.isZero()) return false;
    const auto valid_gt = [](const GT& value) { return mcl::bn::isValidGT(value); };
    for (const auto& step : proof.steps) {
        for (const GT* value : {&step.dory_fold.d1_left, &step.dory_fold.d1_right,
                &step.dory_fold.d2_left, &step.dory_fold.d2_right,
                &step.dory_fold.w1, &step.dory_fold.w2,
                &step.rexp_round.e, &step.rexp_round.f,
                &step.rexp_round.t_left, &step.rexp_round.t_right, &step.u})
            if (!valid_gt(*value)) return false;
    }
    return true;
}

bool valid_inputs(const CRS& crs, const Precomputation& precomputation,
                  const Statement& statement, const Proof& proof) {
    return validate_crs(crs)
        && validate_precomputation_shape(crs, precomputation)
        && validate_precomputation_elements(precomputation)
        && validate_statement_shape(crs, statement)
        && validate_statement_elements(statement)
        && validate_statement_digest(crs, statement)
        && valid_proof(crs, proof);
}

} // namespace

bool verify_reference(const CRS& crs, const Precomputation& p,
                      const Statement& s, const Proof& proof) {
    try {
        if (!valid_inputs(crs, p, s, proof) || crs.lambda[0].isZero()) return false;
        const auto ch = replay_challenges(crs, s, proof);
        Fr rho_inverse = internal::fr_inverse_nonzero(ch.rho[0]);
        GT a0 = internal::gt_product({s.d1_initial,
            internal::gt_power(s.e0, ch.rho[0]),
            internal::gt_power(s.f0, rho_inverse)});
        GT a1 = internal::gt_multiply(
            s.t_left0, internal::gt_power(s.t_right0, ch.rho[0]));
        GT outer = a1;
        GT a2 = internal::gt_multiply(
            p.x[1], internal::gt_power(p.delta2_right[0], rho_inverse));

        for (std::size_t i = 0; i < proof.steps.size(); ++i) {
            const auto& step = proof.steps[i];
            const std::size_t t = i + 2;
            const std::size_t level = i + 1;
            const Fr beta_inverse = internal::fr_inverse_nonzero(ch.beta[i]);
            const Fr alpha_inverse = internal::fr_inverse_nonzero(ch.alpha[i]);
            rho_inverse = internal::fr_inverse_nonzero(ch.rho[t - 1]);
            GT folded0 = internal::gt_product({a0, p.x[level],
                internal::gt_power(a1, beta_inverse),
                internal::gt_power(a2, ch.beta[i]),
                internal::gt_power(step.dory_fold.w1, ch.alpha[i]),
                internal::gt_power(step.dory_fold.w2, alpha_inverse)});
            GT folded1 = internal::gt_product({
                internal::gt_power(step.dory_fold.d1_left, ch.alpha[i]),
                step.dory_fold.d1_right,
                internal::gt_power(p.x[t], internal::fr_multiply(ch.alpha[i], ch.beta[i])),
                internal::gt_power(p.delta1_right[level], ch.beta[i])});
            GT folded2 = internal::gt_product({
                internal::gt_power(step.dory_fold.d2_left, alpha_inverse),
                step.dory_fold.d2_right,
                internal::gt_power(p.x[t], internal::fr_multiply(alpha_inverse, beta_inverse)),
                internal::gt_power(p.delta2_right[level], beta_inverse)});
            GT fresh0 = internal::gt_product({outer,
                internal::gt_power(step.rexp_round.e, ch.rho[t - 1]),
                internal::gt_power(step.rexp_round.f, rho_inverse)});
            GT fresh1 = internal::gt_multiply(step.rexp_round.t_left,
                internal::gt_power(step.rexp_round.t_right, ch.rho[t - 1]));
            GT fresh2 = internal::gt_multiply(p.x[t],
                internal::gt_power(p.delta2_right[level], rho_inverse));
            const Fr gamma_squared = internal::fr_multiply(ch.gamma[i], ch.gamma[i]);
            a0 = internal::gt_product({internal::gt_power(folded0, gamma_squared),
                internal::gt_power(step.u, ch.gamma[i]), fresh0});
            a1 = internal::gt_multiply(internal::gt_power(folded1, ch.gamma[i]), fresh1);
            a2 = internal::gt_multiply(internal::gt_power(folded2, ch.gamma[i]), fresh2);
            outer = fresh1;
        }

        const Fr epsilon_inverse = internal::fr_inverse_nonzero(ch.epsilon);
        const GT lhs = internal::gt_product({a0,
            internal::gt_power(a1, epsilon_inverse),
            internal::gt_power(a2, ch.epsilon), p.x[crs.d]});
        const G1 terminal_g1 = internal::g1_add(proof.phi_final,
            internal::g1_multiply(crs.gamma[0], ch.epsilon));
        const G2 terminal_g2 = internal::g2_add(proof.theta_final,
            internal::g2_multiply(crs.lambda[0], epsilon_inverse));
        GT dory_pairing;
        GT rexp_pairing;
        mcl::bn::pairing(dory_pairing, terminal_g1, terminal_g2);
        mcl::bn::pairing(rexp_pairing, proof.r_final, crs.lambda[0]);
        return lhs == dory_pairing && rexp_pairing == outer;
    } catch (...) {
        return false;
    }
}

} // namespace rexpbf

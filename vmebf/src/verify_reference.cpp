#include "vme_ibf/verify_reference.hpp"

#include "vme_ibf/group_utils.hpp"
#include "internal/crypto.hpp"
#include "internal/protocol.hpp"
#include "internal/verification.hpp"

namespace vme_ibf {
namespace {

GT product(std::initializer_list<GT> values) {
    GT result;
    result.setOne();
    for (const GT& value : values) result = gt_mul(result, value);
    return result;
}

} // namespace

ReferenceVerificationTrace verify_reference_diagnostic(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof) {
    ReferenceVerificationTrace trace;
    try {
        if (!internal::validate_verification_objects(
                crs, precomputation, statement, proof, true)) return trace;
        internal::ProtocolChallenges challenges;
        if (!internal::replay_protocol(
                crs, precomputation, statement, proof, challenges)) return trace;
        trace.rho = challenges.rho;
        trace.beta = challenges.beta;
        trace.alpha = challenges.alpha;
        trace.gamma = challenges.gamma;
        trace.epsilon = challenges.epsilon;

        std::vector<DoryTargetState> fresh(crs.d + 1);
        GT outer = precomputation.pairing_x[0];
        for (std::size_t j = 0; j < crs.d; ++j) {
            const RexpClaims claim = j == 0
                ? internal::initial_rexp_claim(precomputation)
                : proof.rexp_claims[j - 1];
            const Fr rho_inverse = inverse_nonzero(challenges.rho[j]);
            fresh[j + 1].D0 = product({outer,
                gt_pow(claim.E, challenges.rho[j]),
                gt_pow(claim.F, rho_inverse)});
            fresh[j + 1].D1 = gt_mul(
                claim.TL, gt_pow(claim.TR, challenges.rho[j]));
            fresh[j + 1].D2 = gt_mul(
                precomputation.pairing_x[j + 1],
                gt_pow(precomputation.delta2R[j], rho_inverse));
            outer = fresh[j + 1].D1;
        }

        DoryTargetState aggregate;
        aggregate.D0 = gt_pow(
            precomputation.pairing_LLprime, challenges.q);
        mcl::bn::pairing(aggregate.D1, crs.L, statement.X);
        mcl::bn::pairing(aggregate.D2, proof.R, crs.Lprime);
        for (std::size_t t = 1; t <= crs.d; ++t) {
            const std::size_t k = t - 1;
            const auto& fold = proof.dory_folds[k];
            const Fr beta_inverse = inverse_nonzero(challenges.beta[k]);
            const Fr alpha_inverse = inverse_nonzero(challenges.alpha[k]);
            DoryTargetState folded;
            folded.D0 = product({aggregate.D0,
                precomputation.pairing_x[k],
                gt_pow(aggregate.D1, beta_inverse),
                gt_pow(aggregate.D2, challenges.beta[k]),
                gt_pow(fold.W1, challenges.alpha[k]),
                gt_pow(fold.W2, alpha_inverse)});
            folded.D1 = product({
                gt_pow(fold.D1L, challenges.alpha[k]), fold.D1R,
                gt_pow(precomputation.pairing_x[k + 1],
                       internal::fr_multiply(
                           challenges.alpha[k], challenges.beta[k])),
                gt_pow(precomputation.delta1R[k], challenges.beta[k])});
            folded.D2 = product({gt_pow(fold.D2L, alpha_inverse), fold.D2R,
                gt_pow(precomputation.pairing_x[k + 1],
                       internal::fr_multiply(alpha_inverse, beta_inverse)),
                gt_pow(precomputation.delta2R[k], beta_inverse)});
            const Fr& gamma = challenges.gamma[t];
            aggregate.D0 = product({
                gt_pow(folded.D0, internal::fr_multiply(gamma, gamma)),
                gt_pow(proof.batch_U[k], gamma), fresh[t].D0});
            aggregate.D1 = gt_mul(gt_pow(folded.D1, gamma), fresh[t].D1);
            aggregate.D2 = gt_mul(gt_pow(folded.D2, gamma), fresh[t].D2);
        }

        const Fr epsilon_inverse = inverse_nonzero(challenges.epsilon);
        const GT lhs = product({aggregate.D0,
            gt_pow(aggregate.D1, epsilon_inverse),
            gt_pow(aggregate.D2, challenges.epsilon),
            precomputation.pairing_x[crs.d]});
        const G1 terminal_g1 = g1_add(
            proof.PhiFinal, g1_pow(crs.G[0], challenges.epsilon));
        const G2 terminal_g2 = g2_add(
            proof.ThetaFinal, g2_pow(crs.H[0], epsilon_inverse));
        GT dory_pairing;
        mcl::bn::pairing(dory_pairing, terminal_g1, terminal_g2);
        GT dory_inverse;
        GT::inv(dory_inverse, dory_pairing);
        trace.dory_residual = gt_mul(lhs, dory_inverse);
        trace.dory_accepted = lhs == dory_pairing;

        GT rexp_pairing;
        mcl::bn::pairing(rexp_pairing, proof.R, crs.H[0]);
        GT outer_inverse;
        GT::inv(outer_inverse, outer);
        trace.rexp_residual = gt_mul(rexp_pairing, outer_inverse);
        trace.rexp_accepted = rexp_pairing == outer;
        trace.final_aggregate = aggregate;
        trace.final_rexp_d1 = outer;
        trace.accepted = trace.dory_accepted && trace.rexp_accepted;
        return trace;
    } catch (...) {
        return {};
    }
}

bool verify_reference(const VmeIbfCRS& crs,
                      const VmeIbfPrecomputation& precomputation,
                      const VmeIbfStatement& statement,
                      const VmeIbfProof& proof) {
    return verify_reference_diagnostic(
        crs, precomputation, statement, proof).accepted;
}

} // namespace vme_ibf

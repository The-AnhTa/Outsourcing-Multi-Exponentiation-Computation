#include "vme_ibf/phase2.hpp"

#include "vme_ibf/group_utils.hpp"
#include "internal/crypto.hpp"
#include "internal/protocol.hpp"

#include <stdexcept>

namespace vme_ibf {
namespace {

GT product(std::initializer_list<GT> values) {
    GT result;
    result.setOne();
    for (const GT& value : values) result = gt_mul(result, value);
    return result;
}

bool valid_claim(const RexpClaims& claim) {
    return internal::valid_gt(claim.E) && internal::valid_gt(claim.F)
        && internal::valid_gt(claim.TL) && internal::valid_gt(claim.TR);
}

bool valid_fresh_instance(const FreshDoryInstance& fresh,
                          std::size_t expected_size) {
    if (fresh.Phi.size() != expected_size
        || fresh.Theta.size() != expected_size
        || !internal::valid_gt(fresh.D0)
        || !internal::valid_gt(fresh.D1)
        || !internal::valid_gt(fresh.D2)) return false;
    for (const auto& point : fresh.Phi)
        if (!valid_g1(point)) return false;
    for (const auto& point : fresh.Theta)
        if (!valid_g2(point)) return false;
    return true;
}

void validate_phase2_inputs(const VmeIbfCRS& crs,
                            const VmeIbfPrecomputation& precomputation,
                            const Phase1Result& phase1) {
    if (!validate_crs(crs)
        || !validate_precomputation_shape(crs, precomputation)
        || !validate_precomputation_elements(precomputation)
        || !audit_precomputation(crs, precomputation)
        || !validate_statement_shape(crs, phase1.statement)
        || !validate_statement_elements(phase1.statement)
        || !validate_statement_digest(crs, phase1.statement)
        || phase1.rho.size() != crs.d || phase1.r.size() != crs.d
        || phase1.dynamic_claims.size() != crs.d - 1
        || phase1.fresh.size() != crs.d + 1 || !valid_g1(phase1.R))
        throw std::invalid_argument("invalid phase-2 input");
    for (const auto& claim : phase1.dynamic_claims)
        if (!valid_claim(claim))
            throw std::invalid_argument("invalid phase-1 claim");
    for (std::size_t level = 1; level <= crs.d; ++level)
        if (!valid_fresh_instance(phase1.fresh[level], crs.n >> level))
            throw std::invalid_argument("invalid fresh Dory instance");
    Digest after_r;
    const auto replayed = replay_rho(crs, precomputation, phase1.statement,
                                     phase1.dynamic_claims, phase1.R, &after_r);
    if (replayed != phase1.rho || after_r != phase1.transcript_after_R)
        throw std::invalid_argument("phase-1 transcript continuation mismatch");
    for (std::size_t j = 0; j < crs.d; ++j)
        if (phase1.r[crs.d - j - 1] != phase1.rho[j])
            throw std::invalid_argument("phase-1 challenge ordering mismatch");
}

DoryInstanceState initial_instance(const VmeIbfCRS& crs,
                                   const VmeIbfPrecomputation& precomputation,
                                   const Phase1Result& phase1,
                                   std::span<const Fr> weights,
                                   const Fr& q) {
    DoryInstanceState instance;
    instance.witness.Phi.reserve(crs.n);
    instance.witness.Theta.reserve(crs.n);
    for (std::size_t i = 0; i < crs.n; ++i) {
        instance.witness.Phi.push_back(g1_pow(crs.L, phase1.statement.x[i]));
        instance.witness.Theta.push_back(g2_pow(crs.Lprime, weights[i]));
    }
    instance.target.D0 = gt_pow(precomputation.pairing_LLprime, q);
    mcl::bn::pairing(instance.target.D1, crs.L, phase1.statement.X);
    mcl::bn::pairing(instance.target.D2, phase1.R, crs.Lprime);
    return instance;
}

DoryFoldProof first_fold_message(const VmeIbfCRS& crs,
                                 const DoryInstanceState& aggregate,
                                 std::size_t half) {
    DoryFoldProof proof;
    const auto phi = std::span(aggregate.witness.Phi);
    const auto theta = std::span(aggregate.witness.Theta);
    const auto g = std::span(crs.G).first(half);
    const auto h = std::span(crs.H).first(half);
    proof.D1L = pairing_product(phi.first(half), h);
    proof.D1R = pairing_product(phi.subspan(half, half), h);
    proof.D2L = pairing_product(g, theta.first(half));
    proof.D2R = pairing_product(g, theta.subspan(half, half));
    return proof;
}

void compressed_witness(const VmeIbfCRS& crs,
                        const DoryInstanceState& aggregate,
                        const Fr& beta, const Fr& beta_inverse,
                        std::vector<G1>& phi, std::vector<G2>& theta) {
    const std::size_t dimension = aggregate.witness.Phi.size();
    phi.clear();
    theta.clear();
    phi.reserve(dimension);
    theta.reserve(dimension);
    for (std::size_t i = 0; i < dimension; ++i) {
        phi.push_back(g1_add(aggregate.witness.Phi[i],
                             g1_pow(crs.G[i], beta)));
        theta.push_back(g2_add(aggregate.witness.Theta[i],
                               g2_pow(crs.H[i], beta_inverse)));
    }
}

DoryInstanceState fold_instance(const VmeIbfPrecomputation& precomputation,
                                const DoryInstanceState& aggregate,
                                DoryFoldProof& proof,
                                std::span<const G1> compressed_phi,
                                std::span<const G2> compressed_theta,
                                std::size_t round,
                                const Fr& beta, const Fr& beta_inverse,
                                const Fr& alpha, const Fr& alpha_inverse) {
    const std::size_t half = compressed_phi.size() / 2;
    proof.W1 = pairing_product(compressed_phi.first(half),
                               compressed_theta.subspan(half, half));
    proof.W2 = pairing_product(compressed_phi.subspan(half, half),
                               compressed_theta.first(half));
    DoryInstanceState folded;
    folded.witness.Phi.reserve(half);
    folded.witness.Theta.reserve(half);
    for (std::size_t i = 0; i < half; ++i) {
        folded.witness.Phi.push_back(g1_add(
            g1_pow(compressed_phi[i], alpha), compressed_phi[half + i]));
        folded.witness.Theta.push_back(g2_add(
            g2_pow(compressed_theta[i], alpha_inverse),
            compressed_theta[half + i]));
    }
    folded.target.D0 = product({aggregate.target.D0,
        precomputation.pairing_x[round],
        gt_pow(aggregate.target.D1, beta_inverse),
        gt_pow(aggregate.target.D2, beta),
        gt_pow(proof.W1, alpha), gt_pow(proof.W2, alpha_inverse)});
    folded.target.D1 = product({gt_pow(proof.D1L, alpha), proof.D1R,
        gt_pow(precomputation.pairing_x[round + 1],
               internal::fr_multiply(alpha, beta)),
        gt_pow(precomputation.delta1R[round], beta)});
    folded.target.D2 = product({gt_pow(proof.D2L, alpha_inverse), proof.D2R,
        gt_pow(precomputation.pairing_x[round + 1],
               internal::fr_multiply(alpha_inverse, beta_inverse)),
        gt_pow(precomputation.delta2R[round], beta_inverse)});
    return folded;
}

DoryInstanceState batch_instances(const DoryInstanceState& folded,
                                  const FreshDoryInstance& fresh,
                                  const GT& u, const Fr& gamma) {
    const Fr gamma_squared = internal::fr_multiply(gamma, gamma);
    DoryInstanceState result;
    result.target.D0 = product({gt_pow(folded.target.D0, gamma_squared),
                                gt_pow(u, gamma), fresh.D0});
    result.target.D1 = gt_mul(gt_pow(folded.target.D1, gamma), fresh.D1);
    result.target.D2 = gt_mul(gt_pow(folded.target.D2, gamma), fresh.D2);
    result.witness.Phi.reserve(fresh.Phi.size());
    result.witness.Theta.reserve(fresh.Theta.size());
    for (std::size_t i = 0; i < fresh.Phi.size(); ++i) {
        result.witness.Phi.push_back(g1_add(
            g1_pow(folded.witness.Phi[i], gamma), fresh.Phi[i]));
        result.witness.Theta.push_back(g2_add(
            g2_pow(folded.witness.Theta[i], gamma), fresh.Theta[i]));
    }
    return result;
}

} // namespace

Phase2Result prove_phase2(const VmeIbfCRS& crs,
                          const VmeIbfPrecomputation& precomputation,
                          const Phase1Result& phase1) {
    validate_phase2_inputs(crs, precomputation, phase1);
    Phase2Result result;
    result.proof.dory_folds.reserve(crs.d);
    result.proof.batch_U.reserve(crs.d);
    result.challenges.beta.resize(crs.d);
    result.challenges.alpha.resize(crs.d);
    result.challenges.gamma.resize(crs.d + 1);
    result.challenges.gamma[0].clear();
    result.folded.resize(crs.d + 1);
    result.aggregate.resize(crs.d + 1);

    Transcript transcript = Transcript::resume(phase1.transcript_after_R);
    const auto weights = tensor_vector(phase1.r);
    if (weights.size() != crs.n) throw std::logic_error("tensor length mismatch");
    result.q = inner_product(weights, phase1.statement.x);
    internal::absorb_vme_initial(transcript, crs.d, crs.n, result.q);
    DoryInstanceState aggregate = initial_instance(
        crs, precomputation, phase1, weights, result.q);
    result.initial_instance = aggregate;
    result.aggregate[0] = aggregate;

    std::vector<G1> compressed_phi;
    std::vector<G2> compressed_theta;
    for (std::size_t t = 1; t <= crs.d; ++t) {
        const std::size_t round = t - 1;
        const std::size_t dimension = crs.n >> round;
        const std::size_t half = dimension / 2;
        if (aggregate.witness.Phi.size() != dimension
            || aggregate.witness.Theta.size() != dimension)
            throw std::logic_error("aggregate dimension mismatch");

        DoryFoldProof fold_proof = first_fold_message(crs, aggregate, half);
        internal::absorb_dory_beta(
            transcript, round, dimension, fold_proof);
        const Fr beta = internal::derive_beta(transcript, round);
        const Fr beta_inverse = inverse_nonzero(beta);
        result.challenges.beta[round] = beta;
        compressed_witness(crs, aggregate, beta, beta_inverse,
                           compressed_phi, compressed_theta);

        fold_proof.W1 = pairing_product(
            std::span(compressed_phi).first(half),
            std::span(compressed_theta).subspan(half, half));
        fold_proof.W2 = pairing_product(
            std::span(compressed_phi).subspan(half, half),
            std::span(compressed_theta).first(half));
        internal::absorb_dory_alpha(
            transcript, round, dimension, fold_proof);
        const Fr alpha = internal::derive_alpha(transcript, round);
        const Fr alpha_inverse = inverse_nonzero(alpha);
        result.challenges.alpha[round] = alpha;
        DoryInstanceState folded = fold_instance(
            precomputation, aggregate, fold_proof, compressed_phi,
            compressed_theta, round, beta, beta_inverse, alpha, alpha_inverse);
        result.folded[t] = folded;
        result.proof.dory_folds.push_back(fold_proof);

        const FreshDoryInstance& fresh = phase1.fresh[t];
        if (fresh.Phi.size() != half || fresh.Theta.size() != half)
            throw std::logic_error("fresh dimension mismatch");
        const GT u = gt_mul(pairing_product(folded.witness.Phi, fresh.Theta),
                            pairing_product(fresh.Phi, folded.witness.Theta));
        result.proof.batch_U.push_back(u);
        internal::absorb_batch_u(transcript, t, half, u);
        const Fr gamma = internal::derive_gamma(transcript, t);
        result.challenges.gamma[t] = gamma;
        aggregate = batch_instances(folded, fresh, u, gamma);
        result.aggregate[t] = aggregate;
    }

    if (aggregate.witness.Phi.size() != 1
        || aggregate.witness.Theta.size() != 1)
        throw std::logic_error("Dory did not terminate");
    result.proof.PhiFinal = aggregate.witness.Phi[0];
    result.proof.ThetaFinal = aggregate.witness.Theta[0];
    internal::absorb_dory_final(
        transcript, result.proof.PhiFinal, result.proof.ThetaFinal);
    result.challenges.epsilon = internal::derive_epsilon(transcript, crs.d);
    result.final_transcript_digest = transcript.digest();
    result.final_aggregate_target = aggregate.target;
    return result;
}

} // namespace vme_ibf

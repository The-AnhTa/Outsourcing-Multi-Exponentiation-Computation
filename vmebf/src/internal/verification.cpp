#include "verification.hpp"

#include "crypto.hpp"
#include "vme_ibf/batch_inversion.hpp"
#include "vme_ibf/setup.hpp"

#include <array>
#include <stdexcept>

namespace vme_ibf::internal {
namespace {

SymbolicGtExpression gt(GtAtomKind kind, std::size_t index = 0) {
    return gt_atom({kind, index});
}

SymbolicGtExpression pairing(PairingAtomKind kind) {
    return pairing_atom({kind});
}

void multiply(SymbolicGtExpression& destination,
              const SymbolicGtExpression& source) {
    multiply_in_place(destination, source);
}

void multiply_power(SymbolicGtExpression& destination,
                    const SymbolicGtExpression& source, const Fr& scalar) {
    multiply(destination, powered(source, scalar));
}

bool valid_proof(const VmeIbfCRS& crs, const VmeIbfProof& proof) {
    if (proof.rexp_claims.size() != crs.d - 1
        || proof.dory_folds.size() != crs.d
        || proof.batch_U.size() != crs.d
        || !valid_g1(proof.R) || !valid_g1(proof.PhiFinal)
        || !valid_g2(proof.ThetaFinal)) return false;
    for (const auto& claim : proof.rexp_claims)
        for (const GT* value : {&claim.E, &claim.F, &claim.TL, &claim.TR})
            if (!valid_gt(*value)) return false;
    for (const auto& fold : proof.dory_folds)
        for (const GT* value : {&fold.D1L, &fold.D1R, &fold.D2L,
                                &fold.D2R, &fold.W1, &fold.W2})
            if (!valid_gt(*value)) return false;
    for (const auto& value : proof.batch_U)
        if (!valid_gt(value)) return false;
    return true;
}

} 

bool validate_verification_objects(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof,
    bool audit_values) {
    try {
        return validate_crs(crs)
            && validate_precomputation_shape(crs, precomputation)
            && validate_precomputation_elements(precomputation)
            && validate_statement_shape(crs, statement)
            && validate_statement_elements(statement)
            && validate_statement_digest(crs, statement)
            && valid_proof(crs, proof)
            && (!audit_values || audit_precomputation(crs, precomputation));
    } catch (...) {
        return false;
    }
}

bool build_verification_equations(
    const VmeIbfCRS& crs, const VmeIbfPrecomputation& precomputation,
    const VmeIbfStatement& statement, const VmeIbfProof& proof,
    VerificationEquations& output) {
    try {
        output = {};
        const auto start = Clock::now();
        if (!replay_protocol(crs, precomputation, statement, proof,
                             output.challenges)) return false;
        const auto after_transcript = Clock::now();

        std::vector<Fr> values;
        values.reserve(3 * crs.d + 1);
        values.insert(values.end(), output.challenges.rho.begin(),
                      output.challenges.rho.end());
        values.insert(values.end(), output.challenges.beta.begin(),
                      output.challenges.beta.end());
        values.insert(values.end(), output.challenges.alpha.begin(),
                      output.challenges.alpha.end());
        values.push_back(output.challenges.epsilon);
        const auto inverses = batch_invert_nonzero(values);
        const auto after_inversion = Clock::now();
        const std::size_t beta_offset = crs.d;
        const std::size_t alpha_offset = 2 * crs.d;
        const std::size_t epsilon_offset = 3 * crs.d;

        std::vector<std::array<SymbolicGtExpression, 3>> fresh(crs.d + 1);
        auto outer = gt(GtAtomKind::PairingX, 0);
        for (std::size_t j = 0; j < crs.d; ++j) {
            const auto e = j ? gt(GtAtomKind::RexpE, j)
                             : gt(GtAtomKind::Delta1R, 0);
            const auto f = j ? gt(GtAtomKind::RexpF, j)
                             : gt(GtAtomKind::Delta2R, 0);
            const auto tl = j ? gt(GtAtomKind::RexpTL, j)
                              : gt(GtAtomKind::PairingX, 1);
            const auto tr = j ? gt(GtAtomKind::RexpTR, j)
                              : gt(GtAtomKind::Delta1R, 0);
            fresh[j + 1][0] = outer;
            multiply_power(fresh[j + 1][0], e, output.challenges.rho[j]);
            multiply_power(fresh[j + 1][0], f, inverses[j]);
            fresh[j + 1][1] = tl;
            multiply_power(fresh[j + 1][1], tr, output.challenges.rho[j]);
            fresh[j + 1][2] = gt(GtAtomKind::PairingX, j + 1);
            multiply_power(fresh[j + 1][2], gt(GtAtomKind::Delta2R, j),
                           inverses[j]);
            outer = fresh[j + 1][1];
        }

        auto a0 = powered(gt(GtAtomKind::PairingLLprime), output.challenges.q);
        auto a1 = pairing(PairingAtomKind::InitialLX);
        auto a2 = pairing(PairingAtomKind::InitialRLprime);
        for (std::size_t t = 1; t <= crs.d; ++t) {
            const std::size_t k = t - 1;
            const Fr& beta = output.challenges.beta[k];
            const Fr& beta_inverse = inverses[beta_offset + k];
            const Fr& alpha = output.challenges.alpha[k];
            const Fr& alpha_inverse = inverses[alpha_offset + k];
            const Fr& gamma = output.challenges.gamma[t];
            auto b0 = a0;
            multiply(b0, gt(GtAtomKind::PairingX, k));
            multiply_power(b0, a1, beta_inverse);
            multiply_power(b0, a2, beta);
            multiply_power(b0, gt(GtAtomKind::DoryW1, k), alpha);
            multiply_power(b0, gt(GtAtomKind::DoryW2, k), alpha_inverse);
            auto b1 = powered(gt(GtAtomKind::DoryD1L, k), alpha);
            multiply(b1, gt(GtAtomKind::DoryD1R, k));
            multiply_power(b1, gt(GtAtomKind::PairingX, k + 1),
                           fr_multiply(alpha, beta));
            multiply_power(b1, gt(GtAtomKind::Delta1R, k), beta);
            auto b2 = powered(gt(GtAtomKind::DoryD2L, k), alpha_inverse);
            multiply(b2, gt(GtAtomKind::DoryD2R, k));
            multiply_power(b2, gt(GtAtomKind::PairingX, k + 1),
                           fr_multiply(alpha_inverse, beta_inverse));
            multiply_power(b2, gt(GtAtomKind::Delta2R, k), beta_inverse);
            a0 = powered(b0, fr_multiply(gamma, gamma));
            multiply_power(a0, gt(GtAtomKind::BatchU, t), gamma);
            multiply(a0, fresh[t][0]);
            a1 = powered(b1, gamma);
            multiply(a1, fresh[t][1]);
            a2 = powered(b2, gamma);
            multiply(a2, fresh[t][2]);
            normalize(a0);
            normalize(a1);
            normalize(a2);
        }
        const auto after_recurrence = Clock::now();

        const Fr& epsilon = output.challenges.epsilon;
        const Fr& epsilon_inverse = inverses[epsilon_offset];
        output.terminal_dory_g1 = g1_add(
            proof.PhiFinal, g1_pow(crs.G[0], epsilon));
        output.terminal_dory_g2 = g2_add(
            proof.ThetaFinal, g2_pow(crs.H[0], epsilon_inverse));
        output.dory = a0;
        multiply_power(output.dory, a1, epsilon_inverse);
        multiply_power(output.dory, a2, epsilon);
        multiply(output.dory, gt(GtAtomKind::PairingX, crs.d));
        multiply_power(output.dory, pairing(PairingAtomKind::TerminalDory),
                       fr_minus_one());
        output.rexp = pairing(PairingAtomKind::TerminalRexp);
        multiply_power(output.rexp, outer, fr_minus_one());
        const auto after_terminal = Clock::now();

        output.timings.transcript_ms = milliseconds(start, after_transcript);
        output.timings.batch_inversion_ms = milliseconds(
            after_transcript, after_inversion);
        output.timings.recurrence_ms = milliseconds(
            after_inversion, after_recurrence);
        output.timings.terminal_assembly_ms = milliseconds(
            after_recurrence, after_terminal);
        return true;
    } catch (...) {
        output = {};
        return false;
    }
}

const GT& resolve_gt_atom(GtAtomId id,
                          const VmeIbfPrecomputation& p,
                          const VmeIbfProof& proof) {
    switch (id.kind) {
    case GtAtomKind::PairingX: return p.pairing_x[id.index];
    case GtAtomKind::Delta1R: return p.delta1R[id.index];
    case GtAtomKind::Delta2R: return p.delta2R[id.index];
    case GtAtomKind::PairingLLprime: return p.pairing_LLprime;
    case GtAtomKind::RexpE: return proof.rexp_claims[id.index - 1].E;
    case GtAtomKind::RexpF: return proof.rexp_claims[id.index - 1].F;
    case GtAtomKind::RexpTL: return proof.rexp_claims[id.index - 1].TL;
    case GtAtomKind::RexpTR: return proof.rexp_claims[id.index - 1].TR;
    case GtAtomKind::DoryD1L: return proof.dory_folds[id.index].D1L;
    case GtAtomKind::DoryD1R: return proof.dory_folds[id.index].D1R;
    case GtAtomKind::DoryD2L: return proof.dory_folds[id.index].D2L;
    case GtAtomKind::DoryD2R: return proof.dory_folds[id.index].D2R;
    case GtAtomKind::DoryW1: return proof.dory_folds[id.index].W1;
    case GtAtomKind::DoryW2: return proof.dory_folds[id.index].W2;
    case GtAtomKind::BatchU: return proof.batch_U[id.index - 1];
    }
    throw std::logic_error("unknown GT atom");
}

PairingInputs resolve_pairing_atom(
    PairingAtomId id, const VmeIbfCRS& crs,
    const VmeIbfStatement& statement, const VmeIbfProof& proof,
    const VerificationEquations& equations) {
    switch (id.kind) {
    case PairingAtomKind::InitialLX: return {crs.L, statement.X};
    case PairingAtomKind::InitialRLprime: return {proof.R, crs.Lprime};
    case PairingAtomKind::TerminalDory:
        return {equations.terminal_dory_g1, equations.terminal_dory_g2};
    case PairingAtomKind::TerminalRexp: return {proof.R, crs.H[0]};
    }
    throw std::logic_error("unknown pairing atom");
}

} 

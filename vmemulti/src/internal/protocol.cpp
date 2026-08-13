#include "protocol.hpp"

#include "vme_ibf/phase1.hpp"

#include <array>

namespace vme_ibf::internal {

RexpClaims initial_rexp_claim(const VmeIbfPrecomputation& p) {
    return {p.delta1R[0], p.delta2R[0], p.pairing_x[1], p.delta1R[0]};
}

void absorb_rexp_claim(Transcript& transcript, std::size_t round,
                       std::size_t dimension, const RexpClaims& claim) {
    std::array<Bytes, 6> fields{encode_u64_be(round), encode_u64_be(dimension),
        serialize(claim.E), serialize(claim.F), serialize(claim.TL),
        serialize(claim.TR)};
    transcript.absorb("VME.BF.G2/REXP-G1-CLAIMS/V2", fields);
}

Fr derive_rho(Transcript& transcript, std::size_t round) {
    return transcript.challenge_nonzero("VME.BF.G2/RHO/V2", round);
}

void absorb_r(Transcript& transcript, const G1& r) {
    transcript.absorb("VME.BF.G2/R-G1/V2", serialize(r));
}

void absorb_vme_initial(Transcript& transcript, std::size_t d,
                        std::size_t n, const Fr& q) {
    std::array<Bytes, 3> fields{encode_u64_be(d), encode_u64_be(n), serialize(q)};
    transcript.absorb("VME.BF.G2/VME-INITIAL/V2", fields);
}

void absorb_dory_beta(Transcript& transcript, std::size_t round,
                      std::size_t dimension, const DoryFoldProof& proof) {
    std::array<Bytes, 6> fields{encode_u64_be(round), encode_u64_be(dimension),
        serialize(proof.D1L), serialize(proof.D1R), serialize(proof.D2L),
        serialize(proof.D2R)};
    transcript.absorb("VME.BF.G2/DORY-BETA-MESSAGE/V2", fields);
}

Fr derive_beta(Transcript& transcript, std::size_t round) {
    return transcript.challenge_nonzero("VME.BF.G2/BETA/V2", round);
}

void absorb_dory_alpha(Transcript& transcript, std::size_t round,
                       std::size_t dimension, const DoryFoldProof& proof) {
    std::array<Bytes, 4> fields{encode_u64_be(round), encode_u64_be(dimension),
        serialize(proof.W1), serialize(proof.W2)};
    transcript.absorb("VME.BF.G2/DORY-ALPHA-MESSAGE/V2", fields);
}

Fr derive_alpha(Transcript& transcript, std::size_t round) {
    return transcript.challenge_nonzero("VME.BF.G2/ALPHA/V2", round);
}

void absorb_batch_u(Transcript& transcript, std::size_t round,
                    std::size_t half, const GT& value) {
    std::array<Bytes, 3> fields{encode_u64_be(round), encode_u64_be(half),
        serialize(value)};
    transcript.absorb("VME.BF.G2/BATCH-U/V2", fields);
}

Fr derive_gamma(Transcript& transcript, std::size_t round) {
    return transcript.challenge_nonzero("VME.BF.G2/GAMMA/V2", round);
}

void absorb_dory_final(Transcript& transcript, const G1& phi, const G2& theta) {
    std::array<Bytes, 2> fields{serialize(phi), serialize(theta)};
    transcript.absorb("VME.BF.G2/DORY-FINAL/V2", fields);
}

Fr derive_epsilon(Transcript& transcript, std::size_t d) {
    return transcript.challenge_nonzero("VME.BF.G2/EPSILON/V2", d);
}

Fr derive_eta(Transcript& transcript, std::size_t d) {
    return transcript.challenge_nonzero("VME.BF.G2/ETA/V2", d);
}

bool replay_protocol(const VmeIbfCRS& crs, const VmeIbfPrecomputation& p,
                     const VmeIbfStatement& statement, const VmeIbfProof& proof,
                     ProtocolChallenges& output, const Digest* start_state) {
    if (proof.rexp_claims.size() != crs.d - 1
        || proof.dory_folds.size() != crs.d || proof.batch_U.size() != crs.d)
        return false;
    Transcript transcript = start_state
        ? Transcript::resume(*start_state) : Transcript(statement.digest);
    output = {};
    output.rho.resize(crs.d);
    for (std::size_t j = 0; j < crs.d; ++j) {
        const RexpClaims claim = j == 0
            ? initial_rexp_claim(p) : proof.rexp_claims[j - 1];
        absorb_rexp_claim(transcript, j, crs.n >> j, claim);
        output.rho[j] = derive_rho(transcript, j);
    }
    absorb_r(transcript, proof.R);
    std::vector<Fr> reversed(crs.d);
    for (std::size_t j = 0; j < crs.d; ++j)
        reversed[crs.d - j - 1] = output.rho[j];
    const auto weights = tensor_vector(reversed);
    if (weights.size() != crs.n) return false;
    output.q = inner_product(weights, statement.x);
    absorb_vme_initial(transcript, crs.d, crs.n, output.q);
    output.beta.resize(crs.d);
    output.alpha.resize(crs.d);
    output.gamma.resize(crs.d + 1);
    output.gamma[0].clear();
    for (std::size_t t = 1; t <= crs.d; ++t) {
        const std::size_t k = t - 1;
        const std::size_t dimension = crs.n >> k;
        absorb_dory_beta(transcript, k, dimension, proof.dory_folds[k]);
        output.beta[k] = derive_beta(transcript, k);
        absorb_dory_alpha(transcript, k, dimension, proof.dory_folds[k]);
        output.alpha[k] = derive_alpha(transcript, k);
        absorb_batch_u(transcript, t, dimension / 2, proof.batch_U[k]);
        output.gamma[t] = derive_gamma(transcript, t);
    }
    absorb_dory_final(transcript, proof.PhiFinal, proof.ThetaFinal);
    output.epsilon = derive_epsilon(transcript, crs.d);
    output.after_epsilon = transcript.digest();
    return true;
}

} // namespace vme_ibf::internal

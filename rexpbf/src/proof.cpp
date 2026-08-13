#include "rexpbf/prove.hpp"

#include "rexpbf/serialization.hpp"
#include "rexpbf/setup.hpp"
#include "rexpbf/transcript.hpp"
#include "internal/protocol.hpp"
#include "internal/proof.hpp"

#include <stdexcept>
#include <utility>

namespace rexpbf {
namespace {

ChallengeTrace replay_challenges_impl(
    const CRS& crs,
    const Statement& statement,
    const Proof& proof,
    TranscriptMetrics* metrics) {
    if (proof.steps.size() != crs.d - 1) {
        throw std::invalid_argument("invalid proof replay input");
    }

    Transcript transcript(internal::transcript_domain, metrics);
    internal::initialize_protocol_transcript(transcript, crs, statement);
    ChallengeTrace trace;
    trace.rho.push_back(internal::derive_initial_rho(transcript, crs));

    for (std::size_t t = 2; t <= crs.d; ++t) {
        const auto& step = proof.steps[t - 2];
        const std::size_t dimension = crs.n >> (t - 1);
        const std::size_t half = dimension / 2;
        const std::size_t round = t - 1;

        internal::append_round_metadata(
            transcript, "REXP-BF-G1-DORY-FIRST-V1", t, round,
            dimension, half);
        transcript.append_gt("D1Left", step.dory_fold.d1_left);
        transcript.append_gt("D1Right", step.dory_fold.d1_right);
        transcript.append_gt("D2Left", step.dory_fold.d2_left);
        transcript.append_gt("D2Right", step.dory_fold.d2_right);
        trace.beta.push_back(transcript.challenge_nonzero_fr(
            "REXP-BF-G1-DORY-BETA-V1", t - 1));

        internal::append_round_metadata(
            transcript, "REXP-BF-G1-DORY-SECOND-V1", t, round,
            dimension, half);
        transcript.append_gt("W1", step.dory_fold.w1);
        transcript.append_gt("W2", step.dory_fold.w2);
        trace.alpha.push_back(transcript.challenge_nonzero_fr(
            "REXP-BF-G1-DORY-ALPHA-V1", t - 1));

        internal::append_round_metadata(
            transcript, "REXP-BF-G1-REXP-ROUND-V1", t, round,
            dimension, half);
        transcript.append_gt("E", step.rexp_round.e);
        transcript.append_gt("F", step.rexp_round.f);
        transcript.append_gt("TLeft", step.rexp_round.t_left);
        transcript.append_gt("TRight", step.rexp_round.t_right);
        trace.rho.push_back(transcript.challenge_nonzero_fr(
            "REXP-BF-G1-RHO-V1", round));

        internal::append_round_metadata(
            transcript, "REXP-BF-G1-BATCH-V1", t, round,
            dimension, half);
        transcript.append_gt("U", step.u);
        trace.gamma.push_back(transcript.challenge_nonzero_fr(
            "REXP-BF-G1-GAMMA-V1", t));
    }

    internal::append_round_metadata(
        transcript, "REXP-BF-G1-DORY-FINAL-V1", crs.d, crs.d - 1, 1, 1);
    transcript.append_g1("phi_final", proof.phi_final);
    transcript.append_g2("theta_final", proof.theta_final);
    trace.epsilon = transcript.challenge_nonzero_fr(
        "REXP-BF-G1-DORY-EPSILON-V1", crs.d);
    transcript.append_g1("R", proof.r_final);
    trace.final_digest = transcript.digest();
    return trace;
}

} // namespace

ChallengeTrace replay_challenges(
    const CRS& crs,
    const Statement& statement,
    const Proof& proof,
    TranscriptMetrics* metrics) {
    if (!validate_crs(crs)
        || !validate_statement_shape(crs, statement)
        || !validate_statement_elements(statement)
        || !validate_statement_digest(crs, statement)) {
        throw std::invalid_argument("invalid proof replay input");
    }
    return replay_challenges_impl(crs, statement, proof, metrics);
}

namespace internal {

ChallengeTrace replay_challenges_prevalidated(
    const CRS& crs,
    const Statement& statement,
    const Proof& proof,
    TranscriptMetrics* metrics) {
    return replay_challenges_impl(crs, statement, proof, metrics);
}

} // namespace internal

std::vector<std::uint8_t> serialize_proof_payload(const Proof& proof) {
    std::vector<std::uint8_t> result;
    const auto append_gt = [&result](const GT& value) {
        const auto bytes = serialize_gt(value);
        result.insert(result.end(), bytes.begin(), bytes.end());
    };
    for (const auto& step : proof.steps) {
        append_gt(step.dory_fold.d1_left);
        append_gt(step.dory_fold.d1_right);
        append_gt(step.dory_fold.d2_left);
        append_gt(step.dory_fold.d2_right);
        append_gt(step.dory_fold.w1);
        append_gt(step.dory_fold.w2);
        append_gt(step.rexp_round.e);
        append_gt(step.rexp_round.f);
        append_gt(step.rexp_round.t_left);
        append_gt(step.rexp_round.t_right);
        append_gt(step.u);
    }
    for (const auto& bytes : {serialize_g1(proof.phi_final),
                              serialize_g2(proof.theta_final),
                              serialize_g1(proof.r_final)}) {
        result.insert(result.end(), bytes.begin(), bytes.end());
    }
    return result;
}

} // namespace rexpbf

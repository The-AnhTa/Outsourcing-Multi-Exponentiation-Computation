#pragma once
#include "rexpbf/types.hpp"

namespace rexpbf {
struct DoryFoldMessage { GT d1_left, d1_right, d2_left, d2_right, w1, w2; };
struct RexpRoundMessage { GT e, f, t_left, t_right; };
struct BatchFoldStep {
    DoryFoldMessage dory_fold;
    RexpRoundMessage rexp_round;
    GT u;
};
struct Proof {
    std::vector<BatchFoldStep> steps;
    G1 phi_final;
    G2 theta_final;
    G1 r_final;
};
struct ProveResult {
    Proof proof;
    std::vector<Fr> r;
};
struct ChallengeTrace {

    std::vector<Fr> rho, beta, alpha, gamma;
    Fr epsilon;
    Digest32 final_digest{};
};
}

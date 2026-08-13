#pragma once

#include "vme_ibf/proof.hpp"
#include "vme_ibf/transcript.hpp"

namespace vme_ibf::internal {

struct ProtocolChallenges {
    std::vector<Fr> rho;
    std::vector<Fr> beta;
    std::vector<Fr> alpha;
    std::vector<Fr> gamma;
    Fr epsilon;
    Fr q;
    Digest after_epsilon{};
};

RexpClaims initial_rexp_claim(const VmeIbfPrecomputation& precomputation);
void absorb_rexp_claim(Transcript&, std::size_t round, std::size_t dimension,
                       const RexpClaims&);
Fr derive_rho(Transcript&, std::size_t round);
void absorb_r(Transcript&, const G1&);
void absorb_vme_initial(Transcript&, std::size_t d, std::size_t n, const Fr& q);
void absorb_dory_beta(Transcript&, std::size_t round, std::size_t dimension,
                      const DoryFoldProof&);
Fr derive_beta(Transcript&, std::size_t round);
void absorb_dory_alpha(Transcript&, std::size_t round, std::size_t dimension,
                       const DoryFoldProof&);
Fr derive_alpha(Transcript&, std::size_t round);
void absorb_batch_u(Transcript&, std::size_t round, std::size_t half,
                    const GT& value);
Fr derive_gamma(Transcript&, std::size_t round);
void absorb_dory_final(Transcript&, const G1&, const G2&);
Fr derive_epsilon(Transcript&, std::size_t d);
Fr derive_eta(Transcript&, std::size_t d);

bool replay_protocol(const VmeIbfCRS&, const VmeIbfPrecomputation&,
                     const VmeIbfStatement&, const VmeIbfProof&,
                     ProtocolChallenges&);

} // namespace vme_ibf::internal

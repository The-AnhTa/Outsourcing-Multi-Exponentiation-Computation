#pragma once
#include "vpip_bf/types.hpp"
#include <span>
#include <string>

namespace vpip_bf {
class ScalarRng { public: virtual ~ScalarRng() = default; virtual Fr random_fr() = 0; };
class DeterministicRng final : public ScalarRng {
 public: explicit DeterministicRng(std::string seed); Fr random_fr() override;
 private: std::string seed_; std::uint64_t counter_{};
};
class SecureRng final : public ScalarRng { public: Fr random_fr() override; };
SetupResult setup(std::size_t d, const std::vector<G2>& public_lambda,
                  const std::vector<G1>& public_X, ScalarRng& rng);
SetupResult setup(std::size_t d, const std::vector<G2>& public_lambda, ScalarRng& rng);
SetupResult setup_parameters(std::size_t d, const std::vector<G2>& public_lambda,
                             const std::vector<G1>& public_X, ScalarRng& rng);
VpipBfPrecomputation precompute(const VpipBfCRS& crs);
Digest compute_crs_digest(const VpipBfCRS& crs);
Digest compute_precomputation_digest(const VpipBfCRS&, const VpipBfPrecomputation&);
Digest compute_statement_input_digest(const VpipBfCRS&, std::span<const G1> X);
}

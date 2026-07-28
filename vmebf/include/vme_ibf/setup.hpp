#pragma once
#include "vme_ibf/types.hpp"
#include <span>
#include <string>

namespace vme_ibf {
class ScalarRng { public: virtual ~ScalarRng() = default; virtual Fr random_fr() = 0; };
class DeterministicRng final : public ScalarRng {
 public: explicit DeterministicRng(std::string seed); Fr random_fr() override;
 private: std::string seed_; std::uint64_t counter_{};
};
class SecureRng final : public ScalarRng { public: Fr random_fr() override; };
SetupResult setup(std::size_t d, const std::vector<G2>& public_H, ScalarRng& rng);
Digest compute_crs_digest(const VmeIbfCRS& crs);
Digest compute_statement_input_digest(const VmeIbfCRS&, std::span<const Fr> x);
}

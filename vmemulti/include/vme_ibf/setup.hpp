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
bool validate_crs_shape(const VmeIbfCRS&);
bool validate_crs_elements(const VmeIbfCRS&);
bool validate_crs_digest(const VmeIbfCRS&);
bool validate_crs(const VmeIbfCRS&);
bool validate_precomputation_shape(const VmeIbfCRS&, const VmeIbfPrecomputation&);
bool validate_precomputation_elements(const VmeIbfPrecomputation&);
bool audit_precomputation(const VmeIbfCRS&, const VmeIbfPrecomputation&);
bool validate_statement_shape(const VmeIbfCRS&, const VmeIbfStatement&);
bool validate_statement_elements(const VmeIbfStatement&);
bool validate_statement_digest(const VmeIbfCRS&, const VmeIbfStatement&);
}

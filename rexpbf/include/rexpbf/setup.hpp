#pragma once
#include "rexpbf/types.hpp"
#include <span>

namespace rexpbf {
struct SetupBreakdown {
    double crs_generation_ms{}, crs_digest_ms{}, precomputation_ms{};
    double h_generation_ms{}, statement_gt_ms{}, statement_digest_ms{}, total_ms{};
};
SetupResult setup(const SetupConfig& config);
SetupResult setup_with_breakdown(const SetupConfig& config, SetupBreakdown& breakdown);
std::span<const G1> gamma_level(const CRS& crs, std::size_t level);
std::span<const G2> lambda_level(const CRS& crs, std::size_t level);
bool validate_crs(const CRS& crs);
bool validate_crs_shape(const CRS& crs);
bool validate_crs_points(const CRS& crs);
bool validate_crs_digest(const CRS& crs);
bool validate_precomputation_shape(const CRS& crs, const Precomputation& precomputation);
bool validate_precomputation_elements(const Precomputation& precomputation);
bool validate_statement_shape(const CRS& crs, const Statement& statement);
bool validate_statement_elements(const Statement& statement);
bool validate_statement_digest(const CRS& crs, const Statement& statement);
bool audit_precomputation(const CRS& crs, const Precomputation& precomputation);
bool audit_statement(const CRS& crs, const Statement& statement);
}

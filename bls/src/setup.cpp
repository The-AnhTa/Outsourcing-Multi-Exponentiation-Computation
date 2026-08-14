#include "blsagg/protocol.hpp"

#include "blsagg/transcript.hpp"
#include "internal/validation.hpp"

#include <limits>
#include <stdexcept>
#include <type_traits>

namespace blsagg {
namespace {

template<class Group>
Group hash_basis_point(std::string_view domain, std::string_view seed,
                       std::size_t index = 0) {
  Bytes input;
  append_frame(input, domain);
  append_frame(input, seed);
  append_frame(input, encode_u64(index));
  Group point;
  if constexpr (std::is_same_v<Group, G1>)
    mcl::bn::hashAndMapToG1(point, input.data(), input.size());
  else
    mcl::bn::hashAndMapToG2(point, input.data(), input.size());
  return point;
}

PublicParameters generate_parameters(std::size_t d, AggregationMode mode,
                                     std::string_view seed) {
  if (d < 1 || d >= std::numeric_limits<std::size_t>::digits)
    throw std::invalid_argument("invalid d");
  if (mode != AggregationMode::BasicDistinct &&
      mode != AggregationMode::Augmented)
    throw std::invalid_argument("invalid aggregation mode");

  PublicParameters parameters;
  parameters.d = d;
  parameters.k = std::size_t{1} << d;
  parameters.mode = mode;
  parameters.H = hash_basis_point<G2>("bls-agg-bf/H", seed);
  parameters.L = hash_basis_point<G1>("bls-agg-bf/L", seed);
  parameters.Lprime = hash_basis_point<G2>("bls-agg-bf/Lprime", seed);
  parameters.Gamma.reserve(parameters.k);
  parameters.Lambda.reserve(parameters.k);
  for (std::size_t i = 0; i < parameters.k; ++i) {
    parameters.Gamma.push_back(
        hash_basis_point<G1>("bls-agg-bf/Gamma", seed, i));
    parameters.Lambda.push_back(
        hash_basis_point<G2>("bls-agg-bf/Lambda", seed, i));
  }
  parameters.digest = internal::parameter_digest(parameters);
  return parameters;
}

void build_prefix_chains(const PublicParameters& parameters,
                         Precomputation& out) {
  out.gamma_chain.push_back(parameters.Gamma);
  out.lambda_chain.push_back(parameters.Lambda);
  for (std::size_t level = 1; level <= parameters.d; ++level) {
    const auto size = parameters.k >> level;
    out.gamma_chain.emplace_back(out.gamma_chain[level - 1].begin(),
                                 out.gamma_chain[level - 1].begin() + size);
    out.lambda_chain.emplace_back(out.lambda_chain[level - 1].begin(),
                                  out.lambda_chain[level - 1].begin() + size);
  }
}

void build_pairing_tables(const PublicParameters& parameters,
                          Precomputation& out) {
  out.X.reserve(parameters.d + 1);
  for (std::size_t level = 0; level <= parameters.d; ++level)
    out.X.push_back(direct_pairing_product(out.gamma_chain[level],
                                           out.lambda_chain[level]));

  for (std::size_t level = 0; level < parameters.d; ++level) {
    const auto size = parameters.k >> level;
    const auto half = size / 2;
    out.delta1L.push_back(direct_pairing_product(
        std::span(out.gamma_chain[level]).first(half),
        out.lambda_chain[level + 1]));
    out.delta1R.push_back(direct_pairing_product(
        std::span(out.gamma_chain[level]).subspan(half, half),
        out.lambda_chain[level + 1]));
    out.delta2L.push_back(direct_pairing_product(
        out.gamma_chain[level + 1],
        std::span(out.lambda_chain[level]).first(half)));
    out.delta2R.push_back(direct_pairing_product(
        out.gamma_chain[level + 1],
        std::span(out.lambda_chain[level]).subspan(half, half)));
  }
}

void build_initial_rexp_claims(const PublicParameters& parameters,
                               Precomputation& out) {
  const auto half = parameters.k / 2;
  out.g1_round0 = {
      direct_pairing_product(std::span(parameters.Gamma).subspan(half, half),
                             std::span(parameters.Lambda).first(half)),
      direct_pairing_product(std::span(parameters.Gamma).first(half),
                             std::span(parameters.Lambda).subspan(half, half)),
      out.delta1L[0], out.delta1R[0]};
  out.g2_round0 = {out.g1_round0.F, out.g1_round0.E, out.delta2L[0],
                   out.delta2R[0]};
}

}  // namespace

SetupResult setup(std::size_t d, AggregationMode mode, std::string_view seed) {
  initialize();
  SetupResult result;
  result.pp = generate_parameters(d, mode, seed);
  result.aux = precompute(result.pp);
  return result;
}

Precomputation precompute(const PublicParameters& parameters) {
  if (!internal::valid_public_parameters(parameters))
    throw std::invalid_argument("invalid public parameters");
  Precomputation result;
  build_prefix_chains(parameters, result);
  build_pairing_tables(parameters, result);
  build_initial_rexp_claims(parameters, result);
  result.digest = internal::precomputation_digest(parameters, result);
  return result;
}

}  // namespace blsagg

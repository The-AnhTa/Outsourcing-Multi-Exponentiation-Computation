#include "vme_ibf/vmemulti.hpp"
#include "vme_ibf/serialization.hpp"
#include "vme_ibf/verify_online.hpp"
#include <array>
#include <chrono>
#include <limits>
#include <stdexcept>

namespace vme_ibf {
namespace {
constexpr std::string_view kDomain = "vmemulti/vme.mb/fs/v1";
constexpr std::string_view kBinding = "VMEMULTI/VME.MB/STATEMENT/V1";
constexpr std::string_view kAggregateChallenge =
    "VMEMULTI/VME.MB/AGGREGATE-GAMMA/V1";

bool canonical_fr(const Fr& value) {
  const Bytes encoded = serialize(value);
  Fr decoded;
  return decoded.deserialize(encoded.data(), encoded.size()) == encoded.size() &&
         decoded == value && serialize(decoded) == encoded;
}

Digest binding_digest_unchecked(const VmeIbfCRS& c,
                                const Digest& precomputation_digest,
                                const VmeMultiStatement& statement) {
  Bytes binding;
  append_frame(binding, kDomain);
  append_frame(binding, "BN254/mcl-v3.00");
  append_frame(binding, c.digest);
  append_frame(binding, precomputation_digest);
  append_frame(binding, encode_u64_be(c.d));
  append_frame(binding, encode_u64_be(c.n));
  append_frame(binding, encode_u64_be(statement.x_instances.size()));
  for (std::size_t i = 0; i < statement.x_instances.size(); ++i) {
    append_frame(binding, encode_u64_be(i));
    for (const auto& scalar : statement.x_instances[i])
      append_frame(binding, serialize(scalar));
    append_frame(binding, serialize(statement.X_instances[i]));
  }
  Bytes framed;
  append_frame(framed, kBinding);
  append_frame(framed, binding);
  return sha256(framed);
}

Transcript bound_transcript(const VmeIbfCRS& c,
                            const VmeIbfPrecomputation& p,
                            const VmeMultiStatement& statement) {
  return Transcript(compute_vmemulti_binding_digest(c, p, statement));
}

std::pair<std::vector<Fr>, std::vector<Fr>> aggregate_z_and_scalars(
    const std::vector<std::vector<Fr>>& x_instances, const Fr& gamma) {
  const std::size_t n = x_instances.front().size();
  std::vector<Fr> z(n), scalars(x_instances.size());
  for (auto& scalar : z) scalar.clear();
  Fr rho;
  rho = 1;
  for (std::size_t i = 0; i < x_instances.size(); ++i) {
    if (x_instances[i].size() != n)
      throw std::invalid_argument("aggregation dimension mismatch");
    scalars[i] = rho;
    for (std::size_t j = 0; j < n; ++j) {
      Fr term;
      Fr::mul(term, rho, x_instances[i][j]);
      Fr::add(z[j], z[j], term);
    }
    Fr::mul(rho, rho, gamma);
  }
  return {std::move(z), std::move(scalars)};
}
}

bool validate_vmemulti_statement(const VmeIbfCRS& c,
                                 const VmeMultiStatement& statement) {
  if (c.d < 1 || c.d >= std::numeric_limits<std::size_t>::digits ||
      c.n != (std::size_t{1} << c.d) || c.G.size() != c.n ||
      c.H.size() != c.n || statement.x_instances.empty() ||
      statement.x_instances.size() != statement.X_instances.size())
    return false;
  for (const auto& x : statement.x_instances) {
    if (x.size() != c.n) return false;
    for (const auto& scalar : x)
      if (!canonical_fr(scalar)) return false;
  }
  for (const auto& X : statement.X_instances)
    if (!valid_g2(X)) return false;
  return true;
}

Digest compute_vmemulti_binding_digest(const VmeIbfCRS& c,
                                       const VmeIbfPrecomputation& p,
                                       const VmeMultiStatement& statement) {
  if (!validate_vmemulti_statement(c, statement) ||
      compute_crs_digest(c) != c.digest || !validate_precomputation(c, p))
    throw std::invalid_argument("invalid vmemulti statement, CRS, or precomputation");
  return binding_digest_unchecked(c, sha256(serialize_precomputation(c, p)),
                                  statement);
}

Fr derive_vmemulti_aggregate_gamma(const VmeIbfCRS& c,
                                   const VmeIbfPrecomputation& p,
                                   const VmeMultiStatement& statement,
                                   Digest* transcript_after_gamma) {
  Transcript transcript = bound_transcript(c, p, statement);
  Fr gamma = transcript.challenge_nonzero(kAggregateChallenge, 0);
  if (transcript_after_gamma) *transcript_after_gamma = transcript.digest();
  return gamma;
}

AggregatedInstance aggregate_instances(
    const std::vector<std::vector<Fr>>& x_instances,
    const std::vector<G2>& X_instances, const Fr& gamma) {
  if (x_instances.empty() || x_instances.size() != X_instances.size() ||
      gamma.isZero())
    throw std::invalid_argument("invalid aggregation input");
  const std::size_t n = x_instances.front().size();
  if (n == 0) throw std::invalid_argument("empty exponent vector");
  AggregatedInstance out;
  auto [z, scalars] = aggregate_z_and_scalars(x_instances, gamma);
  out.z = std::move(z);
  out.Y = g2_multiexp(X_instances, scalars);
  return out;
}

VmeMultiProverTrace prove_vmemulti_with_trace(
    const VmeIbfCRS& c, const VmeIbfPrecomputation& p,
    const VmeMultiStatement& statement) {
  Digest after_gamma;
  VmeMultiProverTrace out;
  out.aggregate_gamma =
      derive_vmemulti_aggregate_gamma(c, p, statement, &after_gamma);
  out.aggregate = aggregate_instances(statement.x_instances,
                                      statement.X_instances,
                                      out.aggregate_gamma);
  out.phase1 = prove_phase1_core(c, p, out.aggregate.z, out.aggregate.Y,
                                 Transcript::resume(after_gamma));
  out.phase2 = prove_phase2(c, p, out.phase1);
  return out;
}

VmeIbfProof prove_vmemulti(const VmeIbfCRS& c,
                           const VmeIbfPrecomputation& p,
                           const VmeMultiStatement& statement) {
  auto trace = prove_vmemulti_with_trace(c, p, statement);
  return assemble_public_proof(trace.phase1, trace.phase2);
}

ReferenceVerificationTrace verify_vmemulti_diagnostic(
    const VmeIbfCRS& c, const VmeIbfPrecomputation& p,
    const VmeMultiStatement& statement, const VmeIbfProof& proof) {
  try {
    Digest after_gamma;
    Fr gamma = derive_vmemulti_aggregate_gamma(c, p, statement, &after_gamma);
    auto aggregate =
        aggregate_instances(statement.x_instances, statement.X_instances, gamma);
    VmeIbfStatement core_statement;
    core_statement.x = std::move(aggregate.z);
    core_statement.X = aggregate.Y;
    VmeIbfStatementInput input{core_statement.x, {}};
    input.digest = compute_statement_input_digest(c, input.x);
    core_statement.digest = compute_statement_digest(c, input, core_statement.X);
    return verify_reference_core_diagnostic(
        c, p, core_statement, proof, Transcript::resume(after_gamma));
  } catch (...) {
    return {};
  }
}

bool verify_vmemulti(const VmeIbfCRS& c, const VmeIbfPrecomputation& p,
                     const VmeMultiStatement& statement,
                     const VmeIbfProof& proof) {
  ValidatedVmeMultiInputs inputs;
  return prepare_validated_vmemulti_inputs(c, p, statement, proof, inputs) &&
         verify_vmemulti_online(inputs);
}

CombinedVerificationTrace verify_vmemulti_combined_diagnostic(
    const VmeIbfCRS& c, const VmeIbfPrecomputation& p,
    const VmeMultiStatement& statement, const VmeIbfProof& proof) {
  try {
    Digest after_gamma;
    const Fr gamma =
        derive_vmemulti_aggregate_gamma(c, p, statement, &after_gamma);
    auto aggregate =
        aggregate_instances(statement.x_instances, statement.X_instances, gamma);
    VmeIbfStatement core_statement;
    core_statement.x = std::move(aggregate.z);
    core_statement.X = aggregate.Y;
    VmeIbfStatementInput input{core_statement.x, {}};
    input.digest = compute_statement_input_digest(c, input.x);
    core_statement.digest = compute_statement_digest(c, input, core_statement.X);
    core_statement.transcript_state = after_gamma;
    core_statement.has_transcript_state = true;
    return verify_deferred_combined_with_trace(c, p, core_statement, proof);
  } catch (...) {
    return {};
  }
}

bool prepare_validated_vmemulti_inputs(
    const VmeIbfCRS& c, const VmeIbfPrecomputation& p,
    const VmeMultiStatement& statement, const VmeIbfProof& proof,
    ValidatedVmeMultiInputs& out) {
  out = {};
  try {
    if (!validate_vmemulti_statement(c, statement)) return false;
    const Digest precomputation_digest =
        sha256(serialize_precomputation(c, p));
    const Digest binding =
        binding_digest_unchecked(c, precomputation_digest, statement);
    Transcript transcript(binding);
    const Fr gamma = transcript.challenge_nonzero(kAggregateChallenge, 0);
    auto aggregate =
        aggregate_instances(statement.x_instances, statement.X_instances, gamma);
    VmeIbfStatement core_statement;
    core_statement.x = std::move(aggregate.z);
    core_statement.X = aggregate.Y;
    VmeIbfStatementInput input{core_statement.x, {}};
    input.digest = compute_statement_input_digest(c, input.x);
    core_statement.digest = compute_statement_digest(c, input, core_statement.X);
    ValidatedVerificationInputs core_inputs;
    if (!prepare_validated_verification_inputs(c, p, core_statement, proof,
                                               core_inputs))
      return false;
    out = {&c, &p, &statement, &proof, precomputation_digest};
    return true;
  } catch (...) {
    out = {};
    return false;
  }
}

VmeMultiOnlineVerificationTrace verify_vmemulti_online_with_trace(
    const ValidatedVmeMultiInputs& inputs) {
  using Clock = std::chrono::steady_clock;
  auto ms = [](auto a, auto b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
  };
  VmeMultiOnlineVerificationTrace out;
  const auto total_start = Clock::now();
  try {
    if (!inputs.crs || !inputs.precomputation || !inputs.statement ||
        !inputs.proof)
      return out;
    const auto& c = *inputs.crs;
    const auto& statement = *inputs.statement;
    auto start = Clock::now();
    const Digest binding = binding_digest_unchecked(
        c, inputs.precomputation_binding_digest, statement);
    auto end = Clock::now();
    out.transcript_prefix_ms = ms(start, end);

    start = Clock::now();
    Transcript transcript(binding);
    const Fr gamma = transcript.challenge_nonzero(kAggregateChallenge, 0);
    const Digest after_gamma = transcript.digest();
    end = Clock::now();
    out.aggregate_gamma_ms = ms(start, end);

    start = Clock::now();
    auto [z, scalars] = aggregate_z_and_scalars(statement.x_instances, gamma);
    end = Clock::now();
    out.aggregate_z_ms = ms(start, end);

    start = Clock::now();
    G2 Y = g2_multiexp(statement.X_instances, scalars);
    end = Clock::now();
    out.aggregate_Y_ms = ms(start, end);

    start = Clock::now();
    VmeIbfStatement core_statement;
    core_statement.x = std::move(z);
    core_statement.X = Y;
    core_statement.transcript_state = after_gamma;
    core_statement.has_transcript_state = true;
    ValidatedVerificationInputs core_inputs{
        &c, inputs.precomputation, &core_statement, inputs.proof};
    end = Clock::now();
    out.core_input_construction_ms = ms(start, end);
    out.core = verify_online_with_trace(core_inputs);
    out.accepted = out.core.accepted;
    out.total_ms = ms(total_start, Clock::now());
    return out;
  } catch (...) {
    out.total_ms = ms(total_start, Clock::now());
    return out;
  }
}

bool verify_vmemulti_online(const ValidatedVmeMultiInputs& inputs) {
  return verify_vmemulti_online_with_trace(inputs).accepted;
}

}

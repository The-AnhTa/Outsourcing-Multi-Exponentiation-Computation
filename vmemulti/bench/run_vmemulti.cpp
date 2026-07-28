#include "vme_ibf/serialization.hpp"
#include "vme_ibf/vmemulti.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace vme_ibf;

namespace {
std::size_t number(const char* text, const char* name) {
  const std::string value(text);
  if (value.empty() || value[0] == '-')
    throw std::invalid_argument(std::string("invalid ") + name);
  std::size_t parsed = 0;
  const auto result = std::stoull(value, &parsed);
  if (parsed != value.size() ||
      result > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument(std::string("invalid ") + name);
  return static_cast<std::size_t>(result);
}

std::vector<G2> public_h(std::size_t n, std::size_t d) {
  std::vector<G2> out;
  out.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto tag = "vmemulti/single-run/H/" + std::to_string(d) + "/" +
                     std::to_string(i);
    G2 point;
    mcl::bn::hashAndMapToG2(point, tag.data(), tag.size());
    out.push_back(point);
  }
  return out;
}

VmeMultiStatement make_statement(const VmeIbfCRS& crs, std::size_t m,
                                 DeterministicRng& rng) {
  VmeMultiStatement out;
  out.x_instances.reserve(m);
  out.X_instances.reserve(m);
  for (std::size_t i = 0; i < m; ++i) {
    std::vector<Fr> x;
    x.reserve(crs.n);
    for (std::size_t j = 0; j < crs.n; ++j)
      x.push_back(rng.random_fr());
    out.X_instances.push_back(g2_multiexp(crs.H, x));
    out.x_instances.push_back(std::move(x));
  }
  return out;
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 3)
      throw std::invalid_argument("expected dimension and instance count");

    initialize();
    const std::size_t d = number(argv[1], "d");
    const std::size_t m = number(argv[2], "m");
    if (d < 1 || d > 12 || m < 1)
      throw std::invalid_argument("unsupported parameters");
    const std::size_t n = std::size_t{1} << d;

    DeterministicRng setup_rng("vmemulti/single-run/setup/d=" +
                               std::to_string(d));
    auto fixture = setup(d, public_h(n, d), setup_rng);
    DeterministicRng statement_rng(
        "vmemulti/single-run/statement/d=" + std::to_string(d) + "/m=" +
        std::to_string(m));
    const auto statement =
        make_statement(fixture.crs, m, statement_rng);
    const auto proof =
        prove_vmemulti(fixture.crs, fixture.precomp, statement);

    ValidatedVmeMultiInputs validated;
    if (!prepare_validated_vmemulti_inputs(
            fixture.crs, fixture.precomp, statement, proof, validated))
      throw std::runtime_error("input validation failed");

    const auto verify_start = std::chrono::steady_clock::now();
    const bool accepted = verify_vmemulti_online(validated);
    const auto verify_end = std::chrono::steady_clock::now();
    if (!accepted) throw std::runtime_error("verification failed");

    const double verification_ms =
        std::chrono::duration<double, std::milli>(
            verify_end - verify_start)
            .count();
    const auto proof_bytes = serialize_proof(d, proof).size();
    const auto crs_bytes = serialize_crs(fixture.crs).size();

    std::cout << std::fixed << std::setprecision(3)
              << "Verification time: " << verification_ms << " ms\n"
              << "Proof size: " << proof_bytes << " bytes\n"
              << "CRS size: " << crs_bytes << " bytes\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

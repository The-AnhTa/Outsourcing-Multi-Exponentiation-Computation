#include "vpip_bf/group_utils.hpp"
#include "vpip_bf/protocol.hpp"
#include "vpip_bf/serialization.hpp"
#include "vpip_bf/verify_reference.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace vpip_bf;
using Clock = std::chrono::steady_clock;

namespace {
std::size_t dimension(const char* text) {
  const std::string value(text);
  if (value.empty() || value[0] == '-')
    throw std::invalid_argument("invalid dimension");
  std::size_t parsed = 0;
  const auto result = std::stoull(value, &parsed);
  if (parsed != value.size() || result < 1 || result > 12 ||
      result >= std::numeric_limits<std::size_t>::digits)
    throw std::invalid_argument("dimension must be in [1,12]");
  return static_cast<std::size_t>(result);
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 2)
      throw std::invalid_argument("expected one dimension");
    const std::size_t d = dimension(argv[1]);
    const std::size_t n = std::size_t{1} << d;

    initialize();
    const std::string lambda_tag =
        "vpipbf-single-run-lambda/d=" + std::to_string(d);
    G2 lambda_base;
    mcl::bn::hashAndMapToG2(
        lambda_base, lambda_tag.data(), lambda_tag.size());
    std::vector<G2> lambda;
    lambda.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      Fr scalar;
      scalar = static_cast<int>(i + 1);
      lambda.push_back(g2_pow(lambda_base, scalar));
    }

    const std::string statement_tag =
        "vpipbf-single-run-X/d=" + std::to_string(d);
    G1 statement_base;
    mcl::bn::hashAndMapToG1(
        statement_base, statement_tag.data(), statement_tag.size());
    std::vector<G1> x;
    x.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
      Fr scalar;
      scalar = static_cast<int>(i + 3);
      x.push_back(g1_pow(statement_base, scalar));
    }

    DeterministicRng rng(
        "vpipbf-single-run/d=" + std::to_string(d));
    auto setup = Setup(d, lambda, x, rng);
    auto precomputation = Precompute(setup.crs);
    auto result =
        Prove(setup.crs, precomputation, setup.statement_input);
    const auto proof_size = serialize_proof(d, result.proof).size();
    const auto crs_size = serialize_crs(setup.crs).size();

    const auto direct_start = Clock::now();
    const GT direct_pairing_product =
        pairing_product(setup.statement_input.X, setup.crs.H);
    const auto direct_stop = Clock::now();
    if (direct_pairing_product != result.statement.C)
      throw std::runtime_error("direct pairing product mismatch");

    auto validated = PrevalidateVerificationInputs(
        setup.crs, precomputation, result.statement, result.proof);
    if (!validated)
      throw std::runtime_error("input validation failed");

    reset_verification_core_call_count_for_testing();
    const auto verify_start = Clock::now();
    const bool accepted = VerifyOnline(*validated);
    const auto verify_stop = Clock::now();
    if (!accepted || verification_core_call_count_for_testing() != 1)
      throw std::runtime_error("verification failed");

    const double verification_ms =
        std::chrono::duration<double, std::milli>(
            verify_stop - verify_start)
            .count();
    const double direct_pairing_ms =
        std::chrono::duration<double, std::milli>(
            direct_stop - direct_start)
            .count();

    std::cout << std::fixed << std::setprecision(3)
              << "Verification time: " << verification_ms << " ms\n"
              << "Proof size: " << proof_size << " bytes\n"
              << "CRS size: " << crs_size << " bytes\n"
              << "Direct pairing product time: " << direct_pairing_ms
              << " ms\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

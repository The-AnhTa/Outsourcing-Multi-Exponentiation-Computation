#include "bp/bp.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace bp;
using Clock = std::chrono::steady_clock;

namespace {
std::size_t dimension(const char* text) {
  const std::string value(text);
  if (value.empty() || value[0] == '-')
    throw std::invalid_argument("invalid dimension");
  std::size_t parsed = 0;
  const auto result = std::stoull(value, &parsed);
  if (parsed != value.size() || result < 1 || result > 20 ||
      result >= std::numeric_limits<std::size_t>::digits)
    throw std::invalid_argument("dimension must be in [1,20]");
  return static_cast<std::size_t>(result);
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 2)
      throw std::invalid_argument("expected one dimension");
    const std::size_t d = dimension(argv[1]);
    const std::size_t n = std::size_t{1} << d;

    const std::string seed_text =
        "bp-single-run/d=" + std::to_string(d);
    const Bytes seed(seed_text.begin(), seed_text.end());
    const PublicParams parameters = Setup(n, seed);

    std::vector<Scalar> x(n), y(n);
    for (std::size_t i = 0; i < n; ++i) {
      x[i].setByCSPRNG();
      y[i].setByCSPRNG();
    }
    const Group statement = commit(parameters, x, y);

    const auto prove_start = Clock::now();
    const Proof proof = Prove(parameters, statement, x, y);
    const auto prove_stop = Clock::now();

    const auto verify_start = Clock::now();
    const bool accepted = Verify(parameters, statement, proof);
    const auto verify_stop = Clock::now();
    if (!accepted)
      throw std::runtime_error("verification failed");

    const auto proof_size =
        serialize_proof(parameters, proof).size();
    std::size_t crs_size = serialize(parameters.K).size();
    for (const auto& point : parameters.G)
      crs_size += serialize(point).size();
    for (const auto& point : parameters.H)
      crs_size += serialize(point).size();

    const double prover_ms =
        std::chrono::duration<double, std::milli>(
            prove_stop - prove_start)
            .count();
    const double verification_ms =
        std::chrono::duration<double, std::milli>(
            verify_stop - verify_start)
            .count();

    std::cout << std::fixed << std::setprecision(3)
              << "Verification time: " << verification_ms << " ms\n"
              << "Prover time: " << prover_ms << " ms\n"
              << "Proof size: " << proof_size << " bytes\n"
              << "CRS size: " << crs_size << " bytes\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

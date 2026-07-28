#include "pinkas.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

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

pinkas::Scalar random_nonzero_scalar() {
  pinkas::Scalar scalar;
  do {
    scalar.setByCSPRNG();
  } while (scalar.isZero());
  return scalar;
}

std::vector<pinkas::Group> random_bases(std::size_t n) {
  pinkas::Group generator;
  constexpr char domain[] = "pinkas-single-run-bn254-g2-generator";
  mcl::bn::hashAndMapToG2(generator, domain, sizeof(domain) - 1);

  std::vector<pinkas::Group> bases(n);
  for (auto& base : bases)
    pinkas::Group::mul(base, generator, random_nonzero_scalar());
  return bases;
}

pinkas::ScalarMatrix random_scalars(std::size_t k, std::size_t n) {
  pinkas::ScalarMatrix scalars(k, std::vector<pinkas::Scalar>(n));
  for (auto& row : scalars)
    for (auto& scalar : row)
      scalar.setByCSPRNG();
  return scalars;
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 3)
      throw std::invalid_argument("expected dimension and instance count");

    const std::size_t d = number(argv[1], "d");
    const std::size_t k = number(argv[2], "k");
    if (d < 1 || d > 30 || k < 1 ||
        d >= std::numeric_limits<std::size_t>::digits)
      throw std::invalid_argument("unsupported parameters");

    constexpr std::size_t window_width = 8;
    const std::size_t n = std::size_t{1} << d;
    const auto parameters = pinkas::setup(window_width);
    const auto bases = random_bases(n);
    const auto scalars = random_scalars(k, n);
    const auto prover_result = pinkas::prove(parameters, bases, scalars);
    const auto encoded_proof =
        pinkas::serialize_proof(parameters, prover_result.proof);

    pinkas::ValidatedInputs validated;
    if (!pinkas::prepare_validated_inputs(
            parameters, bases, scalars, prover_result.Y,
            prover_result.proof, validated))
      throw std::runtime_error("input validation failed");

    const auto verify_start = std::chrono::steady_clock::now();
    const bool accepted = pinkas::verify_online_prevalidated(validated);
    const auto verify_end = std::chrono::steady_clock::now();
    if (!accepted) throw std::runtime_error("verification failed");

    const double verification_ms =
        std::chrono::duration<double, std::milli>(
            verify_end - verify_start)
            .count();
    std::cout << std::fixed << std::setprecision(3)
              << "Verification time: " << verification_ms << " ms\n"
              << "Proof size: " << encoded_proof.size() << " bytes\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

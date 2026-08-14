#include "pinkas/pinkas.hpp"

#include "benchmark_utils.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <limits>

int main(int argc, char** argv) {
  try {
    if (argc != 3)
      throw std::invalid_argument("expected dimension and instance count");

    const std::size_t d = benchmark_utils::number(argv[1], "d");
    const std::size_t k = benchmark_utils::number(argv[2], "k");
    if (d < 1 || d > 30 || k < 1 ||
        d >= std::numeric_limits<std::size_t>::digits)
      throw std::invalid_argument("unsupported parameters");

    constexpr std::size_t window_width = 8;
    const std::size_t n = std::size_t{1} << d;
    const auto parameters = pinkas::setup(window_width);
    const auto bases = benchmark_utils::random_bases(
        n, "pinkas-single-run-bn254-g2-generator");
    const auto scalars = benchmark_utils::random_scalars(k, n);
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

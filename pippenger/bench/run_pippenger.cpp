#include "pippenger/pippenger.hpp"

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
    pippenger::initialize();
    const std::size_t n = std::size_t{1} << d;
    const auto bases = benchmark_utils::random_bases(
        n, "pippenger-single-run-bn254-g2-generator");
    const auto scalars = benchmark_utils::random_scalars(k, n);

    const auto start = std::chrono::steady_clock::now();
    const auto result =
        pippenger::multi_pippenger(bases, scalars, window_width);
    const auto end = std::chrono::steady_clock::now();
    if (result.size() != k)
      throw std::runtime_error("unexpected result count");

    const double running_ms =
        std::chrono::duration<double, std::milli>(end - start).count();
    std::cout << std::fixed << std::setprecision(3)
              << "Running time: " << running_ms << " ms\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

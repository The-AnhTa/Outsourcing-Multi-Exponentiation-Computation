#include "pippenger.hpp"

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

pippenger::Scalar random_nonzero_scalar() {
  pippenger::Scalar scalar;
  do {
    scalar.setByCSPRNG();
  } while (scalar.isZero());
  return scalar;
}

std::vector<pippenger::Group> random_bases(std::size_t n) {
  pippenger::Group generator;
  constexpr char domain[] = "pippenger-single-run-bn254-g2-generator";
  mcl::bn::hashAndMapToG2(generator, domain, sizeof(domain) - 1);

  std::vector<pippenger::Group> bases(n);
  for (auto& base : bases)
    pippenger::Group::mul(base, generator, random_nonzero_scalar());
  return bases;
}

pippenger::ScalarMatrix random_scalars(std::size_t k, std::size_t n) {
  pippenger::ScalarMatrix scalars(
      k, std::vector<pippenger::Scalar>(n));
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
    pippenger::initialize();
    const std::size_t n = std::size_t{1} << d;
    const auto bases = random_bases(n);
    const auto scalars = random_scalars(k, n);

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

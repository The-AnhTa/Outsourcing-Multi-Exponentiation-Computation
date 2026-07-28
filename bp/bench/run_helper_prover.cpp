#include "bp/helper_prover.hpp"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace bp;

namespace {

std::size_t parse_dimension(const char* text) {
  const std::string value(text);
  if (value.empty() || value.front() == '-')
    throw std::invalid_argument("invalid dimension");
  std::size_t parsed = 0;
  const auto result = std::stoull(value, &parsed);
  if (parsed != value.size() || result < 1 || result > 20 ||
      result >= std::numeric_limits<std::size_t>::digits)
    throw std::invalid_argument("dimension must be in [1,20]");
  return static_cast<std::size_t>(result);
}

std::size_t g1_bytes(const G1& point) {
  unsigned char encoded[256]{};
  const std::size_t size = point.serialize(encoded, sizeof(encoded));
  if (size == 0) throw std::runtime_error("G1 serialization failed");
  return size;
}

std::size_t crs_payload_bytes(const HpPublicParams& parameters) {


  std::size_t size = serialize(parameters.bp.K).size();
  for (const auto& point : parameters.bp.G) size += serialize(point).size();
  for (const auto& point : parameters.bp.H) size += serialize(point).size();
  for (const auto& point : parameters.vme.auxiliary_G)
    size += g1_bytes(point);
  size += g1_bytes(parameters.vme.L);
  size += serialize(parameters.vme.Lprime).size();
  return size;
}

}

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::invalid_argument("expected one dimension");
    const std::size_t d = parse_dimension(argv[1]);
    const std::size_t n = std::size_t{1} << d;

    const std::string seed_text =
        "helper-prover-single-run/d=" + std::to_string(d);
    const Bytes seed(seed_text.begin(), seed_text.end());
    const HpSetupResult setup = setup_hp(n, seed);
    const HpInstance instance = generate_hp_instance(setup.pp);

    const HpProof helper =
        prove_hp(setup.pp, setup.helper_precomp, instance.Z, instance.x,
                 instance.y);
    const Bytes wire = serialize_hp_proof(setup.pp, helper);

    HpVerifyTimings timings;
    const auto result =
        verify_hp(setup.pp, setup.client_precomp, instance.Z, instance.x,
                  instance.y, helper, &timings);
    if (!result || !Verify(setup.pp.bp, instance.Z, *result))
      throw std::runtime_error("helper-prover protocol failed");

    const double prover_ms =
        timings.total_ms - timings.parse_and_validation_ms;
    const std::size_t crs_size = crs_payload_bytes(setup.pp);

    std::cout << std::fixed << std::setprecision(3)
              << "Prover time: " << prover_ms << " ms\n"
              << "Proof size: " << wire.size() << " bytes\n"
              << "CRS size: " << crs_size << " bytes\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

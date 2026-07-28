#include "blsagg/protocol.hpp"
#include "blsagg/serialization.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace blsagg;
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

Statement make_statement(const PublicParameters& parameters) {
  Statement statement;
  std::vector<Fr> secret_keys;
  secret_keys.reserve(parameters.k);
  for (std::size_t i = 0; i < parameters.k; ++i) {
    const auto message =
        "bls-single-run-message-" + std::to_string(i);
    statement.messages.emplace_back(message.begin(), message.end());
    Fr secret;
    secret = static_cast<int>(i + 101);
    secret_keys.push_back(secret);
    G2 public_key;
    G2::mul(public_key, parameters.H, secret);
    statement.public_keys.push_back(public_key);
  }

  const auto message_points = hash_messages(parameters, statement);
  statement.sigma_agg.clear();
  for (std::size_t i = 0; i < parameters.k; ++i) {
    G1 signature;
    G1::mul(signature, message_points[i], secret_keys[i]);
    G1::add(
        statement.sigma_agg, statement.sigma_agg, signature);
  }
  return statement;
}
}

int main(int argc, char** argv) {
  try {
    if (argc != 2)
      throw std::invalid_argument("expected one dimension");
    const std::size_t d = dimension(argv[1]);

    initialize();
    const auto setup_result = setup(
        d, AggregationMode::BasicDistinct,
        "bls-single-run/d=" + std::to_string(d));
    const auto statement = make_statement(setup_result.pp);
    const auto proof =
        prove(setup_result.pp, setup_result.aux, statement);
    const auto proof_wire =
        serialize_proof(setup_result.pp, proof);
    const auto crs_wire =
        serialize_public_parameters(setup_result.pp);

    const auto validated =
        deserialize_and_validate_proof(proof_wire, setup_result.pp);
    const auto context = prepare_verifier_context(
        setup_result.pp, setup_result.aux, statement);
    if (!validated || !context)
      throw std::runtime_error("input preparation failed");

    const auto verify_start = Clock::now();
    const bool accepted = verify_online(*context, *validated);
    const auto verify_stop = Clock::now();
    if (!accepted)
      throw std::runtime_error("proof verification failed");

    const auto direct_start = Clock::now();
    const bool direct_accepted =
        direct_bls_verify(setup_result.pp, statement);
    const auto direct_stop = Clock::now();
    if (!direct_accepted)
      throw std::runtime_error("direct verification failed");

    const double verification_ms =
        std::chrono::duration<double, std::milli>(
            verify_stop - verify_start)
            .count();
    const double direct_verification_ms =
        std::chrono::duration<double, std::milli>(
            direct_stop - direct_start)
            .count();

    std::cout << std::fixed << std::setprecision(3)
              << "Verification time: " << verification_ms << " ms\n"
              << "Proof size: " << proof_wire.size() << " bytes\n"
              << "CRS size: " << crs_wire.size() << " bytes\n"
              << "Direct BLS aggregate verification time: "
              << direct_verification_ms << " ms\n";
    return EXIT_SUCCESS;
  } catch (...) {
    return EXIT_FAILURE;
  }
}

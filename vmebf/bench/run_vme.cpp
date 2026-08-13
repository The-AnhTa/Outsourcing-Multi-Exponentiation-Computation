#include "vme_ibf/vme_ibf.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<vme_ibf::G2> public_h(std::size_t n, std::size_t d) {
    std::vector<vme_ibf::G2> result;
    result.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto seed =
            "VME.BF.G2/SINGLE-RUN/H/" + std::to_string(d) + "/" +
            std::to_string(i);
        vme_ibf::G2 point;
        mcl::bn::hashAndMapToG2(point, seed.data(), seed.size());
        result.push_back(point);
    }
    return result;
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("expected one dimension");
        }

        const auto d = static_cast<std::size_t>(std::stoull(argv[1]));
        if (d < 1 || d > 12) {
            throw std::invalid_argument("unsupported dimension");
        }

        vme_ibf::initialize();
        const auto n = std::size_t{1} << d;
        const auto h = public_h(n, d);
        vme_ibf::DeterministicRng rng(
            "VME.BF.G2/SINGLE-RUN/V1/d=" + std::to_string(d));
        const auto setup_result = vme_ibf::setup(d, h, rng);
        const auto phase1 =
            vme_ibf::prove_phase1(setup_result.crs, setup_result.precomp,
                                  setup_result.statement_input);
        const auto phase2 =
            vme_ibf::prove_phase2(setup_result.crs, setup_result.precomp,
                                  phase1);
        const auto proof =
            vme_ibf::assemble_public_proof(phase1, phase2);

        const auto inputs = vme_ibf::validate_verification_inputs(
            setup_result.crs, setup_result.precomp, phase1.statement, proof);
        if (!inputs) {
            throw std::runtime_error("invalid verification inputs");
        }

        const auto verify_start = std::chrono::steady_clock::now();
        const bool accepted = vme_ibf::verify_online(*inputs);
        const auto verify_end = std::chrono::steady_clock::now();
        if (!accepted) {
            throw std::runtime_error("verification failed");
        }

        const double verification_ms =
            std::chrono::duration<double, std::milli>(
                verify_end - verify_start)
                .count();
        const auto proof_bytes =
            vme_ibf::serialize_proof(d, proof).size();
        const auto crs_bytes =
            vme_ibf::serialize_crs(setup_result.crs).size();

        std::cout << std::fixed << std::setprecision(3)
                  << "Verification time: " << verification_ms << " ms\n"
                  << "Proof size: " << proof_bytes << " bytes\n"
                  << "CRS size: " << crs_bytes << " bytes\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

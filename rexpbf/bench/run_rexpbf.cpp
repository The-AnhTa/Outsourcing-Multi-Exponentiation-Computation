#include "rexpbf/prove.hpp"
#include "rexpbf/serialization.hpp"
#include "rexpbf/setup.hpp"
#include "rexpbf/verify.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::vector<std::uint8_t> bytes(const std::string& value) {
    return {value.begin(), value.end()};
}
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("expected one dimension");
        }

        const auto d = static_cast<std::size_t>(std::stoull(argv[1]));
        const auto instance = rexpbf::setup(
            {d, bytes("REXPBF-SINGLE-RUN-V1|CRS"),
             bytes("REXPBF-SINGLE-RUN-V1|H")});
        const auto prove_result =
            rexpbf::prove(instance.crs, instance.precomputation,
                          instance.statement, instance.prover_input);
        const auto validated = rexpbf::validate_verification_inputs(
            instance.crs, instance.precomputation, instance.statement,
            prove_result.proof);
        if (!validated) {
            throw std::runtime_error("verification input validation failed");
        }

        const auto verify_start = std::chrono::steady_clock::now();
        const bool accepted =
            rexpbf::verify_prevalidated(*validated);
        const auto verify_end = std::chrono::steady_clock::now();
        if (!accepted) {
            throw std::runtime_error("verification failed");
        }

        const double verification_ms =
            std::chrono::duration<double, std::milli>(
                verify_end - verify_start)
                .count();
        const auto proof_bytes =
            rexpbf::serialize_proof_wire(prove_result.proof, d,
                                         instance.crs.n)
                .size();
        const auto crs_bytes =
            rexpbf::serialize_crs_wire(instance.crs).size();

        std::cout << std::fixed << std::setprecision(3)
                  << "Verification time: " << verification_ms << " ms\n"
                  << "Proof size: " << proof_bytes << " bytes\n"
                  << "CRS size: " << crs_bytes << " bytes\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "run_rexpbf: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

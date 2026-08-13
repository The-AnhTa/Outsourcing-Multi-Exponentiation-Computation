#include "rexp/rexp.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        if (argc != 2) {
            throw std::invalid_argument("expected one dimension");
        }

        const auto d = static_cast<std::size_t>(std::stoull(argv[1]));
        const auto raw_crs =
            rexp::GenerateRawCRS(d, "rexp-single-run-seed-v1");
        const auto params = rexp::PreparePublicParameters(raw_crs);
        const auto raw_statement = rexp::GenerateRawStatement(params);
        const auto statement =
            rexp::PrepareStatement(params, raw_statement);
        const rexp::RexpProverInput input{raw_statement.H};
        const auto proof = rexp::Prove(params, statement, input);
        const auto validated = rexp::ValidateRexpProof(proof, d);

        const auto verify_start = std::chrono::steady_clock::now();
        const bool accepted =
            rexp::VerifyOptimized(params, statement, validated);
        const auto verify_end = std::chrono::steady_clock::now();
        if (!accepted) {
            throw std::runtime_error("verification failed");
        }

        const double verification_ms =
            std::chrono::duration<double, std::milli>(
                verify_end - verify_start)
                .count();
        const auto proof_bytes = rexp::SerializeRexpProof(proof, d).size();
        const auto crs_bytes = rexp::SerializeRexpCRS(raw_crs).size();

        std::cout << std::fixed << std::setprecision(3)
                  << "Verification time: " << verification_ms << " ms\n"
                  << "Proof size: " << proof_bytes << " bytes\n"
                  << "CRS size: " << crs_bytes << " bytes\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "run_rexp failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "run_rexp failed: unknown error\n";
        return EXIT_FAILURE;
    }
}

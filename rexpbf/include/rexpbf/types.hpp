#pragma once
#include <mcl/bn.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace rexpbf {
using Digest32 = std::array<std::uint8_t, 32>;
using Fr = mcl::bn::Fr;
using G1 = mcl::bn::G1;
using G2 = mcl::bn::G2;
using GT = mcl::bn::GT;

struct CRS {
    std::size_t d{};
    std::size_t n{};
    std::vector<G1> gamma;
    std::vector<G2> lambda;
    Digest32 digest{};
};
struct Precomputation {
    Digest32 crs_digest{};
    std::vector<GT> x;
    std::vector<GT> delta1_right;
    std::vector<GT> delta2_right;
    std::uint64_t pairing_terms{};
};
struct Statement {
    Digest32 crs_digest{};
    std::vector<G1> h;
    GT d1_initial;
    GT e0;
    GT f0;
    GT t_left0;
    GT t_right0;
    Digest32 digest{};
    std::uint64_t pairing_terms{};
};
struct ProverInput { std::vector<G1> h; };
struct SetupResult {
    CRS crs;
    Precomputation precomputation;
    Statement statement;
    ProverInput prover_input;
};
struct SetupConfig {
    std::size_t d{};
    std::vector<std::uint8_t> crs_seed;
    std::vector<std::uint8_t> instance_seed;
};
}

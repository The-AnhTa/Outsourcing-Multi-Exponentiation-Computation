#pragma once

#include <mcl/bn.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <vector>

namespace rexp {

using Fr = mcl::bn::Fr;
using G1 = mcl::bn::G1;
using G2 = mcl::bn::G2;
using GT = mcl::bn::Fp12;
using Digest = std::array<std::uint8_t, 32>;
using NonzeroScalarSampler = std::function<Fr()>;

struct DoryCRS {
    std::size_t d = 0;
    std::size_t n = 0;
    std::vector<G1> Gamma;
    std::vector<G2> Lambda;
    Digest digest{};
};

struct DoryPrecomputation {
    std::vector<GT> X;
    std::vector<GT> Delta1R;
    std::vector<GT> Delta2R;
    Digest crs_digest{};
    std::size_t pairing_product_terms = 0;
};

struct DoryStatement {
    GT D0;
    GT D1;
    GT D2;
};

struct DoryWitness {
    std::vector<G1> Phi;
    std::vector<G2> Theta;
};

struct SetupResult {
    DoryCRS crs;
    DoryPrecomputation precomp;
    DoryStatement statement;
    DoryWitness witness;
};

struct BatchSetupResult {
    DoryCRS crs;
    DoryPrecomputation precomp;
    std::vector<DoryStatement> statements;
    std::vector<DoryWitness> witnesses;
};


void initialize();
DoryCRS GenerateDoryCRS(std::size_t d, std::string_view crs_seed);
DoryPrecomputation PrepareDoryPrecomputation(const DoryCRS& crs);


SetupResult SetupDory(
    std::size_t d,
    std::string_view crs_seed,
    NonzeroScalarSampler sampler = {});



BatchSetupResult SetupDoryBatch(
    std::size_t d,
    std::size_t batch_size,
    std::string_view crs_seed,
    NonzeroScalarSampler sampler = {});

Digest ComputeDoryCRSDigest(const DoryCRS& crs);
bool ValidateCRS(const DoryCRS& crs);
bool ValidatePrecomputation(
    const DoryCRS& crs,
    const DoryPrecomputation& precomp);




std::vector<std::uint8_t> SerializeCRS(const DoryCRS& crs);
std::vector<std::uint8_t> SerializePrecomputation(
    const DoryPrecomputation& precomp);
std::vector<std::uint8_t> SerializeStatement(const DoryStatement& statement);

struct SetupSizes {
    std::size_t core_crs_bytes = 0;
    std::size_t precomputation_bytes = 0;
    std::size_t statement_bytes = 0;
};

SetupSizes MeasureSetupSizes(const SetupResult& setup);
SetupSizes MeasureSetupSizes(const BatchSetupResult& setup);
std::string DigestHex(const Digest& digest);

}

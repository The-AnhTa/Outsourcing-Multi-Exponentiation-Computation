# Rexp BN254 module

This directory implements the random multi-exponentiation protocol and its
embedded Dory argument over MCL's BN254 pairing groups.

## Public API

`include/rexp/rexp.hpp` is the compatibility umbrella. New code may include a
smaller header according to its responsibility:

- `types.hpp` — raw, prepared, and setup data types.
- `proof.hpp` — proof and validated-proof types.
- `setup.hpp` — CRS, public-parameter, and statement preparation.
- `prove.hpp` — Rexp proving.
- `verify.hpp` — reference, prepared, validated, and optimized verification.
- `serialization.hpp` — mathematical payload and canonical wire encoding.
- `metrics.hpp` — validation and verification instrumentation.

The existing Dory headers remain available for direct use of the embedded
argument.

## Object boundaries

- Raw objects are transport or construction inputs and must be prepared before
  repeated use.
- Prepared public parameters and statements validate their inputs and cache
  pairings used by proving and verification.
- `ValidatedRexpProof` is created only after shape, group, and GT subgroup
  checks. Optimized verification accepts this type to avoid repeating those
  checks.
- `PreparedVerifier` owns immutable copies of its Dory CRS and precomputation;
  it does not borrow caller storage.

Boolean verification entry points return `false` for malformed or rejected
proofs. Setup, preparation, serialization, deserialization, and explicit
validation throw `std::exception` subclasses with diagnostic messages.

## Source layout

- `src/internal/crypto.*` contains common canonical encoding, hashing, group
  arithmetic, pairing-product, dimension, and byte-reader utilities.
- `src/internal/dory_transcript.*` and `rexp_transcript.*` contain all
  domain-separated transcript state transitions.
- `src/dory_setup.cpp`, `dory_prove.cpp`, `dory.cpp`, and
  `dory_serialization.cpp` separate Dory setup, proving/batching, verification,
  and serialization.
- `src/rexp_setup.cpp`, `rexp_prove.cpp`, `rexp.cpp`,
  `rexp_validation.cpp`, and `rexp_serialization.cpp` provide the equivalent
  Rexp separation.

The reference verifiers deliberately remain independent from the optimized
symbolic multi-exponentiation paths so tests can use them as correctness
oracles.

## Wire and transcript compatibility

Domain-separation labels, big-endian integer encoding, element framing, and
canonical MCL encodings are compatibility-sensitive. Change them only as a
versioned protocol update. Tests freeze a Fiat–Shamir challenge, a deterministic
Dory proof hash, Rexp payload sizes, and malformed-wire behavior.

## Build and test

From this directory with a configured C++ toolchain:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

On Windows, `run_rexp.ps1` locates Visual Studio Build Tools, configures the
Release build when necessary, and runs one benchmark:

```powershell
.\run_rexp.ps1 --d 10
```

The production statement generator uses MCL's CSPRNG. Tests use an injected,
deterministic nonzero scalar sampler only for reproducible fixtures.

## Profiling and optimization

Define `REXP_ENABLE_PROFILING` to populate per-Dory Rexp metrics. Optimizations
must preserve the transcript and wire fixtures and should be benchmarked in
isolated changes. Pairing products operate directly on checked vector slices,
and rvalue proof validation is available to avoid copying large proofs.

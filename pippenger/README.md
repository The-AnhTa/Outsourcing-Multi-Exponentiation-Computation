# Multi-instance Pippenger and Pinkas verification

This directory contains two C++20 libraries over mcl BN254 G2:

- `pippenger` evaluates one or more multi-scalar multiplications sharing the
  same bases.
- `pinkas` constructs and verifies the Pinkas proof used to validate those
  outputs.

The public headers are `pippenger/pippenger.hpp` and `pinkas/pinkas.hpp`.
The legacy flat headers remain as forwarding headers for source compatibility.

## Safety and compatibility

- Window widths are restricted to `1..maximum_window_width()` (currently 16)
  before bucket allocation. The benchmark uses width 8.
- Bounded MSM rejects scalars that exceed the declared unsigned exponent bit
  length.
- Pinkas public parameters, dimensions, curve points, proof shape, and canonical
  encodings are validated at their API boundaries.
- `ValidatedInputs` owns its validated state; caller mutation or destruction
  cannot invalidate online verification.
- Valid proofs retain protocol identifier `PinkasFS128BN-v2` and the existing
  serialized wire format. Tests lock deterministic transcript and proof hashes.
- `deserialize_proof_detailed` reports structured decode failures. The existing
  `deserialize_proof` boolean API remains available.

## Build and test

From a Visual Studio x64 developer shell:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

mcl v3.00 is fetched with a pinned SHA-256 hash. Generated `build*` directories
are ignored by Git.

## Benchmarks

```powershell
.\run_pippenger.ps1 --d 8 --k 2
.\run_pinkas.ps1 --d 8 --k 2
.\run_pinkas_all.ps1
```

`d` defines `n = 2^d` shared bases and `k` is the number of instances. The
executables are benchmarks, while correctness and rejection behavior are
covered by `pippenger_tests`.

## Source organization

- `include/`: public APIs and compatibility forwarding headers.
- `src/internal/`: private constants, transcript, arithmetic, and validation
  declarations.
- `src/`: initialization, MSM core, protocol, prover, verifier, validation,
  transcript, and serialization implementations.
- `bench/`: benchmark executables and shared benchmark helpers.
- `tests/`: deterministic compatibility, correctness, boundary, tamper, and
  malformed-encoding tests.

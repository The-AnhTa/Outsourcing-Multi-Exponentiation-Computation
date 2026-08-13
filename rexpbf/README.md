# REXP-BF

`rexpbf` is the BN254 implementation of the batched recursive
multi-exponentiation outsourcing protocol in this repository. The refactored
library keeps protocol construction, transcript replay, wire encoding,
validation, reference verification, and optimized verification in separate
modules.

## Build and test

The build requires CMake 3.20 or newer, a C++20 compiler, and network access on
the first configure so CMake can fetch the pinned mcl v3.00 dependency.

On Windows with Visual Studio Build Tools and Ninja installed:

```powershell
./run_rexpbf.ps1 --d 3
```

For an explicit build and test run from a Visual Studio developer shell:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

To use an existing mcl checkout without downloading it, configure with
`-DFETCHCONTENT_SOURCE_DIR_MCL=<path-to-mcl>`.

## Library flow

Include `rexpbf/rexpbf.hpp` for the complete public API, or include individual
headers for narrower dependencies.

1. Call `setup` with a dimension and domain-specific seeds.
2. Call `prove` with the returned CRS, precomputation, statement, and prover
   input.
3. Call `verify` for the safe validation-plus-verification path.
4. Use `verify_reference` only as an independent correctness oracle.

`verify_prevalidated` accepts an owning `ValidatedVerificationInputs` token.
Construct the token with `validate_verification_inputs`; this makes the trust
boundary explicit and prevents later mutation of the validated objects.

## Compatibility and tests

The test suite checks deterministic setup, reference/optimized agreement,
tamper rejection, canonical wire round trips, malformed input rejection,
pairing and multi-exponentiation helpers, profiling invariants, and size
baselines for dimensions 1 through 4. Golden SHA-256 values pin the proof wire
encoding and Fiat-Shamir transcript, so refactors cannot silently change the
protocol.

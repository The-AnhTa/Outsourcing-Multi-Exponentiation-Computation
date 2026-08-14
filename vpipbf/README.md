# VPIP.BF

This directory contains the BN254 implementation of the VPIP.BF protocol.
It provides deterministic and secure setup, proving, canonical wire encoding,
and reference and optimized verification APIs.

## Build and test

From an MSVC x64 developer command prompt:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Run the benchmark for dimension `d` with:

```powershell
.\build\run_vpip.exe 4
```

or use the wrapper:

```powershell
.\run_vpip.ps1 --d 4
```

The public statement contains `n = 2^d` G1 points and their pairing-product
commitment. Verification first validates the CRS, audited precomputation,
statement, and proof before replaying the Fiat-Shamir transcript.

The wire format is version 2. Refactoring must preserve transcript domains,
challenge ordering, canonical encodings, and established proof and CRS sizes.

## Verification boundary

The ordinary verification APIs validate and audit every supplied object.
`PrevalidateVerificationInputs` returns an owning token, so subsequent source
mutation or destruction cannot invalidate online verification. Transcript and
unchecked equation helpers live under `src/internal` and are not public API.

The test suite locks deterministic protocol hashes and wire sizes, compares
the direct reference and optimized symbolic verification paths, and exercises
malformed phase boundaries, serialization failures, and proof tampering.

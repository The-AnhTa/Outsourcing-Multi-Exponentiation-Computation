# BLS aggregate proof

This directory implements the non-zero-knowledge BLS aggregate outsourcing
protocol over BN254. It supports distinct-message Basic aggregation and
public-key-augmented aggregation.

The public C++20 API is in `include/blsagg`, implementation code is in `src`,
deterministic regression tests are in `tests`, and `bench/run_bls.cpp` is the
single-run benchmark. Setup/precomputation, transcript handling, canonical
serialization, centralized validation, and shared cryptographic helpers are
separate modules; protocol proving and verification remain in `protocol.cpp`.

## Build and test

Use CMake 3.20 or newer. MCL 3.00 is pinned by URL and SHA-256 digest.

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

Set `-DBLS_BUILD_BENCHMARKS=OFF` for a library-and-tests-only build, or
`-DBUILD_TESTING=OFF` when embedding only the library.

## Run

The PowerShell runner accepts `d` from 1 through 12 and either `basic` (the
default) or `augmented` mode. The aggregate contains `k = 2^d` messages and
public keys.

```powershell
.\run_bls.ps1 --d 2
.\run_bls.ps1 --d 2 --mode augmented
```

## API safety

`verify_safe` validates parameters, precomputation, statement, and proof in a
single call. Repeated verification should prepare a `ValidatedVerifierContext`
and deserialize a `ValidatedProof`, then use a named online verification
strategy. Versioned wire decoders reject non-canonical group elements,
dimension mismatches, malformed lengths, and trailing data. Decoder failures
clear their destination objects. A `ValidatedProof` is immutable and bound to
both its canonical wire representation and parameter digest; online verifier
entry points enforce that binding.

Named verifier strategies expose sequential MSM, parallel MSM, symbolic GT,
and split-G2 MSM execution. All strategies use one ordered transcript replay
implementation, and the regression suite compares their challenges and
terminal targets. Transcript behavior is additionally protected by a golden
challenge fixture.

The small `proof.hpp`, `prove.hpp`, and `verify.hpp` headers are compatibility
forwarders to `protocol.hpp`; new code may include `protocol.hpp` directly.

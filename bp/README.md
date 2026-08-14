# BP protocols

This directory implements the base inner-product proof (`BP`) and its
helper-prover (`HP`) and helper-verifier (`HV`) outsourcing variants. The
public C++20 API is under `include/bp`; implementation-only protocol machinery
is under `src`; deterministic regression tests are under `tests`; and the
`bench` programs report timings and wire sizes.

## Build and test

The project uses CMake 3.20 or newer and fetches the pinned MCL 3.00 source.
On Windows, run these commands from a Visual Studio x64 developer shell:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build
ctest --test-dir build --output-on-failure
```

Set `-DBP_BUILD_BENCHMARKS=OFF` for a library-and-tests-only build. Set
`-DBUILD_TESTING=OFF` when embedding only the library.

## Benchmarks

The PowerShell entry points validate `d` in the range 1 through 20 and build a
Release executable on first use. The vector length is `n = 2^d`.

```powershell
.\run_bp.ps1 --d 4
.\run_helper_prover.ps1 --d 4
.\run_helper_verifier.ps1 --d 4
```

## Compatibility and validation

Proof and statement encodings are versioned and decoded canonically. Public
parameters, proof dimensions, transcript bindings, and precomputations are
validated before cryptographic work. The tests cover deterministic positive
flows, reference/optimized agreement, serialization round trips, cache reuse,
and rejection of malformed or tampered inputs.

`HvPreparedStatementCache` intentionally exposes only `ready()` and `clear()`;
its prepared scalars and relation point are private protocol state. Reuse is
bound to both the CRS digest and the complete statement context.

## Implementation map

- `src/internal/protocol_utils.*`: checked framing, canonical decoding,
  arithmetic, dimension, and wire-size primitives.
- `src/internal/protocol_validation.*`: shared BP/HP/HV structural validation.
- `src/internal/bp_transcript.*`: BP Fiat-Shamir transcript and deterministic
  generator derivation.
- `src/helper_prover.cpp` and `src/helper_verifier.cpp`: outsourcing workflow
  orchestration and versioned wire codecs.
- `src/hp_vme_internal.cpp`: explicit REXP/Dory proving phases and separate
  reference/optimized VME terminal verifiers.

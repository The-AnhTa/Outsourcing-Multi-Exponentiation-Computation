# VME.bf

This directory contains the BLS12-381 implementation of the VME.bf protocol. It
provides deterministic setup and proving, reference and optimized verification,
canonical wire serialization, and regression tests.

## Build and test

Requirements are CMake 3.20 or newer, a C++20 compiler, and internet access for
CMake's pinned `mcl` dependency on the first configuration.

On Windows with Visual Studio C++ Build Tools and Ninja:

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

The helper script configures and builds the benchmark automatically:

```powershell
.\run_vme.ps1 --d 3
```

Supported benchmark dimensions are 1 through 12. Runtime and memory grow with
`2^d`, so begin with a small dimension.

## Library entry points

Include `vme_ibf/vme_ibf.hpp` for the complete public API. A normal flow is:

1. Call `initialize()` once.
2. Create the public G2 vector and call `setup()`.
3. Call `prove_phase1()`, `prove_phase2()`, and `assemble_public_proof()`.
4. Establish the verification trust boundary with
   `validate_verification_inputs()`.
5. Pass the validated object to `verify_online()`.

`verify_reference()` provides an independent arithmetic implementation for
cross-checking. The deferred and combined APIs expose diagnostic traces intended
for testing and performance analysis.

The serialization API rejects malformed, non-canonical, wrong-curve, and
digest-inconsistent objects. Treat deserialized data as untrusted until its
corresponding validation function succeeds.

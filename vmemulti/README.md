# VMEMulti

This directory implements the multi-instance verifiable multi-exponentiation
protocol over BN254. It binds multiple public exponent vectors and outputs,
derives one Fiat–Shamir aggregation challenge, reduces them to one VME instance,
and proves that aggregate with a proof whose size is independent of the number
of instances.

## Build and test

Requirements are CMake 3.20 or newer, a C++20 compiler, and internet access for
the pinned `mcl` dependency during the first configuration.

```powershell
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
ctest --test-dir build --output-on-failure
```

The Windows helper configures and builds the benchmark automatically:

```powershell
.\run_vmemulti.ps1 --d 3 --m 4
```

Dimensions 1–12 and instance counts 1–1024 are supported. Runtime and memory
grow with both `2^d` and the instance count, so begin with small parameters.

## Public flow

Include `vme_ibf/vmemulti.hpp` and:

1. Call `initialize()` once.
2. Build the public G2 vector and call `setup()`.
3. Construct a rectangular `VmeMultiStatement` whose outputs correspond to its
   exponent vectors.
4. Call `prove_vmemulti()`.
5. Call `validate_vmemulti_inputs()` to establish an owning trust boundary.
6. Pass the validated value to `verify_vmemulti_online()`.

`verify_vmemulti_diagnostic()` is the independent reference verifier.
`verify_vmemulti_combined_diagnostic()` and the online trace API expose detailed
evaluation and aggregation timings for testing and performance work.

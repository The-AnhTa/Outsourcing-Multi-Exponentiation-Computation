# MultiExp

This repository contains research prototype implementations of outsourcing
protocols and algorithms for group multi-exponentiation computation.

The code is intended for research and benchmarking. It has not been audited
for production use.

## Requirements

- 64-bit Windows
- PowerShell 5.1 or later
- CMake 3.20 or later
- Ninja
- Visual Studio C++ Build Tools with the x64 compiler
- Internet access for the first build

The run scripts locate Visual Studio through `vswhere.exe`. During the first
build, CMake downloads the pinned MCL v3.00 dependency and verifies its
SHA-256 hash. Generated files are stored in a project-local `build` directory,
which is ignored by Git.

## Projects

| Folder | Contents |
| --- | --- |
| `bls` | Batched-verification outsourcing for aggregate BLS signatures (`bls.agg.bf`). |
| `bp` | Bulletproof inner-product protocol, helper-prover protocol, and helper-verifier protocol. |
| `pippenger` | Multi-instance Pippenger multi-exponentiation and the Fiat-Shamir Pinkas protocol with 128-bit challenges. |
| `rexp` | Recursive exponentiation protocol. |
| `rexpbf` | Batch-and-fold recursive exponentiation protocol (`rexp.bf`). |
| `vmebf` | Batch-and-fold verifiable multi-exponentiation protocol (`vme.bf`). |
| `vmemulti` | Multi-instance verifiable multi-exponentiation protocol. |
| `vpipbf` | Batch-and-fold verifiable pairing inner-product protocol. |

## Running

Each PowerShell script configures and builds its required executable when
needed, then runs the protocol once. There are no benchmark warm-up runs or
repeated iterations.

The parameter `d` defines the vector length:

```text
n = 2^d
```

Run the commands from the corresponding project folder.

### Aggregate BLS

```powershell
cd bls
.\run_bls.ps1 --d 10
```

The script reports protocol verification time, proof size, CRS size, and
direct aggregate-BLS verification time.

### Bulletproof and helper protocols

```powershell
cd bp
.\run_bp.ps1 --d 10
.\run_helper_prover.ps1 --d 10
.\run_helper_verifier.ps1 --d 10
```

`run_bp.ps1` reports verification time, prover time, proof size, and CRS size.
The helper-prover and helper-verifier scripts report their relevant online
time, proof size, and CRS size.

### Pippenger

```powershell
cd pippenger
.\run_pippenger.ps1 --d 10 --k 1
```

Here, `k` is the number of instances. The script reports running time.

### Pinkas

```powershell
cd pippenger
.\run_pinkas.ps1 --d 10 --k 1
```

The script reports verification time and proof size.

To run all combinations of `d = 8, 9, 10, 11, 12` and
`k = 1, 2, 4, 8`, with a 10-second pause after every run:

```powershell
.\run_pinkas_all.ps1
```

### Recursive exponentiation

```powershell
cd rexp
.\run_rexp.ps1 --d 10
```

The script reports verification time, proof size, and CRS size.

### Batch-and-fold recursive exponentiation

```powershell
cd rexpbf
.\run_rexpbf.ps1 --d 10
```

The script reports verification time, proof size, and CRS size.

### Batch-and-fold verifiable multi-exponentiation

```powershell
cd vmebf
.\run_vme.ps1 --d 10
```

The script reports verification time, proof size, and CRS size.

### Multi-instance verifiable multi-exponentiation

```powershell
cd vmemulti
.\run_vmemulti.ps1 --d 10 --m 4
```

Here, `m` is the number of instances. The script reports verification time,
proof size, and CRS size.

### Verifiable pairing inner product

```powershell
cd vpipbf
.\run_vpip.ps1 --d 10
```

The script reports verification time, proof size, CRS size, and direct pairing
product time.

## Clean rebuild

Delete the `build` directory inside a project, then run its PowerShell script
again. The script will configure the project and rebuild its executable from
scratch.

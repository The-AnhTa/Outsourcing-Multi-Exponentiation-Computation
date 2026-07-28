### Outsourcing Protocols for Group Multi-Exponentiation Computation

This repository contains research prototype implementations of outsourcing
protocols and algorithms for group multi-exponentiation computation.
The code is intended for research and benchmarking. 

## Requirements

- CMake >= 3.20
- Ninja
- Visual Studio C++ Build Tools with the x64 compiler

The run scripts locate Visual Studio through `vswhere.exe`. During the first
build, CMake downloads the pinned MCL v3.00 dependency.

## Projects

| Folder | Contents |
| --- | --- |
| `bls` | Batched verification outsourcing for aggregate BLS signatures. |
| `bp` | Bulletproof inner-product protocol, helper-prover protocol, and helper-verifier protocol. |
| `pippenger` | Multi-instance Pippenger multi-exponentiation and the Pinkas protocol with 128-bit challenges. |
| `rexp` | Random multi-exponentiation protocol. |
| `rexpbf` | Batch-and-fold random multi-exponentiation protocol. |
| `vmebf` | Batch-and-fold verifiable multi-exponentiation protocol. |
| `vmemulti` | Multi-instance verifiable multi-exponentiation protocol. |
| `vpipbf` | Batch-and-fold verifiable pairing inner-product protocol. |

## Running

Each PowerShell script .ps1 configures and builds its required executable when
needed, then runs the protocol once. 

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
The helper-prover and helper-verifier scripts report the outsourcer verification time, proof size, and CRS size.

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


```powershell
.\run_pinkas_all.ps1
```

### Random multi-exponentiation

```powershell
cd rexp
.\run_rexp.ps1 --d 10
```

The script reports verification time, proof size, and CRS size.

### Batch-and-fold random multi-exponentiation

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

Here, `m` is the number of instances. The script reports verification time, proof size, and CRS size.

### Verifiable pairing inner product

```powershell
cd vpipbf
.\run_vpip.ps1 --d 10
```

The script reports verification time, proof size, CRS size, and direct pairing
product time.

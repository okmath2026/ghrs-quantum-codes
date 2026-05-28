# GHRS Quantum Code Certificates

This repository contains exact verification code for the paper
*Hamming Quantum MDS Codes from Generalized Hyperderivative Reed--Solomon
Codes*.

The main library file is:

- `ghrs_quantum_certificates.magma`

It defines routines for:

- finite-field GHRS generator-matrix construction;
- relative Hamming distance computation by rank;
- rational `s=2` Hermite rank certificates;
- exceptional prime set computation;
- prime-field square-class and density checks.

Two short Magma driver files are provided:

- `s2_bad_primes.magma`, for the rational `s=2` certificate computation;
- `ghrs_search.magma`, for the finite-field search.

For the five-point family in the paper, run in Magma:

```magma
load "s2_bad_primes.magma";
```

This returns the `210` rank conditions, `420` determinant/minor witnesses, the
finite exceptional prime set, the five values `-2S_j`, the square-class rank,
and the first prime-field hits.

The Magma files are the verification artifacts intended for submission.

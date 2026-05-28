# GHRS Quantum Code Certificates

This repository contains exact verification code for the paper
*Hamming Quantum MDS Codes from Generalized Hyperderivative Reed--Solomon
Codes*.

The main reproducibility file is:

- `ghrs_quantum_certificates.magma`

It contains:

- finite-field GHRS generator-matrix construction;
- relative Hamming distance computation by rank;
- rational `s=2` Hermite rank certificates;
- exceptional prime set computation;
- prime-field square-class and density checks.

For the five-point family in the paper, run in Magma:

```magma
load "ghrs_quantum_certificates.magma";
PrintS2Certificate([3,6,12,17,24] : PrimeHits := 10);
```

This returns the `210` rank conditions, `420` determinant/minor witnesses, the
finite exceptional prime set, the five values `-2S_j`, the square-class rank,
and the first prime-field hits.

The Python files are optional development companions; the Magma file is the
single-file verification artifact intended for submission.

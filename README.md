# Dense Matrix QR Eigenvalue Approximation

This repository contains a real-arithmetic shifted QR algorithm for approximating all eigenvalues of a real square matrix. It first reduces the matrix to Hessenberg form, then performs shifted QR steps. Small subdiagonal entries split the matrix into independent blocks; 1-by-1 and 2-by-2 blocks are read directly, with 2-by-2 blocks allowing complex-conjugate eigenvalue pairs.

Eigenvalues are printed as `(real, imaginary)`. The implementation uses only the C++ standard library and the matrix utilities included here.

## Build and run

With MinGW g++ installed at `C:\MinGW\bin\g++.exe`, open the repository in VS Code and select **Terminal → Run Task**. The tasks build and run the 25-by-25 Sampling matrix, the 24-by-24 Hanowa matrix, or the small verification example.

To build and run from a terminal at the repository root:

```text
g++ -std=c++17 -I include src/eigentest.cpp src/QR_based_Approximation.cpp src/factorizations.cpp src/NLA_FUNCTIONS.cpp src/MATRIX.cpp -o eigenvalue.exe
eigenvalue.exe data/eigentest_sampling_25.txt
```

Pass a different whitespace-separated square matrix file as the program argument, for example:

```text
eigenvalue.exe data/eigentest_hanowa_24.txt
```

Each input file has one matrix row per line, with numeric entries separated by spaces. No header or dimensions line is used.

## Included test matrices

- `data/eigentest.txt` — the 7-by-7 example used with the verification program.
- `data/eigentest_sampling_25.txt` — the 25-by-25 dense MATLAB Gallery Sampling matrix. Its exact eigenvalues are 0 through 24; these eigenvalues are ill-conditioned, making this a demanding numerical test.
- `data/eigentest_hanowa_24.txt` — a dense orthogonal similarity of the 24-by-24 Hanowa matrix with `alpha = -1`. Its exact eigenvalues are `-1 ± k i` for `k = 1, ..., 12`.

The test matrix definitions and spectra are documented in the [MathWorks Gallery reference](https://www.mathworks.com/help/matlab/ref/gallery.html), under `sampling` and `hanowa`.

## Verification example

`src/eigenvalue_varification.cpp` checks the 7-by-7 example by comparing `trace(A^k)` with the sum of the listed eigenvalues raised to powers `k = 1, ..., 7`. The example eigenvalues are rounded, so its tolerance accounts for display rounding.

## Numerical limitation

This is a compact educational implementation using a single shift. Some matrices may not converge within the iteration limit, and the ill-conditioned Sampling matrix can show noticeable eigenvalue error. A production dense nonsymmetric solver typically uses an implicit multishift QR method.

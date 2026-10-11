#include "matrix.h"
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using ComplexPair = std::array<double, 2>; // {real part, imaginary part}

ComplexPair multiply(const ComplexPair &a, const ComplexPair &b)
{
    return {a[0] * b[0] - a[1] * b[1],
            a[0] * b[1] + a[1] * b[0]};
}

double trace(const Matrix &A)
{
    double result = 0.0;
    for (unsigned i = 0; i < A.getRows(); ++i)
        result += A(i, i);
    return result;
}

int main()
{
    const Matrix A("data/eigentest.txt");
    const std::vector<ComplexPair> eigenvalues = {
        {-8.10306, 0.0},
        {-10.8857, 0.0},
        {1.03435, 0.0},
        {12.8752, 0.0},
        {9.4564, 2.69583},
        {9.4564, -2.69583},
        {8.16633, 0.0}
    };

    if (A.getRows() != A.getCols() || eigenvalues.size() != A.getRows())
    {
        std::cerr << "The matrix must be square and have one eigenvalue per row.\n";
        return 1;
    }

    const unsigned n = A.getRows();
    Matrix matrixPower(n, n, 0.0);
    for (unsigned i = 0; i < n; ++i)
        matrixPower(i, i) = 1.0;

    std::vector<ComplexPair> eigenvaluePowers(n, {1.0, 0.0});
    bool allPassed = true;

    std::cout << std::setprecision(10)
              << "Compare trace(A^k) with the sum of lambda^k.\n"
              << "The tolerance accounts for eigenvalues rounded in the output.\n\n";

    for (unsigned k = 1; k <= n; ++k)
    {
        matrixPower = matrixPower * A;
        ComplexPair eigenvalueSum = {0.0, 0.0};
        double sumOfMagnitudes = 0.0;

        for (unsigned i = 0; i < n; ++i)
        {
            eigenvaluePowers[i] = multiply(eigenvaluePowers[i], eigenvalues[i]);
            eigenvalueSum[0] += eigenvaluePowers[i][0];
            eigenvalueSum[1] += eigenvaluePowers[i][1];
            sumOfMagnitudes += std::hypot(eigenvaluePowers[i][0], eigenvaluePowers[i][1]);
        }

        const double matrixTrace = trace(matrixPower);
        const double error = std::hypot(matrixTrace - eigenvalueSum[0], eigenvalueSum[1]);
        const double tolerance = 5e-6 * sumOfMagnitudes + 1e-8;
        const bool passed = error <= tolerance;
        allPassed = allPassed && passed;

        std::cout << "k=" << k
                  << "  trace(A^k)=" << matrixTrace
                  << "  sum(lambda^k)=(" << eigenvalueSum[0] << ", " << eigenvalueSum[1] << ")"
                  << "  error=" << error
                  << "  " << (passed ? "PASS" : "FAIL") << '\n';
    }

    std::cout << "\nOverall: " << (allPassed ? "PASS" : "FAIL") << '\n';
    return allPassed ? 0 : 1;
}

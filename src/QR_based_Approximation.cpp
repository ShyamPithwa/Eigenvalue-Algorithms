#include "Eigen_Approximation.h"
#include "NLA_FUNCTIONS.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace
{
void applySimilarity(Matrix &A_current, const Matrix &Q_k, int first, int last)
{
    const int n = static_cast<int>(A_current.getRows());
    const int size = last - first + 1;
    Matrix updated(A_current);

    // H <- Q^T H Q, applying Q only to the selected diagonal block.
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < n; ++j)
        {
            double value = 0.0;
            for (int k = 0; k < size; ++k)
                value += Q_k(k, i) * A_current(first + k, j);
            updated(first + i, j) = value;
        }

    A_current = updated;
    updated = A_current;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < size; ++j)
        {
            double value = 0.0;
            for (int k = 0; k < size; ++k)
                value += A_current(i, first + k) * Q_k(k, j);
            updated(i, first + j) = value;
        }
    A_current = updated;
}

double deflationThreshold(const Matrix &A_current, int i, double tolerance)
{
    const double scale = std::abs(A_current(i, i)) + std::abs(A_current(i + 1, i + 1));
    return tolerance * scale + std::numeric_limits<double>::min();
}

void appendTwoByTwoEigenvalues(const Matrix &A_current, int first, EigenvalueMatrix &values)
{
    const double a = A_current(first, first);
    const double b = A_current(first, first + 1);
    const double c = A_current(first + 1, first);
    const double d = A_current(first + 1, first + 1);
    const double halfTrace = (a + d) / 2.0;
    const double discriminant = (a - d) * (a - d) / 4.0 + b * c;

    if (discriminant >= 0.0)
    {
        const double root = std::sqrt(discriminant);
        values.push_back({halfTrace + root, 0.0});
        values.push_back({halfTrace - root, 0.0});
    }
    else
    {
        const double imaginary = std::sqrt(-discriminant);
        values.push_back({halfTrace, imaginary});
        values.push_back({halfTrace, -imaginary});
    }
}
}

EigenvalueMatrix QR_based_Approximation(const Matrix &A, double tolerance, unsigned maxIterations)
{
    if (A.getRows() != A.getCols())
        throw std::invalid_argument("QR eigenvalue approximation requires a square matrix.");
    if (tolerance <= 0.0 || !std::isfinite(tolerance))
        throw std::invalid_argument("Tolerance must be finite and positive.");

    Householder_Hessenberg reduction(A);
    // A_current is A^(0) in the textbook notation: Q^(0)^T A Q^(0).
    Matrix A_current = reduction.H;
    const int n = static_cast<int>(A_current.getRows());
    unsigned iteration = 0;
    EigenvalueMatrix eigenvalues;

    while (true)
    {
        // Deflate every negligible subdiagonal, splitting H into independent blocks.
        for (int i = 0; i + 1 < n; ++i)
        {
            if (std::abs(A_current(i + 1, i)) <= deflationThreshold(A_current, i, tolerance))
                A_current(i + 1, i) = 0.0;
        }

        int firstUnresolved = -1;
        for (int first = 0; first < n;)
        {
            int last = first;
            while (last + 1 < n && A_current(last + 1, last) != 0.0)
                ++last;
            const int size = last - first + 1;

            if (size > 2)
            {
                firstUnresolved = first;
                break;
            }
            first = last + 1;
        }

        if (firstUnresolved < 0)
        {
            for (int first = 0; first < n;)
            {
                int last = first;
                while (last + 1 < n && A_current(last + 1, last) != 0.0)
                    ++last;
                if (last == first)
                    eigenvalues.push_back({A_current(first, first), 0.0});
                else
                    appendTwoByTwoEigenvalues(A_current, first, eigenvalues);
                first = last + 1;
            }
            return eigenvalues;
        }
        if (iteration >= maxIterations)
            throw std::runtime_error("QR iteration did not converge within maxIterations.");

        int last = firstUnresolved;
        while (last + 1 < n && A_current(last + 1, last) != 0.0)
            ++last;
        const int size = last - firstUnresolved + 1;
        Matrix A_block = Sub_Matrix(A_current, firstUnresolved, last, firstUnresolved, last);
        const double mu_k = A_current(last, last);
        Matrix A_shifted = A_block - identity(size) * mu_k;
        Householder_QR qr(A_shifted);
        const Matrix &Q_k = qr.Q;
        const Matrix &R_k = qr.R;
        Matrix A_next = R_k * Q_k + identity(size) * mu_k;
        applySimilarity(A_current, Q_k, firstUnresolved, last);
        InsertSubMatrix(A_current, A_next, firstUnresolved, firstUnresolved);
        ++iteration;
    }
}

void printEigenvalues(const EigenvalueMatrix &eigenvalues)
{
    for (const auto &value : eigenvalues)
        std::cout << '(' << value[0] << ", " << value[1] << ")\n";
}

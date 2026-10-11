#pragma once

#include "matrix.h"
#include "factorizations.h"
#include <array>
#include <vector>

// Each row stores {real part, imaginary part}.
using EigenvalueMatrix = std::vector<std::array<double, 2>>;

EigenvalueMatrix QR_based_Approximation(const Matrix &A,
                                        double tolerance = 1e-10,
                                        unsigned maxIterations = 1000);

void printEigenvalues(const EigenvalueMatrix &eigenvalues);

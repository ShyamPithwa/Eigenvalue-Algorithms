#include "Eigen_Approximation.h"
#include <string>

int main(int argc, char *argv[])
{
    const std::string inputFile = argc > 1
        ? argv[1]
        : "data/eigentest_sampling_25.txt";
    Matrix A(inputFile);
    double tolerance = 1e-6;
    unsigned maxIterations = 1000;

    try
    {
        EigenvalueMatrix eigenvalues = QR_based_Approximation(A, tolerance, maxIterations);
        std::cout << "Eigenvalues:\n";
        printEigenvalues(eigenvalues);
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}

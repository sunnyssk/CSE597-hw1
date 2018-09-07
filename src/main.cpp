#include <iostream>
#include <cmath>
#include <exception>

#include "matrix.h"

int main (int argc, char** argv) {
    double buffer1[9] = {1, 1, 1, 1, 2, 2, 1, 2, 3};
    double buffer2[9] = {1, 1, 1, 0, 1, 1, 0, 0, 1};
    double buffer3[9] = {1, 0, 0, 1, 1, 0, 1, 1, 1};
    MatD mat(3, 3, buffer1);
    MatD L(3, 3, buffer2), U(3, 3, buffer3);

    std::cout << "mat =\n" << mat << std::endl;
    mat.LUDecomposition(L, U);
    std::cout << "L =\n" << L << std::endl;
    std::cout << "U =\n" << U << std::endl;
    std::cout << "L * U =\n" << L * U << std::endl;
    return 0;
}
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <exception>

#include "matrix.h"

int main (int argc, char** argv) {
    srand(time(NULL));
    double buffer1[9] = {1.0, 4.0, 4.0, 2.0, 4.0, 6.0, 2.0, 2.0, 4.0};
    int size = 6;
    MatD A(size, size);
    A.Rand();
    MatD P(size, size), L(size, size), U(size, size);
    MatD invL(size, size), invU(size, size), invA(size, size);

    std::cout << "A = \n" << A << std::endl;
    A.PLUDecomposition(P, L, U);    
    std::cout << "P = \n" << P << std::endl;
    std::cout << "L = \n" << L << std::endl;
    std::cout << "U = \n" << U << std::endl;
    std::cout << "P * A - L * U = \n" << (P * A - L * U) << std::endl;

    L.LInverse(invL);
    U.UInverse(invU);
    invA = invU * invL * P;
    std::cout << "A * inv(A) = \n" << (A * invA) << std::endl;
    return 0;
}

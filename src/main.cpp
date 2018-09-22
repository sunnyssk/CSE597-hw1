#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <exception>
#include <stdio.h>

#include "matrix.h"
#include "debye.h"

int main (int argc, char** argv) {
    srand(time(NULL));
    int size = 7;
    int nx = 30, ny = 20, nz = 10;
    double debye_length = 1E-7;
    Field3D rhs(nx, ny, nz, 1E-9), potential1(nx, ny, nz, 1E-9), potential2(nx, ny, nz, 1E-9);
    double ee = 1.60217662E-19, e0 = 8.854187817E-12;
    for (int i = 0; i < 30; i++) rhs(rand() % nx, rand() % ny, rand() % nz) = -1 * ee / e0 * 1E27;       // -rho / epsilon_0

    // Direct Solver
    DebyeLUSolver dsolve;
    dsolve.GenerateSolverMatrix(rhs, debye_length);
    std::cout << "System matrix generated." << std::endl;
    dsolve.SolverMatrixDecompose();
    std::cout << "Matrix decomposition completed." << std::endl;
    dsolve.RhsInput(rhs);
    dsolve.LUSolve(potential1);
    std::cout << "Direct solver completed." << std::endl;
    FILE *fout = fopen("output/potential-direct.txt", "w");
    fprintf(fout, "%d %d\n", nx, ny);
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) fprintf(fout, "%lf\n", potential1(i, j, nz / 2));
    fclose(fout);

    // Iterative Solver
    DebyeJacobiSolve(rhs, potential2, debye_length, 1E-10);
    fout = fopen("output/potential-iterative.txt", "w");
    fprintf(fout, "%d %d\n", nx, ny);
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++) fprintf(fout, "%lf\n", potential2(i, j, nz / 2));
    fclose(fout);
    return 0;
}

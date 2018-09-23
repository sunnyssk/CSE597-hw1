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
    if(argc != 2){
        std::cout << "Usage: main [PATH-OF-INPUT]\n"
                  << "The input file should contain one line with three numbers, (nx, ny, nz) of the system.\n"
                  << "Please try again with the correct grammar." << std::endl;
        return 1;
    }
    int nx, ny, nz;
    FILE *fin = fopen(argv[1], "r");
    fscanf(fin, "%d %d %d", &nx, &ny, &nz);
    int n = nx * ny * nz;
    fclose(fin);
    double debye_length = 1E-7;
    double * error_array = new double[MAX_ITER_NUM];
    Field3D rhs(nx, ny, nz, 1E-9), potential1(nx, ny, nz, 1E-9), potential2(nx, ny, nz, 1E-9), potential3(nx, ny, nz, 1E-9);
    double ee = 1.60217662E-19, e0 = 8.854187817E-12;
    for (int i = 0; i < 30; i++) rhs(rand() % (nx - 2) + 1, rand() % (ny - 2) + 1, rand() % (nz - 2) + 1) = -1 * ee / e0 * 1E27;       // -rho / epsilon_0
    std::clock_t start_time;

    // Direct Solver
    std::cout << "===============\n"
              << " Direct Solver \n"
              << "===============\n" << std::endl;
    DebyeLUSolver dsolve;
    start_time = std::clock();
    dsolve.GenerateSolverMatrix(rhs, debye_length);
    std::cout << "System matrix generated. Time: " << double(std::clock() - start_time) / (CLOCKS_PER_SEC / 1000) << " ms." << std::endl;
    start_time = std::clock();
    dsolve.SolverMatrixDecompose();
    std::cout << "Matrix decomposition completed. Time: " << double(std::clock() - start_time) / (CLOCKS_PER_SEC / 1000) << " ms." << std::endl;
    dsolve.RhsInput(rhs);
    start_time = std::clock();
    dsolve.LUSolve(potential1);
    std::cout << "Direct solver completed. Time: " << double(std::clock() - start_time) / (CLOCKS_PER_SEC / 1000) << " ms." << std::endl;
    FILE *fout = fopen("output/potential-direct.txt", "w");
    potential1.WriteField(fout);
    fclose(fout);
    MatD invL(n, n), invU(n, n);
    start_time = std::clock();
    dsolve.Lmat().LInverse(invL);
    dsolve.Umat().UInverse(invU);
    std::cout << "Inverse matrices solved. Time: " << double(std::clock() - start_time) / (CLOCKS_PER_SEC / 1000) << " ms." << std::endl;
    start_time = std::clock();
    MatD potential_vector(n, 1);
    potential_vector.Copy(invU * (invL * (dsolve.Pmat() * dsolve.bvec())));
    std::cout << "Applied direct inverse multiplication. Time: " << double(std::clock() - start_time) / (CLOCKS_PER_SEC / 1000) << " ms." << std::endl;

    // Iterative Solver: Zero Initial Value
    std::cout << "==============================\n"
              << " Iterative Solver: Zero guess \n"
              << "==============================\n" << std::endl;
    DebyeJacobiSolve(rhs, potential2, error_array, debye_length, 1E-10);
    fout = fopen("output/potential-iterative-zero.txt", "w");
    potential2.WriteField(fout);
    fclose(fout);
    fout = fopen("output/error-zero.txt", "w");
    WriteArray(fout, error_array, MAX_ITER_NUM);
    fclose(fout);

    // Iterative Solver: Random Guess
    std::cout << "================================\n"
              << " Iterative Solver: Random guess \n"
              << "================================\n" << std::endl;
    for (int i = 0; i < n; i++) potential3(i) = double(rand()) / RAND_MAX - 0.5;
    DebyeJacobiSolve(rhs, potential3, error_array, debye_length, 1E-10);
    fout = fopen("output/potential-iterative-random.txt", "w");
    potential3.WriteField(fout);
    fclose(fout);
    fout = fopen("output/error-random.txt", "w");
    WriteArray(fout, error_array, MAX_ITER_NUM);
    fclose(fout);

    // Iterative Solver: Better Guess
    std::cout << "================================\n"
              << " Iterative Solver: Better guess \n"
              << "================================\n" << std::endl;
    Field3D potential4(rhs);
    for (int i = 0; i < n; i++) potential4(i) = -rhs(i) / (ee / e0 * 1E27);
    DebyeJacobiSolve(rhs, potential4, error_array, debye_length, 1E-10);
    fout = fopen("output/potential-iterative-better.txt", "w");
    potential4.WriteField(fout);
    fclose(fout);
    fout = fopen("output/error-better.txt", "w");
    WriteArray(fout, error_array, MAX_ITER_NUM);
    fclose(fout);

    delete[] error_array;
    return 0;
}

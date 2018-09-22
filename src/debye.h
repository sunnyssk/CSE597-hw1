#ifndef _DEBYE_H_
#define _DEBYE_H_

#define MAX_ITER_NUM (10000)

#include <iostream>
#include "matrix.h"

class Field3D {
public:
    Field3D () : nx_(0), ny_(0), nz_(0), dx_(0.0), dy_(0.0), dz_(0.0), nl_(0), n_(0), buffer_(nullptr) {}
    Field3D (int nx, int ny, int nz, double dx, double dy, double dz) : nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), nl_(nx * ny), n_(nx * ny * nz), buffer_(nullptr) { buffer_ = new double[n_](); }
    Field3D (int nx, int ny, int nz, double dl) : nx_(nx), ny_(ny), nz_(nz), dx_(dl), dy_(dl), dz_(dl), nl_(nx * ny), n_(nx * ny * nz), buffer_(nullptr) { buffer_ = new double[n_](); }
    Field3D (Field3D const & src) : nx_(src.nx_), ny_(src.ny_), nz_(src.nz_), dx_(src.dx_), dy_(src.dy_), dz_(src.dz_), nl_(src.nl_), n_(src.n_), buffer_(nullptr) { buffer_ = new double[n_]; for (int i = 0; i < n_; i++) buffer_[i] = src.buffer_[i]; }
    ~Field3D () { if(buffer_) delete[] buffer_; }
    
    double const & operator () (int i, int j, int k) const { return buffer_[k * nl_ + j * nx_ + i]; }
    double const & operator () (int ind) const { return buffer_[ind]; }
    double & operator () (int i, int j, int k) { return buffer_[k * nl_ + j * nx_ + i]; }
    double & operator () (int ind) { return buffer_[ind]; }
    int const Nx () const { return nx_; }
    int const Ny () const { return ny_; }
    int const Nz () const { return nz_; }
    int const N () const { return n_; }
    double const Dx () const { return dx_; }
    double const Dy () const { return dy_; }
    double const Dz () const { return dz_; }

    void AssignField (double * buffer) { for (int i = 0; i < n_; i++) buffer_[i] = buffer[i]; }

protected:
    int nx_;
    int ny_;
    int nz_;
    double dx_;
    double dy_;
    double dz_;
    int nl_;                                  // number of points in an XY layer
    int n_;
    double * buffer_;
};

class DebyeLUSolver {
public:
    DebyeLUSolver () : pfield_(nullptr), pAmat_(nullptr), pPmat_(nullptr), pLmat_(nullptr), pUmat_(nullptr) {}
    ~DebyeLUSolver () {if(pfield_) delete pfield_; if(pAmat_) delete pAmat_; if(pPmat_) delete pPmat_; if(pLmat_) delete pLmat_; if(pUmat_) delete pUmat_; }

    void GenerateSolverMatrix (Field3D const & rhs, double debye_length);
    void SolverMatrixDecompose ();
    void RhsInput (Field3D const & field);
    void LUSolve (Field3D & res_container);

protected:
    MatD *pfield_;
    MatD *pAmat_;
    MatD *pPmat_;
    MatD *pLmat_;
    MatD *pUmat_;

    DebyeLUSolver(DebyeLUSolver const & src) {}
};

void ApplDirichletCond (Field3D & field);

void DebyeJacobiSolve (Field3D const & rhs, Field3D & potential, double debye_length, double err_threshold=1E-10);

#endif /* _DEBYE_H_ */
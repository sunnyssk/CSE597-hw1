#include "debye.h"

void ApplDirichletCond (Field3D & field) {
    int nx = field.Nx(), ny = field.Ny(), nz = field.Nz();
    for (int j = 0; j < ny; j++) {
        for (int k = 0; k < nz; k++) {
            field(0, j, k) = 0.0;
            field(nx - 1, j, k) = 0.0;
        }
    }
    for (int k = 0; k < nz; k++) {
        for (int i = 0; i < nx; i++) {
            field(i, 0, k) = 0.0;
            field(i, ny - 1, k) = 0.0;
        }
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            field(i, j, 0) = 0.0;
            field(i, j, nz - 1) = 0.0;
        }
    }
}

int DebyeJacobiSolve (Field3D const & rhs, Field3D & potential, double * iter_err_array, double debye_length, double err_threshold) {
    int nx = rhs.Nx(), ny = rhs.Ny(), nz = rhs.Nz();
    Field3D aux_field(nx, ny, nz, rhs.Dx(), rhs.Dy(), rhs.Dz());
    double dx2 = rhs.Dx(), dy2 = rhs.Dy(), dz2 = rhs.Dz(), ldi2 = 1 / (debye_length * debye_length);
    dx2 = dx2 * dx2;
    dy2 = dy2 * dy2;
    dz2 = dz2 * dz2;
    Field3D *prev = &potential, *next = &aux_field;
    int iter_cnt = 0;
    double errmax = 0.0, denominator = 2 / dx2 + 2 / dy2 + 2 / dz2 + ldi2;

    ApplDirichletCond(potential);
    ApplDirichletCond(aux_field);
    do {
        errmax = 0.0;
        if (iter_cnt++ >= MAX_ITER_NUM - 1) {
            std::cout << "Iteration number exceeds limit! Exit automatically." << std::endl;
            break;
        }
        for (int i = 1; i < nx - 1; i++) {
            for (int j = 1; j < ny - 1; j++) {
                for (int k = 1; k < nz - 1; k++) {
                    (*next)(i, j, k) = (((*prev)(i + 1, j, k) + (*prev)(i - 1, j, k)) / dx2 + ((*prev)(i, j + 1, k) + (*prev)(i, j - 1, k)) / dy2 + ((*prev)(i, j, k + 1) + (*prev)(i, j, k - 1)) / dz2 - rhs(i, j, k)) / denominator;
                    errmax = errmax > fabs((*next)(i, j, k) - (*prev)(i, j, k)) ? errmax : fabs((*next)(i, j, k) - (*prev)(i, j, k));
                }
            }
        }
        Field3D *tmp = prev;
        prev = next;
        next = tmp;
        std::cout << "Iteration round #" << iter_cnt << ", maximum error: " << errmax << "\n";
        iter_err_array[iter_cnt - 1] = errmax;
    } while (errmax >= err_threshold || &potential != next);                    // only exits iteration when the original field is updated
    char *buffer = new char[256];
    MemSizeOutput(buffer);
    delete[] buffer;
    return iter_cnt;
}

void DebyeLUSolver::GenerateSolverMatrix (Field3D const & rhs, double debye_length) {
    int nx = rhs.Nx(), ny = rhs.Ny(), nz = rhs.Nz(), n = nx * ny * nz;
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    double dx = rhs.Dx(), dy = rhs.Dy(), dz = rhs.Dz();
    double dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz, ldi2 = 1 / (debye_length * debye_length);
    if (pAmat_) delete pAmat_;
    pAmat_ = new MatD(ntrim, ntrim);
    MatD &A = *pAmat_;

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++) {
        for (int j = 1; j <= nytrim; j++) {
            for (int k = 1; k <= nztrim; k++) {
                A(TRI(i, j, k), TRI(i, j, k)) = -(2 / dx2 + 2 / dy2 + 2 / dz2 + ldi2);
                if(i > 1) A(TRI(i, j, k), TRI(i - 1, j, k)) = 1 / dx2;
                if(j > 1) A(TRI(i, j, k), TRI(i, j - 1, k)) = 1 / dy2;
                if(k > 1) A(TRI(i, j, k), TRI(i, j, k - 1)) = 1 / dz2;
                if(i < nxtrim) A(TRI(i, j, k), TRI(i + 1, j, k)) = 1 / dx2;
                if(j < nytrim) A(TRI(i, j, k), TRI(i, j + 1, k)) = 1 / dy2;
                if(k < nztrim) A(TRI(i, j, k), TRI(i, j, k + 1)) = 1 / dz2;
            }
        }
    }
    #undef TRI
}

void DebyeLUSolver::LUSolve (Field3D & res_container) {
    int nx = res_container.Nx(), ny = res_container.Ny(), nz = res_container.Nz(), n = nx * ny * nz;
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    MatD res_trim(ntrim, 1), LPb(ntrim, 1), Pb(ntrim, 1);
    Pb = (*pPmat_) * (*pfield_);
    pLmat_->LSolve(Pb, LPb);
    pUmat_->USolve(LPb, res_trim);
    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) res_container(i, j, k) = res_trim(TRI(i, j, k));
    #undef TRI
}

void DebyeLUSolver::RhsInput (Field3D const & field) {
    int nx = field.Nx(), ny = field.Ny(), nz = field.Nz(), n = nx * ny * nz;
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    if (pfield_) delete pfield_;
    pfield_ = new MatD(ntrim, 1);
    MatD &rhs = *pfield_;

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) rhs(TRI(i, j, k)) = field(i, j, k);
    #undef TRI
}

void DebyeLUSolver::SolverMatrixDecompose () {
    int r = pAmat_->Rows(), c = pAmat_->Cols();
    if(pPmat_) delete pPmat_;
    if(pLmat_) delete pLmat_;
    if(pUmat_) delete pUmat_;
    pPmat_ = new MatD(r, c);
    pLmat_ = new MatD(r, c);
    pUmat_ = new MatD(r, c);
    pAmat_->PLUDecomposition(*pPmat_, *pLmat_, *pUmat_);
}

void Field3D::WriteField (FILE * output_file) {
    fprintf(output_file, "%d %d %d\n", nx_, ny_, nz_);
    fprintf(output_file, "%lE %lE %lE\n", dx_, dy_, dz_);
    for (int k = 0; k < nz_; k++)
        for (int j = 0; j < ny_; j++)
            for (int i = 0; i < nx_; i++) fprintf(output_file, "%lE\n", (*this)(i, j, k));
}

void WriteArray (FILE * output_file, double * array, int length) {
    fprintf(output_file, "%d\n", length);
    for (int i = 0; i < length; i++) fprintf(output_file, "%lE\n", array[i]);
}
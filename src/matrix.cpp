#include "matrix.h"

template <typename T> void Matrix<T>::Copy (Matrix<T> const & src) {
    if(!mapping_only_) delete[] data_;      // releases current data array first
    rows_ = src.rows_;
    cols_ = src.cols_;
    mapping_only_ = false;
    int size = rows_ * cols_;
    data_ = new T[size];
    for (int i = 0; i < size; i++) data_[i] = src.data_[i]; 
}

template <typename T> void Matrix<T>::LUDecomposition (Matrix<T> & l_container, Matrix<T> & u_container) {
    if(rows_ != cols_) throw "ERROR: LU decomposition not applying on a square matrix.";
    // else:
    u_container.Copy(*this);
    l_container.Reset(rows_, cols_);
    for (int k = 0; k < rows_ - 1; k++) {
        l_container(k, k) = T(1);
        for (int i = k + 1; i < rows_; i++) {
            T ratio = -u_container(k, k) / u_container(i, k);
            u_container(i, k) = T(0);
            l_container(i, k) = -ratio;
            for(int j = k + 1; j < cols_; j++) u_container(i, j) += ratio * u_container(k, j);
        }
    }
}

template <typename T> void Matrix<T>::Print (std::ostream & output_stream) const {
    for (int i = 0; i < rows_; i++) {
        for (int j = 0; j < cols_ - 1; j++) {
            output_stream << data_[i + j * rows_] << ", ";
        }
        output_stream << data_[i + (cols_ - 1) * rows_] << "\n";
    }
}

template <typename T> void Matrix<T>::Reset (int rows, int cols) {
    if(!mapping_only_) delete[] data_;
    rows_ = rows;
    cols_ = cols;
    mapping_only_ = false;
    data_ = new T[rows_ * cols_];
}

template <typename T> void Matrix<T>::Reshape (int new_rows, int new_cols) {
    if (new_rows * new_cols != rows_ * cols_) throw "ERROR: Matrix size changed during reshape.";
    // else:
    rows_ = new_rows;
    cols_ = new_cols;
}

template <typename T> void Matrix<T>::SwapRows (int row1, int row2) {
    T tmp;
    for (int i = 0; i < cols_; i++) {
        int ind1 = i * rows_ + row1, ind2 = i * rows_ + row2;
        tmp = data_[ind1];
        data_[ind1] = data_[ind2];
        data_[ind2] = tmp;
    }
}

template <typename T> void Matrix<T>::SwapCols (int col1, int col2) {
    T tmp;
    for (int i = 0; i < rows_; i++) {
        int ind1 = col1 * rows_ + i, ind2 = col2 * rows_ + i;
        tmp = data_[ind1];
        data_[ind1] = data_[ind2];
        data_[ind2] = tmp;
    }
}

template <typename T> void Matrix<T>::Multiply (Matrix<T> const & opr1, Matrix<T> const & opr2, Matrix<T> & res) {
    if (opr1.cols_ != opr2.rows_ || opr1.rows_ != res.rows_ || opr2.cols_ != res.cols_) throw "ERROR: Matrix operand dimensions do not match.";
    // else:
    memset(res.data_, 0, sizeof(T) * res.Size());
    for (int i = 0; i < opr1.rows_; i++)
        for (int j = 0; j < opr2.cols_; j++)
            for (int k = 0; k < opr1.cols_; k++) res(i, j) += opr1(i, k) * opr2(k, j);
}
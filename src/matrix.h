#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <iostream>
#include <cmath>
#include <cstring>
#include <exception>

template <typename T> class Matrix {
public:
    Matrix (int rows, int cols) : rows_(rows), cols_(cols), mapping_only_(false), data_(nullptr) { data_ = new T[rows_ * cols_]; }
    Matrix (int rows, int cols, T* buffer, bool mapping_only=false) : rows_(rows), cols_(cols), mapping_only_(mapping_only), data_(nullptr) {
        int size = rows_ * cols_;
        data_ = new T[size];
        if (mapping_only) data_ = buffer;
        else for (int i = 0; i < size; i++) data_[i] = buffer[i];
    }
    Matrix(Matrix const & src) : rows_(src.rows_), cols_(src.cols_), mapping_only_(false), data_(nullptr) {
        int size = rows_ * cols_;
        data_ = new T[size];
        for(int i = 0; i < size; i++) data_[i] = src(i);
    }
    ~Matrix () { if (!mapping_only_) delete[] data_; }

    inline int Rows () const { return rows_; }
    inline int Cols () const { return cols_; }
    inline int Size () const { return rows_ * cols_; }
    inline T& operator () (int index) { return data_[index]; }
    inline T& operator () (int row, int col) { return data_[col * rows_ + row]; }
    const inline T& operator () (int index) const { return data_[index]; }
    const inline T& operator () (int row, int col) const { return data_[col * rows_ + row]; }

    // Copies the info from src to current object.
    void Copy(Matrix const & src);
    Matrix & operator = (Matrix const & src) { this->Copy(src); return *this; }
    // Resets a matrix with the assigned shape.
    void Reset (int rows, int cols);
    // Reshapes matrix to a new shape. Number of elements must be invariant.
    void Reshape (int new_rows, int new_cols);
    // Swaps two rows of a matrix.
    void SwapRows (int row1, int row2);
    // Swaps two columns of a matrix.
    void SwapCols (int col1, int col2);

    // Prints a matrix to std::out.
    void Print (std::ostream & output_stream=std::cout) const;
    friend std::ostream & operator << (std::ostream & output_stream, Matrix matrix) {
        matrix.Print(output_stream);
        return output_stream;
    }

    // Carries out LU decomposition on a matrix and saves the result into matrices l_container and u_container.
    void LUDecomposition (Matrix & l_container, Matrix & u_container);

    // Matrix multiplications.
    static void Multiply (Matrix const & opr1, Matrix const & opr2, Matrix & res);
    Matrix operator * (Matrix const & opr2) { Matrix result(this->Rows(), opr2.Cols()); Multiply((*this), opr2, result); return result; }

protected:
    int rows_;             // rows of matrix
    int cols_;             // columns of matrix
    bool mapping_only_;    // whether data in the matrix is only a mapping of existing array or not
    T* data_;              // data buffer

    Matrix ();             // forbids empty constructor
};

template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;

typedef Matrix<int> MatI;
typedef Matrix<float> MatF;
typedef Matrix<double> MatD;

#endif /* _MATRIX_H_ */
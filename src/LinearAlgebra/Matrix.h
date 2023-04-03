#pragma once

#include "LinearAlgebra/NumericArray.h"

template <typename T>
class Matrix final : public NumericArray<T> {
 public:
  Matrix(int rows, int cols) : NumericArray<T>(0, cols, rows) {}

  T& operator()(int i, int j) { return this->buffer_[i * cols() + j]; }
  T operator()(int i, int j) const { return this->buffer_[i * cols() + j]; }

  int rows() const { return this->stride_; }
  int cols() const { return this->size_; }
};

void gemv(double alpha, const Matrix<double>& A, const NumericArray<double>& x, double beta,
          NumericArray<double>& y);

void ger(double alpha, const NumericArray<double>& x, const NumericArray<double>& y,
         Matrix<double>& A);

void gemm(double alpha, const Matrix<double>& a, const Matrix<double>& b, double beta,
          Matrix<double>& c);

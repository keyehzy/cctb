#pragma once

#include "LinearAlgebra/NumericArray.h"

template <typename T>
class Matrix final : public NumericArray<T> {
 public:
  Matrix(int rows, int cols) : NumericArray<T>(0, cols, rows) {}

  T& operator()(int i, int j) { return this->buffer_[j * rows() + i]; }
  T operator()(int i, int j) const { return this->buffer_[j * rows() + i]; }

  int rows() const { return this->stride_; }
  int cols() const { return this->size_; }
};

template <typename T>
ALWAYS_INLINE void gemv(T alpha, const Matrix<T>& a, const NumericArray<T>& x, T beta,
                        NumericArray<T>& y) {
  for (int i = 0; i < a.rows(); ++i) {
    double result = 0;
    for (int j = 0; j < a.cols(); ++j) {
      result += a(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
}

template <typename T>
ALWAYS_INLINE void ger(T alpha, const NumericArray<T>& x, const NumericArray<T>& y, Matrix<T>& a) {
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      a(i, j) += alpha * x[i] * y[j];
    }
  }
}

template <typename T>
ALWAYS_INLINE void gemm(T alpha, const Matrix<T>& a, const Matrix<T>& b, T beta, Matrix<T>& c) {
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < b.cols(); ++j) {
      double result = 0;
      for (int k = 0; k < a.cols(); ++k) {
        result += a(i, k) * b(k, j);
      }
      c(i, j) = alpha * result + beta * c(i, j);
    }
  }
}

#pragma once

#include "LinearAlgebra/NumericArray.h"
#include "LinearAlgebra/Vector.h"

class Matrix final : public NumericArray<double> {
 public:
  Matrix(int rows, int cols) : NumericArray<double>(rows * cols), rows_(rows), cols_(cols) {}

  double& operator()(int i, int j) { return contiguous_->at(i * cols_ + j); }
  double operator()(int i, int j) const { return contiguous_->at(i * cols_ + j); }

  int rows() const { return rows_; }
  int cols() const { return cols_; }

  friend ALWAYS_INLINE void gemv(double alpha, const Matrix& a, const Vector& x, double beta,
                                 Vector& y);
  friend ALWAYS_INLINE void ger(double alpha, const Vector& x, const Vector& y, Matrix& a);
  friend ALWAYS_INLINE void gemm(double alpha, const Matrix& a, const Matrix& b, double beta,
                                 Matrix& c);

 private:
  int rows_;
  int cols_;
};

ALWAYS_INLINE void gemv(double alpha, const Matrix& a, const Vector& x, double beta, Vector& y) {
  for (int i = 0; i < a.rows(); ++i) {
    double result = 0;
    for (int j = 0; j < a.cols(); ++j) {
      result += a(i, j) * x[j];
    }
    y[i] = alpha * result + beta * y[i];
  }
}

ALWAYS_INLINE void ger(double alpha, const Vector& x, const Vector& y, Matrix& a) {
  for (int i = 0; i < a.rows(); ++i) {
    for (int j = 0; j < a.cols(); ++j) {
      a(i, j) += alpha * x[i] * y[j];
    }
  }
}

ALWAYS_INLINE void gemm(double alpha, const Matrix& a, const Matrix& b, double beta, Matrix& c) {
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

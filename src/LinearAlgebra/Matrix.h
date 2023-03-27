#pragma once

#include "LinearAlgebra/NumericArray.h"
#include "LinearAlgebra/Vector.h"

class Matrix final : public NumericArray<double> {
 public:
  Matrix(int rows, int cols) : NumericArray<double>(rows * cols), rows_(rows), cols_(cols) {}

  double& operator()(int i, int j) { return this->data_[i * cols_ + j]; }
  double operator()(int i, int j) const { return this->data_[i * cols_ + j]; }

  ALWAYS_INLINE Vector row(int i) const;
  ALWAYS_INLINE Vector col(int j) const;

  // matrix vector multiplication
  ALWAYS_INLINE void gemv(double alpha, const Vector& x, double beta) const;

  // A := alpha*x*y' + A
  ALWAYS_INLINE void ger(double alpha, const Vector& x, const Vector& y) const;

 private:
  int rows_;
  int cols_;
};

ALWAYS_INLINE Vector Matrix::row(int i) const { 
    return (*this).into(i * this->cols_, this->cols_, 1);
}

ALWAYS_INLINE Vector Matrix::col(int j) const { 
    return (*this).into(j, this->rows_, this->cols_);
}

ALWAYS_INLINE void Matrix::gemv(double alpha, const Vector& x, double beta) const {
  for (int i = 0; i < this->rows_; ++i) {
    double result = 0;
    for (int j = 0; j < this->cols_; ++j) {
      result += this->data_[i * this->cols_ + j] * x[j];
    }
    this->data_[i] = alpha * result + beta * this->data_[i];
  }
}

ALWAYS_INLINE void Matrix::ger(double alpha, const Vector& x, const Vector& y) const {
  for (int i = 0; i < this->rows_; ++i) {
    for (int j = 0; j < this->cols_; ++j) {
      this->data_[i * this->cols_ + j] += alpha * x[i] * y[j];
    }
  }
}

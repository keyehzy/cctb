#pragma once

#include "LinearAlgebra/NumericArray.h"

struct Vector final : public NumericArray<double> {
  Vector(int size) : NumericArray<double>(size) {}
  Vector(int start, int size, int stride) : NumericArray<double>(start, size, stride) {}
  Vector(const Vector& other) : NumericArray<double>(other) {}
  Vector(const NumericArray<double>& other) : NumericArray<double>(other) {}

  ALWAYS_INLINE double dot(const Vector& other) const;
  ALWAYS_INLINE double norm() const;
  ALWAYS_INLINE void scale(double alpha);
  ALWAYS_INLINE double sum() const;
};

ALWAYS_INLINE double Vector::dot(const Vector& other) const {
  double result = 0;
  for (int i = 0; i < size(); ++i) {
    result += (*this)[i] * other[i];
  }
  return result;
}

ALWAYS_INLINE double Vector::norm() const {
  double result = 0;
  for (int i = 0; i < size(); ++i) {
    result += (*this)[i] * (*this)[i];
  }
  return std::sqrt(result);
}

ALWAYS_INLINE void Vector::scale(double alpha) {
  for (int i = 0; i < size(); ++i) {
    (*this)[i] *= alpha;
  }
}

ALWAYS_INLINE double Vector::sum() const {
  double result = 0;
  for (int i = 0; i < size(); ++i) {
    result += (*this)[i];
  }
  return result;
}

ALWAYS_INLINE void axpy(double alpha, const Vector& x, Vector& y) {
  for (int i = 0; i < x.size(); ++i) {
    y[i] += alpha * x[i];
  }
}
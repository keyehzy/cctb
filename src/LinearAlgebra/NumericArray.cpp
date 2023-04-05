#include "LinearAlgebra/NumericArray.h"

#include <complex>

#include "LinearAlgebra/BlasImpl.h"

template <>
double NumericArray<double>::dot(const NumericArray<double>& other) const {
  double result;
  ::dot(*this, other, &result);
  return result;
}

template <>
double NumericArray<double>::norm() const {
  double result;
  nrm2(*this, &result);
  return result;
}

template <>
double NumericArray<double>::sum() const {
  double result;
  asum(*this, &result);
  return result;
}

template <>
void NumericArray<double>::scale(double alpha) {
  scal(alpha, *this);
}

template <>
std::complex<double> NumericArray<std::complex<double>>::dot(
    const NumericArray<std::complex<double>>& other) const {
  std::complex<double> result;
  ::dot(*this, other, &result);
  return result;
}

template <>
double NumericArray<std::complex<double>>::norm() const {
  std::complex<double> result;
  ::dot(*this, *this, &result);
  return result.real();
}

template <>
std::complex<double> NumericArray<std::complex<double>>::sum() const {
  std::complex<double> result;
  asum(*this, &result);
  return result;
}

template <>
void NumericArray<std::complex<double>>::scale(std::complex<double> alpha) {
  scal(alpha, *this);
}

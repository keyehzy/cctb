#pragma once

#include <complex>

#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericArray.h"

void geev(const Matrix<double>& a, NumericArray<std::complex<double>>& w,
          Matrix<std::complex<double>>& v);

void geev(const Matrix<std::complex<double>>& a, NumericArray<std::complex<double>>& w,
          Matrix<std::complex<double>>& v);

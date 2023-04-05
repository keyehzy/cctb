#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>

#include "LinearAlgebra/BlasImpl.h"
#include "LinearAlgebra/LapackImpl.h"
#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericArray.h"
#include "TestMatchers.h"

TEST_CASE("gemv", "[Matrix]") {
  Matrix<double> A(3, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = 4.0;
  A(1, 1) = 5.0;
  A(1, 2) = 6.0;
  A(2, 0) = 7.0;
  A(2, 1) = 8.0;
  A(2, 2) = 9.0;
  NumericArray<double> x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  NumericArray<double> y(3);
  y[0] = 4.0;
  y[1] = 5.0;
  y[2] = 6.0;
  gemv(2.0, A, x, 3.0, y);
  REQUIRE_THAT(y[0], Catch::Matchers::WithinAbs(40.0, 1e-10));
  REQUIRE_THAT(y[1], Catch::Matchers::WithinAbs(79.0, 1e-10));
  REQUIRE_THAT(y[2], Catch::Matchers::WithinAbs(118.0, 1e-10));
}

TEST_CASE("ger", "[Matrix]") {
  Matrix<double> A(3, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = 4.0;
  A(1, 1) = 5.0;
  A(1, 2) = 6.0;
  A(2, 0) = 7.0;
  A(2, 1) = 8.0;
  A(2, 2) = 9.0;
  NumericArray<double> x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  NumericArray<double> y(3);
  y[0] = 4.0;
  y[1] = 5.0;
  y[2] = 6.0;
  ger(2.0, x, y, A);
  REQUIRE_THAT(A(0, 0), Catch::Matchers::WithinAbs(9.0, 1e-10));
  REQUIRE_THAT(A(0, 1), Catch::Matchers::WithinAbs(12.0, 1e-10));
  REQUIRE_THAT(A(0, 2), Catch::Matchers::WithinAbs(15.0, 1e-10));
  REQUIRE_THAT(A(1, 0), Catch::Matchers::WithinAbs(20.0, 1e-10));
  REQUIRE_THAT(A(1, 1), Catch::Matchers::WithinAbs(25.0, 1e-10));
  REQUIRE_THAT(A(1, 2), Catch::Matchers::WithinAbs(30.0, 1e-10));
  REQUIRE_THAT(A(2, 0), Catch::Matchers::WithinAbs(31.0, 1e-10));
  REQUIRE_THAT(A(2, 1), Catch::Matchers::WithinAbs(38.0, 1e-10));
  REQUIRE_THAT(A(2, 2), Catch::Matchers::WithinAbs(45.0, 1e-10));
}

TEST_CASE("gemm", "[Matrix]") {
  Matrix<double> A(3, 3);
  A(0, 0) = 1.0;
  A(0, 1) = 2.0;
  A(0, 2) = 3.0;
  A(1, 0) = 4.0;
  A(1, 1) = 5.0;
  A(1, 2) = 6.0;
  A(2, 0) = 7.0;
  A(2, 1) = 8.0;
  A(2, 2) = 9.0;
  Matrix<double> B(3, 3);
  B(0, 0) = 1.0;
  B(0, 1) = 2.0;
  B(0, 2) = 3.0;
  B(1, 0) = 4.0;
  B(1, 1) = 5.0;
  B(1, 2) = 6.0;
  B(2, 0) = 7.0;
  B(2, 1) = 8.0;
  B(2, 2) = 9.0;
  Matrix<double> C(3, 3);
  C(0, 0) = 1.0;
  C(0, 1) = 2.0;
  C(0, 2) = 3.0;
  C(1, 0) = 4.0;
  C(1, 1) = 5.0;
  C(1, 2) = 6.0;
  C(2, 0) = 7.0;
  C(2, 1) = 8.0;
  C(2, 2) = 9.0;
  gemm(2.0, A, B, 3.0, C);
  REQUIRE_THAT(C(0, 0), Catch::Matchers::WithinAbs(63.0, 1e-10));
  REQUIRE_THAT(C(0, 1), Catch::Matchers::WithinAbs(78.0, 1e-10));
  REQUIRE_THAT(C(0, 2), Catch::Matchers::WithinAbs(93.0, 1e-10));
  REQUIRE_THAT(C(1, 0), Catch::Matchers::WithinAbs(144.0, 1e-10));
  REQUIRE_THAT(C(1, 1), Catch::Matchers::WithinAbs(177.0, 1e-10));
  REQUIRE_THAT(C(1, 2), Catch::Matchers::WithinAbs(210.0, 1e-10));
  REQUIRE_THAT(C(2, 0), Catch::Matchers::WithinAbs(225.0, 1e-10));
  REQUIRE_THAT(C(2, 1), Catch::Matchers::WithinAbs(276.0, 1e-10));
  REQUIRE_THAT(C(2, 2), Catch::Matchers::WithinAbs(327.0, 1e-10));
}

TEST_CASE("dgeev", "[Matrix]") {
  {
    Matrix<double> A(2, 2);
    A(0, 1) = 1.0;
    A(1, 0) = 1.0;
    NumericArray<std::complex<double>> w(2);
    Matrix<std::complex<double>> v(2, 2);
    geev(A, w, v);
    REQUIRE_THAT(w[0], ApproxEqualComplex(std::complex<double>(1.0, 0.0)));
    REQUIRE_THAT(w[1], ApproxEqualComplex(std::complex<double>(-1.0, 0.0)));
    NumericArray<std::complex<double>> x(2), y(2);
    x[0] = v(0, 0);
    x[1] = v(0, 1);
    y[0] = v(1, 0);
    y[1] = v(1, 1);
    REQUIRE_THAT(x.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(y.norm(), Catch::Matchers::WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(x.dot(y), ApproxEqualComplex(std::complex<double>(0.0, 0.0)));
    NumericArray<std::complex<double>> z(2);
    std::complex<double> alpha = 1.0;
    std::complex<double> beta = 0.0;
    Matrix<std::complex<double>> A_c(2, 2);
    A_c(0, 1) = 1.0;
    A_c(1, 0) = 1.0;
    gemv(alpha, A_c, x, beta, z);
    REQUIRE_THAT(z[0], ApproxEqualComplex(w[0] * x[0]));
    REQUIRE_THAT(z[1], ApproxEqualComplex(w[0] * x[1]));
    gemv(alpha, A_c, y, beta, z);
    REQUIRE_THAT(z[0], ApproxEqualComplex(w[1] * y[0]));
    REQUIRE_THAT(z[1], ApproxEqualComplex(w[1] * y[1]));
  }
}

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>
#include <iostream>

#include "LinearAlgebra/BlasImpl.h"
#include "LinearAlgebra/LapackImpl.h"
#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericArray.h"

struct ApproxEqualComplex : Catch::Matchers::MatcherGenericBase {
 public:
  ApproxEqualComplex(std::complex<double> v) : v_(v) {}

  bool match(std::complex<double> const &v) const {
    return std::abs(v.real() - v_.real()) < 1e-10 && std::abs(v.imag() - v_.imag()) < 1e-10;
  }

  std::string describe() const override { return "ApproxEqualComplex"; }

 private:
  std::complex<double> v_;
};

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
  REQUIRE(y[0] == 40.0);
  REQUIRE(y[1] == 79.0);
  REQUIRE(y[2] == 118.0);
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
  REQUIRE(A(0, 0) == 9.0);
  REQUIRE(A(0, 1) == 12.0);
  REQUIRE(A(0, 2) == 15.0);
  REQUIRE(A(1, 0) == 20.0);
  REQUIRE(A(1, 1) == 25.0);
  REQUIRE(A(1, 2) == 30.0);
  REQUIRE(A(2, 0) == 31.0);
  REQUIRE(A(2, 1) == 38.0);
  REQUIRE(A(2, 2) == 45.0);
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
  REQUIRE(C(0, 0) == 63.0);
  REQUIRE(C(0, 1) == 78.0);
  REQUIRE(C(0, 2) == 93.0);
  REQUIRE(C(1, 0) == 144.0);
  REQUIRE(C(1, 1) == 177.0);
  REQUIRE(C(1, 2) == 210.0);
  REQUIRE(C(2, 0) == 225.0);
  REQUIRE(C(2, 1) == 276.0);
  REQUIRE(C(2, 2) == 327.0);
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

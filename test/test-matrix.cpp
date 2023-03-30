#include <catch2/catch_test_macros.hpp>

#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/NumericArray.h"

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

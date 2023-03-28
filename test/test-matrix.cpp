#include <catch2/catch_test_macros.hpp>

#include "LinearAlgebra/Matrix.h"
#include "LinearAlgebra/Vector.h"

TEST_CASE("gemv", "[Matrix]") {
  Matrix A(3, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;
  A(2, 0) = 7;
  A(2, 1) = 8;
  A(2, 2) = 9;
  Vector x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  Vector y(3);
  y[0] = 4;
  y[1] = 5;
  y[2] = 6;
  gemv(2, A, x, 3, y);
  REQUIRE(y[0] == 40);
  REQUIRE(y[1] == 79);
  REQUIRE(y[2] == 118);
}

TEST_CASE("ger", "[Matrix]") {
  Matrix A(3, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;
  A(2, 0) = 7;
  A(2, 1) = 8;
  A(2, 2) = 9;
  Vector x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  Vector y(3);
  y[0] = 4;
  y[1] = 5;
  y[2] = 6;
  ger(2, x, y, A);
  REQUIRE(A(0, 0) == 9);
  REQUIRE(A(0, 1) == 12);
  REQUIRE(A(0, 2) == 15);
  REQUIRE(A(1, 0) == 20);
  REQUIRE(A(1, 1) == 25);
  REQUIRE(A(1, 2) == 30);
  REQUIRE(A(2, 0) == 31);
  REQUIRE(A(2, 1) == 38);
  REQUIRE(A(2, 2) == 45);
}

TEST_CASE("gemm", "[Matrix]") {
  Matrix A(3, 3);
  A(0, 0) = 1;
  A(0, 1) = 2;
  A(0, 2) = 3;
  A(1, 0) = 4;
  A(1, 1) = 5;
  A(1, 2) = 6;
  A(2, 0) = 7;
  A(2, 1) = 8;
  A(2, 2) = 9;
  Matrix B(3, 3);
  B(0, 0) = 1;
  B(0, 1) = 2;
  B(0, 2) = 3;
  B(1, 0) = 4;
  B(1, 1) = 5;
  B(1, 2) = 6;
  B(2, 0) = 7;
  B(2, 1) = 8;
  B(2, 2) = 9;
  Matrix C(3, 3);
  C(0, 0) = 1;
  C(0, 1) = 2;
  C(0, 2) = 3;
  C(1, 0) = 4;
  C(1, 1) = 5;
  C(1, 2) = 6;
  C(2, 0) = 7;
  C(2, 1) = 8;
  C(2, 2) = 9;
  gemm(2, A, B, 3, C);
  REQUIRE(C(0, 0) == 63);
  REQUIRE(C(0, 1) == 78);
  REQUIRE(C(0, 2) == 93);
  REQUIRE(C(1, 0) == 144);
  REQUIRE(C(1, 1) == 177);
  REQUIRE(C(1, 2) == 210);
  REQUIRE(C(2, 0) == 225);
  REQUIRE(C(2, 1) == 276);
  REQUIRE(C(2, 2) == 327);
}
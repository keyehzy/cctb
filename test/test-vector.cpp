#include <catch2/catch_test_macros.hpp>

#include "LinearAlgebra/Vector.h"

TEST_CASE("Vector::dot", "[Vector]") {
  {
    Vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    Vector y(3);
    y[0] = 4;
    y[1] = 5;
    y[2] = 6;
    REQUIRE(x.dot(y) == 32);
  }
}

TEST_CASE("Vector::norm", "[Vector]") {
  {
    Vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    REQUIRE(x.norm() == std::sqrt(14));
  }
}

TEST_CASE("Vector::scale", "[Vector]") {
  {
    Vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    x.scale(2);
    REQUIRE(x[0] == 2);
    REQUIRE(x[1] == 4);
    REQUIRE(x[2] == 6);
  }
}

TEST_CASE("Vector::sum", "[Vector]") {
  {
    Vector x(3);
    x[0] = 1;
    x[1] = 2;
    x[2] = 3;
    REQUIRE(x.sum() == 6);
  }
}

TEST_CASE("axpy", "[Vector]") {
  Vector x(3);
  x[0] = 1;
  x[1] = 2;
  x[2] = 3;
  Vector y(3);
  y[0] = 4;
  y[1] = 5;
  y[2] = 6;
  axpy(2, x, y);
  REQUIRE(y[0] == 6);
  REQUIRE(y[1] == 9);
  REQUIRE(y[2] == 12);
}
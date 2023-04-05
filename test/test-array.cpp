#include <catch2/catch_test_macros.hpp>

#include "LinearAlgebra/BlasImpl.h"

TEST_CASE("axpy", "[NumericArray]") {
  NumericArray<double> x(3);
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  NumericArray<double> y(3);
  y[0] = 4.0;
  y[1] = 5.0;
  y[2] = 6.0;
  axpy(2.0, x, y);
  REQUIRE(y[0] == 6.0);
  REQUIRE(y[1] == 9.0);
  REQUIRE(y[2] == 12.0);
}

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
  REQUIRE_THAT(y[0], Catch::Matchers::WithinAbs(6.0, 1e-10));
  REQUIRE_THAT(y[1], Catch::Matchers::WithinAbs(9.0, 1e-10));
  REQUIRE_THAT(y[2], Catch::Matchers::WithinAbs(12.0, 1e-10));
}

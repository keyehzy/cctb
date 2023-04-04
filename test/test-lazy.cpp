#include <catch2/catch_test_macros.hpp>

#include "LinearAlgebra/Lazy.h"
#include "LinearAlgebra/NumericArray.h"

TEST_CASE("Lazy", "[linear_algebra]") {
  {
    double a = 1;
    double b = 45;
    LazyOp c = AddOp(a, b);
    double d = c;
    REQUIRE(d == 46);
  }

  {
    LazyOp a = AddOp(1, 45);
    LazyOp b = AddOp(3, 12);
    double c = a + b;
    REQUIRE(c == 62);
  }
}

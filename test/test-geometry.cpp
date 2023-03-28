#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "Geometry/Line.h"

TEST_CASE("Line::intercect", "[Line]") {
  {
    Line line1(Vec<double>{0, 0}, Vec<double>{1, 1});
    Line line2(Vec<double>{0, 1}, Vec<double>{1, 0});
    Vec<double> intercect = line1.intercect(line2);
    REQUIRE(intercect == Vec<double>{0.5, 0.5});
  }

  {
    Line line1(Vec<double>{-2.22, -2.53}, Vec<double>{5.04, 2.17});
    Line line2(Vec<double>{-0.76, 4.99}, Vec<double>{6.22, -1.85});
    Vec<double> intercect = line1.intercect(line2);
    REQUIRE_THAT(intercect[0], Catch::Matchers::WithinAbs(3.28, 0.01));
  }
}

TEST_CASE("Line::orthogonal", "[Line]") {
  {
    Line line1(Vec<double>{0, 0}, Vec<double>{1, 1});
    Line line2(Vec<double>{0, 1}, Vec<double>{1, 0});
    REQUIRE(line1.orthogonal(line2));
  }

  {
    Line line1(Vec<double>{1.5, 0.7}, Vec<double>{0.55, -1.3});
    Line line2(Vec<double>{1.025, -0.3}, Vec<double>{0.393421, 0});
    REQUIRE(line1.orthogonal(line2));
  }
}

TEST_CASE("Line::perpendicular_bisector", "[Line]") {
  {
    Line line(Vec<double>{0, 0}, Vec<double>{1, 1});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == Vec<double>{0.5, 0.5});
    REQUIRE(line.orthogonal(bisector));
  }

  {
    Line line(Vec<double>{0, 0}, Vec<double>{2.0 / 3.0, 2.0 / std::sqrt(3.0)});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == line.midpoint());
    REQUIRE(line.orthogonal(bisector));
  }
}

TEST_CASE("Line::midpoint", "[Line]") {
  {
    Line line(Vec<double>{0, 0}, Vec<double>{1, 1});
    REQUIRE(line.midpoint() == Vec<double>{0.5, 0.5});
  }

  {
    Line line(Vec<double>{0, 0}, Vec<double>{2.0 / 3.0, 2.0 / std::sqrt(3.0)});
    REQUIRE(line.midpoint() == Vec<double>{1.0 / 3.0, 1.0 / std::sqrt(3.0)});
  }
}
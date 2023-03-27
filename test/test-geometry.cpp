#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "Geometry/Line.h"

TEST_CASE("Line::intercect", "[Line]") {
  {
    Line line1(Vec<float>{0, 0}, Vec<float>{1, 1});
    Line line2(Vec<float>{0, 1}, Vec<float>{1, 0});
    Vec<float> intercect = line1.intercect(line2);
    REQUIRE(intercect == Vec<float>{0.5f, 0.5f});
  }

  {
    Line line1(Vec<float>{-2.22f, -2.53f}, Vec<float>{5.04f, 2.17f});
    Line line2(Vec<float>{-0.76f, 4.99f}, Vec<float>{6.22f, -1.85f});
    Vec<float> intercect = line1.intercect(line2);
    REQUIRE_THAT(intercect[0], Catch::Matchers::WithinAbs(3.28f, 0.01f));
  }
}

TEST_CASE("Line::orthogonal", "[Line]") {
  {
    Line line1(Vec<float>{0, 0}, Vec<float>{1, 1});
    Line line2(Vec<float>{0, 1}, Vec<float>{1, 0});
    REQUIRE(line1.orthogonal(line2));
  }

  {
    Line line1(Vec<float>{1.5f, 0.7f}, Vec<float>{0.55f, -1.3f});
    Line line2(Vec<float>{1.025f, -0.3f}, Vec<float>{0.393421f, 0});
    REQUIRE(line1.orthogonal(line2));
  }
}

TEST_CASE("Line::perpendicular_bisector", "[Line]") {
  {
    Line line(Vec<float>{0, 0}, Vec<float>{1, 1});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == Vec<float>{0.5f, 0.5f});
    REQUIRE(line.orthogonal(bisector));
  }

  {
    Line line(Vec<float>{0, 0}, Vec<float>{2.0f / 3.0f, 2.0f / std::sqrt(3.0f)});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == line.midpoint());
    REQUIRE(line.orthogonal(bisector));
  }
}

TEST_CASE("Line::midpoint", "[Line]") {
  {
    Line line(Vec<float>{0, 0}, Vec<float>{1, 1});
    REQUIRE(line.midpoint() == Vec<float>{0.5f, 0.5f});
  }

  {
    Line line(Vec<float>{0, 0}, Vec<float>{2.0f / 3.0f, 2.0f / std::sqrt(3.0f)});
    REQUIRE(line.midpoint() == Vec<float>{1.0f / 3.0f, 1.0f / std::sqrt(3.0f)});
  }
}
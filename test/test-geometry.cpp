#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

#include "Geometry/Line.h"

TEST_CASE("Line::intercect", "[Line]") {
  {
    Line line1(NewVec<2>{0, 0}, NewVec<2>{1, 1});
    Line line2(NewVec<2>{0, 1}, NewVec<2>{1, 0});
    NewVec<2> intercect = line1.intercect(line2);
    REQUIRE(intercect == NewVec<2>{0.5, 0.5});
  }

  {
    Line line1(NewVec<2>{-2.22, -2.53}, NewVec<2>{5.04, 2.17});
    Line line2(NewVec<2>{-0.76, 4.99}, NewVec<2>{6.22, -1.85});
    NewVec<2> intercect = line1.intercect(line2);
    REQUIRE_THAT(intercect[0], Catch::Matchers::WithinAbs(3.28, 0.01));
  }
}

TEST_CASE("Line::orthogonal", "[Line]") {
  {
    Line line1(NewVec<2>{0, 0}, NewVec<2>{1, 1});
    Line line2(NewVec<2>{0, 1}, NewVec<2>{1, 0});
    REQUIRE(line1.orthogonal(line2));
  }

  {
    Line line1(NewVec<2>{1.5, 0.7}, NewVec<2>{0.55, -1.3});
    Line line2(NewVec<2>{1.025, -0.3}, NewVec<2>{0.393421, 0});
    REQUIRE(line1.orthogonal(line2));
  }
}

TEST_CASE("Line::perpendicular_bisector", "[Line]") {
  {
    Line line(NewVec<2>{0, 0}, NewVec<2>{1, 1});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == NewVec<2>{0.5, 0.5});
    REQUIRE(line.orthogonal(bisector));
  }

  {
    Line line(NewVec<2>{0, 0}, NewVec<2>{2.0 / 3.0, 2.0 / std::sqrt(3.0)});
    Line bisector = line.perpendicular_bisector();
    REQUIRE(bisector.p1() == line.midpoint());
    REQUIRE(line.orthogonal(bisector));
  }
}

TEST_CASE("Line::midpoint", "[Line]") {
  {
    Line line(NewVec<2>{0, 0}, NewVec<2>{1, 1});
    REQUIRE(line.midpoint() == NewVec<2>{0.5, 0.5});
  }

  {
    Line line(NewVec<2>{0, 0}, NewVec<2>{2.0 / 3.0, 2.0 / std::sqrt(3.0)});
    REQUIRE(line.midpoint() == NewVec<2>{1.0 / 3.0, 1.0 / std::sqrt(3.0)});
  }
}

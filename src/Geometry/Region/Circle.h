#pragma once

#include "Geometry/Point.h"
#include "Geometry/Region/Region.h"

class Circle : public Region<2> {
 public:
  Circle(){};
  Circle(const Point<2>& origin, double radius) : origin_(origin), radius_(radius) {}

  Point<2> origin() const override { return origin_; }
  bool contains(const Point<2>& point) const override;
  std::vector<Point<2>> grid(size_t n) const override;
  double radius() const { return radius_; }

 private:
  Point<2> origin_;
  double radius_;
};

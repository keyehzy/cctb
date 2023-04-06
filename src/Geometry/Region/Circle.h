#pragma once

#include "Geometry/Point.h"
#include "Geometry/Region/Region.h"

class Circle : public Region<2> {
 public:
  Circle(){};
  Circle(const Point<2>& origin, double radius) : origin_(origin), radius_(radius) {}

  virtual Point<2> origin() const override { return origin_; }
  virtual bool contains(const Point<2>& point) const override;
  double radius() const { return radius_; }

 private:
  Point<2> origin_;
  double radius_;
};

#pragma once

#include "Geometry/Point.h"

class Line {
 public:
  Line() {}
  Line(const Point<2>& p1, const Point<2>& p2) : p1_(p1), p2_(p2) {}

  const Point<2>& p1() const { return p1_; }
  const Point<2>& p2() const { return p2_; }

  Point<2> p1_;
  Point<2> p2_;

  Point<2> intercect(const Line& other) const;
  Line perpendicular_bisector() const;
  bool orthogonal(const Line& other) const;
  Point<2> midpoint() const;

  void Print() const;
};

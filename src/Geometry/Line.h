#pragma once

#include "Geometry/Point.h"
#include "Geometry/Region/Circle.h"
#include "Geometry/Vector.h"

class Line {
 public:
  Line() {}
  Line(const Point<2>& p1, const Point<2>& p2) : p1_(p1), p2_(p2) {}
  Line(const Line& other) : p1_(other.p1_), p2_(other.p2_) {}
  Line(const Point<2>& p, const Vector<2>& v) : p1_(p), p2_(p.translated(v)) {}

  const Point<2>& p1() const { return p1_; }
  const Point<2>& p2() const { return p2_; }

  Point<2> intersection_with(const Line& other) const;
  Point<2> intersection_with(const Circle& circle) const;
  Line perpendicular_bisector() const;
  bool orthogonal(const Line& other) const;
  Point<2> midpoint() const;

  void Print() const;

 private:
  Point<2> p1_;
  Point<2> p2_;
};

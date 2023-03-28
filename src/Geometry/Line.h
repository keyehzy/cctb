#pragma once

#include "cctb/vec.h"

class Line {
 public:
  Line() {}
  Line(const Vec<double>& p1, const Vec<double>& p2) : p1_(p1), p2_(p2) {}

  const Vec<double>& p1() const { return p1_; }
  const Vec<double>& p2() const { return p2_; }

  Vec<double> p1_;
  Vec<double> p2_;

  Vec<double> intercect(const Line& other) const;
  Line perpendicular_bisector() const;
  bool orthogonal(const Line& other) const;
  Vec<double> midpoint() const;

  void Print() const;
};

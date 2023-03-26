#pragma once

#include "cctb/vec.h"

class Line {
 public:
  Line() {}
  Line(const Vec<float>& p1, const Vec<float>& p2) : p1_(p1), p2_(p2) {}

  const Vec<float>& p1() const { return p1_; }
  const Vec<float>& p2() const { return p2_; }

  Vec<float> p1_;
  Vec<float> p2_;

  Vec<float> intercect(const Line& other) const;
  Line perpendicular_bisector() const;
  bool orthogonal(const Line& other) const;
  Vec<float> midpoint() const;

  void Print() const;
};

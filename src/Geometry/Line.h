#pragma once

#include "cctb/vec.h"

class Line {
 public:
  Line() {}
  Line(const NewVec<2>& p1, const NewVec<2>& p2) : p1_(p1), p2_(p2) {}

  const NewVec<2>& p1() const { return p1_; }
  const NewVec<2>& p2() const { return p2_; }

  NewVec<2> p1_;
  NewVec<2> p2_;

  NewVec<2> intercect(const Line& other) const;
  Line perpendicular_bisector() const;
  bool orthogonal(const Line& other) const;
  NewVec<2> midpoint() const;

  void Print() const;
};

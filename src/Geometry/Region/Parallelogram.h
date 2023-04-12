#pragma once

#include "Geometry/Region/Region.h"

class Parallelogram : public Region<2> {
 public:
  Parallelogram(const Point<2>& origin, const Vector<2>& u, const Vector<2>& v)
      : origin_(origin), u_(u), v_(v) {}
  virtual ~Parallelogram() {}

  Point<2> origin() const override { return origin_; }
  bool contains(const Point<2>& point) const override;
  Mesh MakeMesh(int n) const override;
  Vector<2> u() const { return u_; }
  Vector<2> v() const { return v_; }

 private:
  Point<2> origin_;
  Vector<2> u_;
  Vector<2> v_;
};

#pragma once

#include <cassert>

#include "Geometry/Point.h"

template <size_t N>
class Region {
 public:
  Region() {}
  virtual ~Region() {}

  virtual Point<N> origin() const = 0;
  virtual bool contains(const Point<N>& point) const = 0;
};

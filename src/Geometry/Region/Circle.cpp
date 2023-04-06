#include "Circle.h"

#include "Geometry/Line.h"

bool Circle::contains(const Point<2>& point) const {
  return (point[0] - origin_[0]) * (point[0] - origin_[0]) +
             (point[1] - origin_[1]) * (point[1] - origin_[1]) <=
         radius_ * radius_;
}

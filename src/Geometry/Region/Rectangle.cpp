#include "Geometry/Region/Rectangle.h"

#include "Geometry/Line.h"

bool Rectangle::contains(const Point<2>& point) const {
  return point[0] >= origin_[0] && point[0] <= origin_[0] + width_ && point[1] >= origin_[1] &&
         point[1] <= origin_[1] + height_;
}

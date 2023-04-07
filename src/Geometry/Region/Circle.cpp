#include "Circle.h"

#include "Geometry/Line.h"

bool Circle::contains(const Point<2>& point) const {
  return (point[0] - origin_[0]) * (point[0] - origin_[0]) +
             (point[1] - origin_[1]) * (point[1] - origin_[1]) <=
         radius_ * radius_;
}

std::vector<Point<2>> Circle::grid(size_t n) const {
  std::vector<Point<2>> grid;
  double dx = 2 * radius_ / (n + 1);
  double dy = 2 * radius_ / (n + 1);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      Point<2> point = {origin_[0] - radius_ + dx * i, origin_[1] - radius_ + dy * j};
      if (contains(point)) {
        grid.push_back(point);
      }
    }
  }
  return grid;
}

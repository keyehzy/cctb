#include "Geometry/Region/Rectangle.h"

#include "Geometry/Line.h"

bool Rectangle::contains(const Point<2>& point) const {
  return point[0] >= origin_[0] && point[0] <= origin_[0] + width_ && point[1] >= origin_[1] &&
         point[1] <= origin_[1] + height_;
}

std::vector<Point<2>> Rectangle::grid(size_t n) const {
  std::vector<Point<2>> grid;
  int grid_size[2] = {10, 10};
  double dx = width_ / static_cast<double>(grid_size[0] - 1);
  double dy = height_ / static_cast<double>(grid_size[1] - 1);
  for (int i = 0; i < grid_size[0]; ++i) {
    for (int j = 0; j < grid_size[1]; ++j) {
      grid.push_back(Point<2>(origin_[0] + i * dx, origin_[1] + j * dy));
    }
  }
  return grid;
}

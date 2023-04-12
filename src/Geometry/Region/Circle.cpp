#include "Circle.h"

#include "Geometry/Line.h"

bool Circle::contains(const Point<2>& point) const {
  return (point[0] - origin_[0]) * (point[0] - origin_[0]) +
             (point[1] - origin_[1]) * (point[1] - origin_[1]) <=
         radius_ * radius_;
}

Mesh Circle::MakeMesh(int n) const {
  Mesh mesh;
  double angle_step = 2 * M_PI / n;
  double radius_step = radius_ / n;
  for (int i = 0; i < n; ++i) {
    double angle = i * angle_step;
    for (int j = 0; j < n; ++j) {
      double radius = j * radius_step;
      Point<2> point = {origin_[0] + radius * cos(angle), origin_[1] + radius * sin(angle)};
      mesh.Add(point);
    }
  }
  return mesh;
}

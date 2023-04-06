#include "Geometry/Line.h"

Point<2> Line::intersection_with(const Line& other) const {
  double x1 = p1_[0];
  double y1 = p1_[1];
  double x2 = p2_[0];
  double y2 = p2_[1];
  double x3 = other.p1()[0];
  double y3 = other.p1()[1];
  double x4 = other.p2()[0];
  double y4 = other.p2()[1];
  double x = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) /
             ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
  double y = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) /
             ((x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4));
  return Point<2>{x, y};
}

Line Line::perpendicular_bisector() const {
  double x1 = p1_[0];
  double y1 = p1_[1];
  double x2 = p2_[0];
  double y2 = p2_[1];
  double x = (x1 + x2) / 2;
  double y = (y1 + y2) / 2;
  return Line(Point<2>{x, y}, Point<2>{x + -y, y + x});
}

bool Line::orthogonal(const Line& other) const {
  double dx = (p2_[0] - p1_[0]);
  double dy = (p2_[1] - p1_[1]);
  double dx2 = (other.p2()[0] - other.p1()[0]);
  double dy2 = (other.p2()[1] - other.p1()[1]);
  return fp_eq(dx * dx2 + dy * dy2, 0.0);
}

Point<2> Line::midpoint() const {
  double x1 = p1_[0];
  double y1 = p1_[1];
  double x2 = p2_[0];
  double y2 = p2_[1];
  double x = (x1 + x2) / 2;
  double y = (y1 + y2) / 2;
  return Point<2>{x, y};
}

void Line::Print() const {
  std::cout << "Line: (" << p1_[0] << ", " << p1_[1] << ") -> (" << p2_[0] << ", " << p2_[1] << ")"
            << std::endl;
}

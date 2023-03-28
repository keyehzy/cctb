#include "Geometry/Line.h"

Vec<double> Line::intercect(const Line& other) const {
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
  return Vec<double>{x, y};
}

Line Line::perpendicular_bisector() const {
  double x1 = p1_[0];
  double y1 = p1_[1];
  double x2 = p2_[0];
  double y2 = p2_[1];
  double x = (x1 + x2) / 2;
  double y = (y1 + y2) / 2;
  return Line(Vec<double>{x, y}, Vec<double>{x + -y, y + x});
}

bool Line::orthogonal(const Line& other) const {
  Vec<double> p1 = p2_ - p1_;
  Vec<double> p2 = other.p2() - other.p1();
  return std::abs(p1.dot(p2)) < 1e-6;
}

Vec<double> Line::midpoint() const {
  double x1 = p1_[0];
  double y1 = p1_[1];
  double x2 = p2_[0];
  double y2 = p2_[1];
  double x = (x1 + x2) / 2;
  double y = (y1 + y2) / 2;
  return Vec<double>{x, y};
}

void Line::Print() const {
  std::cout << "Line: (" << p1_[0] << ", " << p1_[1] << ") -> (" << p2_[0] << ", " << p2_[1] << ")"
            << std::endl;
}
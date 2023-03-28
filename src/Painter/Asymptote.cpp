#include "Painter/Painter.h"

void AsymptotePainter::Prepare() const {
  out_ << "import graph;\n";
  out_ << "size(500,500);\n";
};

void AsymptotePainter::SetAxis(double xmin, double xmax, double ymin, double ymax) const {
  out_ << "draw(( " << xmin << ", " << 0 << ") -- ( " << xmax << ", " << 0 << "));\n";
  out_ << "draw(( " << 0 << ", " << ymin << ") -- ( " << 0 << ", " << ymax << "));\n";
}

void AsymptotePainter::DrawLine(double x1, double y1, double x2, double y2) const {
  out_ << "draw(( " << x1 << ", " << y1 << ") -- ( " << x2 << ", " << y2 << "));\n";
}

void AsymptotePainter::DrawArrow(double x1, double y1, double x2, double y2) const {
  out_ << "draw(( " << x1 << ", " << y1 << ") -- ( " << x2 << ", " << y2 << "), red, Arrow);\n";
}

void AsymptotePainter::DrawDottedLine(double x1, double y1, double x2, double y2) const {
  out_ << "draw(( " << x1 << ", " << y1 << ") -- ( " << x2 << ", " << y2 << "), dotted);\n";
}

void AsymptotePainter::DrawPoint(double x, double y) const {
  out_ << "fill(circle(( " << x << ", " << y << "), 0.1));\n";
}

void AsymptotePainter::DrawText(double x, double y, const std::string& text) const {
  out_ << "label(\"" << text << "\", ( " << x << ", " << y << "), fontsize(20pt));\n";
}

void AsymptotePainter::Finish() const {
  return;  // nothing to do
}
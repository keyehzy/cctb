#include "Painter/Painter.h"

void AsymptotePainter::Prepare() const {
  out_ << "import graph;\n";
  out_ << "size(500,500);\n";
};

void AsymptotePainter::DrawLine(float x1, float y1, float x2, float y2) const {
  out_ << "draw(( " << x1 << ", " << y1 << ") -- ( " << x2 << ", " << y2
       << "));\n";
}

void AsymptotePainter::DrawArrow(float x1, float y1, float x2, float y2) const {
  out_ << "draw(( " << x1 << ", " << y1 << ") -- ( " << x2 << ", " << y2
       << "), red, Arrow);\n";
}

void AsymptotePainter::DrawPoint(float x, float y) const {
  out_ << "fill(circle(( " << x << ", " << y << "), 0.1));\n";
}

void AsymptotePainter::DrawText(float x, float y,
                                const std::string& text) const {
  out_ << "label(\"" << text << "\", ( " << x << ", " << y
       << "), fontsize(20pt));\n";
}

void AsymptotePainter::Finish() const {
  return;  // nothing to do
}
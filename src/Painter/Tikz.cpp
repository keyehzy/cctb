#include "Painter/Painter.h"

void TikzPainter::Prepare() const {
  out_ << "\\documentclass[tikz]{standalone}\n";
  out_ << "\\begin{document}\n";
  out_ << "\\begin{tikzpicture}\n";
}

void TikzPainter::SetAxis(float xmin, float xmax, float ymin,
                          float ymax) const {
  out_ << "\\draw (" << xmin << "," << 0 << ") -- (" << xmax << "," << 0
       << ");\n";
  out_ << "\\draw (" << 0 << "," << ymin << ") -- (" << 0 << "," << ymax
       << ");\n";
}

void TikzPainter::DrawLine(float x1, float y1, float x2, float y2) const {
  out_ << "\\draw (" << x1 << "," << y1 << ") -- (" << x2 << "," << y2
       << ");\n";
}

void TikzPainter::DrawArrow(float x1, float y1, float x2, float y2) const {
  out_ << "\\draw[->, red] (" << x1 << "," << y1 << ") -- (" << x2 << "," << y2
       << ");\n";
}

void TikzPainter::DrawDottedLine(float x1, float y1, float x2, float y2) const {
  out_ << "\\draw[dotted] (" << x1 << "," << y1 << ") -- (" << x2 << "," << y2
       << ");\n";
}

void TikzPainter::DrawPoint(float x, float y) const {
  out_ << "\\filldraw (" << x << "," << y << ") circle (3pt);\n";
}

void TikzPainter::DrawText(float x, float y, const std::string& text) const {
  out_ << "\\node at (" << x << "," << y << ") [font=\\small] {" << text
       << "};\n";
}

void TikzPainter::Finish() const {
  out_ << "\\end{tikzpicture}\n";
  out_ << "\\end{document}\n";
}
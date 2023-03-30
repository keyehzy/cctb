#include "cctb/lattice.h"

#include <ostream>

#include "Geometry/Line.h"

void OneDimensionalLattice::Plot(PainterBackend backend, std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);

  plotter->Prepare();

  for (const Point<1> &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      Point<1> p = site.translated(i * m_a1);
      plotter->DrawPoint(p[0], 0.0);
    }
  }

  for (const Edge &edge : Edges()) {
    for (int i = -1; i <= 1; i++) {
      Point<1> from_position = SiteAt(edge.from_index).translated(i * m_a1);
      Point<1> to_position =
          SiteAt(edge.to_index).translated(i * m_a1 + edge.relative_index[0] * m_a1);
      plotter->DrawLine(from_position[0], 0.0, to_position[0], 0.0);
    }
  }

  plotter->Finish();
}

Matrix<std::complex<double>> OneDimensionalLattice::HoppingMatrix(Vector<1> k) const {
  Matrix<std::complex<double>> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    Point<1> from_position = SiteAt(edge.from_index);
    Point<1> to_position = SiteAt(edge.to_index).translated(edge.relative_index[0] * m_a1);
    Vector<1> r(from_position, to_position);
    double dot = k.dot(r);
    std::complex<double> phase = std::complex<double>(0, dot);
    hamiltonian(edge.from_index, edge.to_index) += edge.weight * std::exp(phase);
    hamiltonian(edge.to_index, edge.from_index) += edge.weight * std::exp(-phase);
  }
  return hamiltonian;
}

void TwoDimensionalLattice::Plot(PainterBackend backend, std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);
  plotter->Prepare();

  for (const Point<2> &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Vector<2> offset = m_a1 * i + m_a2 * j;
        Point<2> p = site.translated(offset);
        plotter->DrawPoint(p[0], p[1]);
      }
    }
  }

  for (const Edge &edge : Edges()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Point<2> from_position = SiteAt(edge.from_index).translated(m_a1 * i + m_a2 * j);
        Point<2> to_position = SiteAt(edge.to_index)
                                   .translated(m_a1 * i + m_a2 * j + m_a1 * edge.relative_index[0] +
                                               m_a2 * edge.relative_index[1]);
        plotter->DrawLine(from_position[0], from_position[1], to_position[0], to_position[1]);
      }
    }
  }

  plotter->DrawText(m_a1[0] / 2, m_a1[1] / 2 + 0.25, "$a_1$");
  plotter->DrawArrow(0, 0, m_a1[0], m_a1[1]);
  plotter->DrawText(m_a2[0] / 2, m_a2[1] / 2 - 0.25, "$a_2$");
  plotter->DrawArrow(0, 0, m_a2[0], m_a2[1]);

  plotter->Finish();
}

void TwoDimensionalLattice::PlotBrillouinZone(PainterBackend backend, std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);
  plotter->Prepare();

  // Draw Axes
  double max_x = std::max(m_b1[0], m_b2[0]);
  double max_y = std::max(m_b1[1], m_b2[1]);
  double max_both = std::max(max_x, max_y);
  plotter->SetAxis(-1.2 * max_both, 1.2 * max_both, -1.2 * max_both, 1.2 * max_both);
  plotter->DrawText(1.2 * max_both, 0.25, "$k_x$");
  plotter->DrawText(0.25, 1.2 * max_both, "$k_y$");

  // Calculate all perpendicular bisectors
  Line b1_line(Point<2>(0, 0), m_b1);
  Line b2_line(Point<2>(0, 0), m_b2);
  Line b3_line(Point<2>(0, 0), m_b1 + m_b2);
  Line b1_line_mirrored(Point<2>(0, 0), m_b1 * -1.0);
  Line b2_line_mirrored(Point<2>(0, 0), m_b2 * -1.0);
  Line b3_line_mirrored(Point<2>(0, 0), (m_b1 + m_b2) * -1.0);
  Line b1_perp = b1_line.perpendicular_bisector();
  Line b2_perp = b2_line.perpendicular_bisector();
  Line b3_perp = b3_line.perpendicular_bisector();
  Line b1_perp_mirrored = b1_line_mirrored.perpendicular_bisector();
  Line b2_perp_mirrored = b2_line_mirrored.perpendicular_bisector();
  Line b3_perp_mirrored = b3_line_mirrored.perpendicular_bisector();

  // Intersect bisectors
  Point<2> k1 = b1_perp.intercect(b2_perp_mirrored);
  Point<2> k2 = b2_perp.intercect(b1_perp_mirrored);
  Point<2> k3 = b1_perp.intercect(b3_perp);
  Point<2> k4 = b2_perp.intercect(b3_perp);
  Point<2> k5 = b1_perp_mirrored.intercect(b3_perp_mirrored);
  Point<2> k6 = b2_perp_mirrored.intercect(b3_perp_mirrored);

  // Get midpoints
  Point<2> b1_mid = b1_line.midpoint();
  Point<2> b2_mid = b2_line.midpoint();

  // Draw lines
  plotter->DrawLine(b1_mid[0], b1_mid[1], k1[0], k1[1]);
  plotter->DrawLine(-b2_mid[0], -b2_mid[1], k1[0], k1[1]);

  plotter->DrawLine(b2_mid[0], b2_mid[1], k2[0], k2[1]);
  plotter->DrawLine(-b1_mid[0], -b1_mid[1], k2[0], k2[1]);

  plotter->DrawLine(b1_mid[0], b1_mid[1], k3[0], k3[1]);
  plotter->DrawLine(b2_mid[0], b2_mid[1], k4[0], k4[1]);
  plotter->DrawLine(-b1_mid[0], -b1_mid[1], k5[0], k5[1]);
  plotter->DrawLine(-b2_mid[0], -b2_mid[1], k6[0], k6[1]);
  plotter->DrawLine(k3[0], k3[1], k4[0], k4[1]);
  plotter->DrawLine(k5[0], k5[1], k6[0], k6[1]);

  // Draw Brillouin Zone
  plotter->DrawArrow(0, 0, m_b1[0], m_b1[1]);
  plotter->DrawArrow(0, 0, m_b2[0], m_b2[1]);
  plotter->DrawDottedLine(m_b1[0], m_b1[1], m_b1[0] + m_b2[0], m_b1[1] + m_b2[1]);
  plotter->DrawDottedLine(m_b2[0], m_b2[1], m_b1[0] + m_b2[0], m_b1[1] + m_b2[1]);

  // Draw labels
  plotter->DrawText(m_b1[0] + 0.35, m_b1[1] + 0.35, "$b_1$");
  plotter->DrawText(m_b2[0] + 0.35, m_b2[1] - 0.35, "$b_2$");

  plotter->DrawText(-.2, .2, "$\\Gamma$");

  plotter->DrawPoint(k3[0], k3[1]);
  plotter->Finish();
}

Matrix<std::complex<double>> TwoDimensionalLattice::HoppingMatrix(Vector<2> k) const {
  Matrix<std::complex<double>> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    Point<2> from_position = SiteAt(edge.from_index);
    Point<2> to_position =
        SiteAt(edge.to_index)
            .translated(m_a1 * edge.relative_index[0] + m_a2 * edge.relative_index[1]);
    Vector<2> r(from_position, to_position);
    double dot = k.dot(r);
    std::complex<double> phase = std::complex<double>(0, dot);
    hamiltonian(edge.from_index, edge.to_index) += edge.weight * std::exp(phase);
    hamiltonian(edge.to_index, edge.from_index) += edge.weight * std::exp(-phase);
  }
  return hamiltonian;
}

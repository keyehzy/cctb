#include "Lattice/Lattice.h"

#include <ostream>
#include <unordered_set>

#include "Geometry/Line.h"
#include "Geometry/Region/Parallelogram.h"

Eigen::MatrixXcd OneDimensionalLattice::HoppingMatrix(Vector<1> k) const {
  Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(size(), size());
  for (int src = 0; src < size(); src++) {
    for (const auto &edge : site(src).edges) {
      Point<1> to_position = site(edge.dst).position.translated(edge.offset[0] * m_a1);
      Vector<1> r(site(src).position, to_position);
      double dot = k.dot(r);
      std::complex<double> phase = std::complex<double>(0, dot);
      H(src, edge.dst) += edge.weight * std::exp(phase);
      H(edge.dst, src) += edge.weight * std::exp(-phase);
    }
  }
  return H;
}

void OneDimensionalLattice::Plot(PainterBackend backend, std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);
  plotter->Prepare();
  for (const auto &site : sites()) {
    for (int i = -1; i <= 1; i++) {
      Point<1> p = site.position.translated(i * m_a1);
      plotter->DrawPoint(p[0], 0.0);
      for (const auto &edge : site.edges) {
        Point<1> to_position = this->site(edge.dst).position.translated(edge.offset[0] * m_a1);
        plotter->DrawLine(p[0], 0.0, to_position[0], 0.0);
      }
    }
  }
  plotter->Finish();
}

void OneDimensionalLattice::PlotBandStructure(std::ostream &out) const {
  out << "\\documentclass[border=10pt]{standalone}" << '\n';
  out << "\\usepackage{pgfplots}" << '\n';
  out << "\\usepackage{pgfplotstable}" << '\n';
  out << "\\pgfplotstableread{" << '\n';

  int n = 25;
  for (int i = 0; i < n; i++) {
    double k = m_b1[0] * i / static_cast<double>(n - 1);
    auto Hk = HoppingMatrix(Vector<1>(k));
    Eigen::Vector<double, 1> w(size());
    Eigen::MatrixXcd v(size(), size());
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(Hk);
    w = solver.eigenvalues();
    v = solver.eigenvectors();

    out << k;
    for (int k = 0; k < size(); k++) {
      out << " " << w[k];
    }
    out << '\n';
  }

  out << "}{\\data}" << '\n';
  out << "\\pgfplotsset{width=7cm,compat=1.18}" << '\n';
  out << "\\begin{document}" << '\n';
  out << "\\begin{tikzpicture}" << '\n';
  out << "\\begin{axis}[" << '\n';
  out << "title={}," << '\n';
  out << "xlabel={$k$}," << '\n';
  out << "ylabel={$E(k)$}," << '\n';
  out << "grid=both," << '\n';
  out << "minor tick num=1," << '\n';
  out << "major grid style={gray!50, very thin}," << '\n';
  out << "minor grid style={gray!50, ultra thin}," << '\n';
  out << "colormap name=colormap/jet," << '\n';
  out << "cycle list={[of colormap]}," << '\n';

  out << "legend style={at={(1.3,0.5)},anchor=east}," << '\n';
  out << "legend entries={";
  for (int i = 0; i < size(); i++) {
    out << "$E_" << i << "$";
    if (i != size() - 1) {
      out << ",";
    }
  }
  out << "}," << '\n';

  out << "]" << '\n';

  for (int i = 0; i < size(); i++) {
    out << "\\addplot+[mark=none, style=ultra thick] table[x index=0,y index=" << i + 1
        << "] {\\data};" << '\n';
  }

  out << "\\end{axis}" << '\n';
  out << "\\end{tikzpicture}" << '\n';
  out << "\\end{document}" << '\n';
}

void TwoDimensionalLattice::Plot(PainterBackend backend, std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);
  plotter->Prepare();
  for (const auto &site : sites()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Vector<2> offset = m_a1 * i + m_a2 * j;
        Point<2> p = site.position.translated(offset);
        plotter->DrawPoint(p[0], p[1]);
      }
    }
  }
  for (int src = 0; src < size(); src++) {
    for (const auto &edge : site(src).edges) {
      for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
          Point<2> from_position = site(src).position.translated(m_a1 * i + m_a2 * j);
          Point<2> to_position = site(edge.dst).position.translated(
              m_a1 * i + m_a2 * j + m_a1 * edge.offset[0] + m_a2 * edge.offset[1]);
          plotter->DrawLine(from_position[0], from_position[1], to_position[0], to_position[1]);
        }
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
  Point<2> k1 = b1_perp.intersection_with(b2_perp_mirrored);
  Point<2> k2 = b2_perp.intersection_with(b1_perp_mirrored);
  Point<2> k3 = b1_perp.intersection_with(b3_perp);
  Point<2> k4 = b2_perp.intersection_with(b3_perp);
  Point<2> k5 = b1_perp_mirrored.intersection_with(b3_perp_mirrored);
  Point<2> k6 = b2_perp_mirrored.intersection_with(b3_perp_mirrored);

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

Eigen::MatrixXcd TwoDimensionalLattice::HoppingMatrix(Vector<2> k) const {
  Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(size(), size());
  for (int src = 0; src < size(); src++) {
    for (const auto &edge : site(src).edges) {
      Point<2> to_position =
          site(edge.dst).position.translated(m_a1 * edge.offset[0] + m_a2 * edge.offset[1]);
      Vector<2> r(site(src).position, to_position);
      double dot = k.dot(r);
      std::complex<double> phase = std::complex<double>(0, dot);
      H(src, edge.dst) += edge.weight * std::exp(phase);
      H(edge.dst, src) += edge.weight * std::exp(-phase);
    }
  }
  return H;
}

void TwoDimensionalLattice::PlotBandStructure(std::ostream &out) const {
  out << "\\documentclass[border=10pt]{standalone}" << '\n';
  out << "\\usepackage{pgfplots}" << '\n';
  out << "\\usepackage{pgfplotstable}" << '\n';
  out << "\\pgfplotstableread{" << '\n';

  int n = 25;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double kx = (m_b1[0] * i + m_b2[0] * j) / static_cast<double>(n - 1);
      double ky = (m_b1[1] * i + m_b2[1] * j) / static_cast<double>(n - 1);
      auto Hk = HoppingMatrix(Vector<2>(kx, ky));
      Eigen::VectorXd w(size());
      Eigen::MatrixXcd v(size(), size());
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(Hk);
      w = solver.eigenvalues();
      v = solver.eigenvectors();

      out << kx << " " << ky;
      for (int k = 0; k < size(); k++) {
        out << " " << w[k];
      }
      out << '\n';
    }
    out << '\n';
  }

  out << "}{\\data}" << '\n';
  out << "\\pgfplotsset{width=7cm,compat=1.18}" << '\n';
  out << "\\begin{document}" << '\n';
  out << "\\begin{tikzpicture}" << '\n';
  out << "\\begin{axis}[" << '\n';
  out << "title={}," << '\n';
  out << "xlabel={$k_x$}," << '\n';
  out << "ylabel={$k_y$}," << '\n';
  out << "zlabel={$E(k)$}," << '\n';
  out << "grid=both," << '\n';
  out << "minor tick num=1," << '\n';
  out << "major grid style={gray!50, very thin}," << '\n';
  out << "minor grid style={gray!50, ultra thin}," << '\n';
  out << "view={25}{25}," << '\n';
  out << "colormap name=viridis," << '\n';
  out << "cycle list={[of colormap]}," << '\n';

  out << "legend style={at={(1.3,0.5)},anchor=east}," << '\n';
  // out << "legend entries={";
  // for (int i = 0; i < size(); i++) {
  //   out << "$E_" << i << "$";
  //   if (i != size() - 1) {
  //     out << ",";
  //   }
  // }
  // out << "}," << '\n';

  out << "]" << '\n';

  out << "\\pgfplotsinvokeforeach{0,...," << size() - 1
      << "}{\\pgfplotscolormapdefinemappedcolor{\\the\\numexpr(#1)*1000/" << (size() - 1) << "}"
      << '\n';
  out << "\\colorlet{leg#1}{mapped color}" << '\n';
  out << "\\addlegendimage{area legend,color=leg#1,fill}" << '\n';
  out << "\\addlegendentry{$E_{#1}$}}" << '\n';

  for (int i = 0; i < size(); i++) {
    out << "\\addplot3[surf, mesh/ordering=y varies, mesh/rows=25, opacity=0.8, point meta=" << i
        << "] table[x index=0,y index=1,z index=" << i + 2 << "] {\\data};" << '\n';
  }

  out << "\\end{axis}" << '\n';
  out << "\\end{tikzpicture}" << '\n';
  out << "\\end{document}" << '\n';
}

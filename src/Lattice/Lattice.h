#pragma once

#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

#include "Painter/Painter.h"

template <typename NumberType, size_t Dimension, size_t Orbitals>
class Node {
 public:
  using ComplexType = std::complex<NumberType>;
  using VectorType = Eigen::Vector<NumberType, Dimension>;
  using MatrixType = Eigen::Matrix<ComplexType, Orbitals, Orbitals>;

  explicit Node(size_t index, const VectorType &position, const MatrixType &onsite)
      : m_index(index), m_position(position), m_onsite(onsite) {}

  size_t index() const { return m_index; }

  const VectorType &position() const { return m_position; }

  const MatrixType &onsite() const { return m_onsite; }

 private:
  size_t m_index;
  VectorType m_position;
  MatrixType m_onsite;
};

template <typename NumberType, size_t Dimension, size_t Orbitals>
class Edge {
 public:
  using ComplexType = std::complex<NumberType>;
  using MatrixType = Eigen::Matrix<ComplexType, Orbitals, Orbitals>;

  explicit Edge(size_t src, size_t dst, const std::array<int, Dimension> &basis_index,
                const MatrixType &hopping)
      : m_src(src), m_dst(dst), m_basis_indices(basis_index), m_hopping(hopping) {}

  size_t src() const { return m_src; }

  size_t dst() const { return m_dst; }

  int basis_index(size_t i) const { return m_basis_indices[i]; }

  const MatrixType &hopping() const { return m_hopping; }

 private:
  size_t m_src;
  size_t m_dst;
  std::array<int, Dimension> m_basis_indices;
  MatrixType m_hopping;
};

template <typename NumberType, size_t Dimension, size_t Orbitals>
class Graph {
 public:
  using NodeType = Node<NumberType, Dimension, Orbitals>;
  using EdgeType = Edge<NumberType, Dimension, Orbitals>;

  void add_node(NodeType node) { m_nodes.push_back(node); }

  void add_edge(EdgeType edge) { m_edges.push_back(edge); }

  size_t size() const { return m_nodes.size(); }

  Eigen::MatrixXi adj_matrix() const {
    Eigen::MatrixXi Adj = Eigen::MatrixXi::Zero(m_nodes.size(), m_nodes.size());
    for (const EdgeType &edge : m_edges) {
      Adj(edge.src(), edge.dst()) = 1;
      Adj(edge.dst(), edge.src()) = 1;
    }
    return Adj;
  }

  const std::vector<NodeType> &nodes() const { return m_nodes; }

  const std::vector<EdgeType> &edges() const { return m_edges; }

  const NodeType &node(size_t i) const { return m_nodes[i]; }

  const EdgeType &edge(size_t i) const { return m_edges[i]; }

 private:
  std::vector<NodeType> m_nodes;
  std::vector<EdgeType> m_edges;
};

template <size_t Dimension, size_t Orbitals>
class Lattice {
 public:
  using NumberType = double;
  using ComplexType = std::complex<NumberType>;
  using NodeType = Node<double, Dimension, Orbitals>;
  using EdgeType = Edge<double, Dimension, Orbitals>;
  using GraphType = Graph<double, Dimension, Orbitals>;
  using VectorType = Eigen::Vector<NumberType, Dimension>;

  explicit Lattice(const std::array<VectorType, Dimension> &basis) : m_basis(basis) {
    // https://en.wikipedia.org/wiki/Reciprocal_lattice#Generalization_of_a_dual_lattice
    Eigen::Matrix<NumberType, Dimension, Dimension> basis_matrix;

    for (size_t i = 0; i < Dimension; ++i) {
      basis_matrix.col(i) = m_basis[i];
    }

    Eigen::Matrix<NumberType, Dimension, Dimension> reciprocal_basis_matrix =
        2 * M_PI * basis_matrix.inverse().transpose();

    for (size_t i = 0; i < Dimension; ++i) {
      m_reciprocal_basis[i] = reciprocal_basis_matrix.col(i);
    }
  }

  Eigen::MatrixXcd hopping_matrix(Eigen::Vector<NumberType, Dimension> k) const {
    Eigen::MatrixXcd H = Eigen::MatrixXcd::Zero(size(), size());

    for (const NodeType &node : m_graph.nodes()) {
      H.block<Orbitals, Orbitals>(node.index() * Orbitals, node.index() * Orbitals) +=
          node.onsite();
    }

    for (const EdgeType &edge : m_graph.edges()) {
      VectorType basis_displacement = VectorType::Zero();

      for (size_t i = 0; i < Dimension; ++i) {
        basis_displacement += edge.basis_index(i) * m_basis[i];
      }

      VectorType src_pos = m_graph.node(edge.src()).position();
      VectorType dst_pos = m_graph.node(edge.dst()).position() + basis_displacement;
      VectorType r = dst_pos - src_pos;
      double k_dot_r = k.dot(r);
      H.block<Orbitals, Orbitals>(edge.src() * Orbitals, edge.dst() * Orbitals) +=
          edge.hopping() * std::exp(ComplexType(0, k_dot_r));
      H.block<Orbitals, Orbitals>(edge.dst() * Orbitals, edge.src() * Orbitals) +=
          edge.hopping().conjugate() * std::exp(ComplexType(0, -k_dot_r));
    }
    return H;
  }

  void add_site(size_t index, const VectorType &position, const Eigen::MatrixXcd &onsite) {
    m_graph.add_node(NodeType(index, position, onsite));
  }

  void add_edge(size_t src, size_t dst, const std::array<int, Dimension> &basis_index,
                const Eigen::MatrixXcd &hopping) {
    m_graph.add_edge(EdgeType(src, dst, basis_index, hopping));
  }

  size_t size() const { return m_graph.size() * Orbitals; }

  size_t orbitals() const { return Orbitals; }

  const GraphType &graph() const { return m_graph; }

  GraphType &graph() { return m_graph; }

  const std::array<VectorType, Dimension> &basis() const { return m_basis; }

  void statistics() {
    std::cout << "Lattice statistics:" << std::endl;
    std::cout << "  Sites: " << m_graph.size() << std::endl;
    std::cout << "  Dimension: " << Dimension << std::endl;
    std::cout << "  Orbitals: " << Orbitals << std::endl;
    std::cout << "  Total size: " << size() << std::endl;

    std::cout << "  Basis vectors:" << std::endl;
    for (size_t i = 0; i < Dimension; ++i) {
      std::cout << "    " << i << ": " << m_basis[i].transpose() << std::endl;
    }

    std::cout << "  Reciprocal basis vectors:" << std::endl;
    for (size_t i = 0; i < Dimension; ++i) {
      std::cout << "    " << i << ": " << m_reciprocal_basis[i].transpose() << std::endl;
    }

    std::cout << "  Sites: " << m_graph.size() << std::endl;

    for (size_t i = 0; i < m_graph.size(); ++i) {
      std::cout << "    " << i << ": " << m_graph.node(i).position().transpose() << std::endl;
    }

    std::cout << "  Hoppings: " << m_graph.edges().size() << std::endl;

    for (size_t i = 0; i < m_graph.edges().size(); ++i) {
      std::cout << "    " << i << ": " << m_graph.edge(i).src() << " -> " << m_graph.edge(i).dst()
                << std::endl;
    }

    std::cout << "  Adjacency matrix:" << std::endl;
    std::cout << m_graph.adj_matrix() << std::endl;

    std::cout << "  Hoping matrix at k = (0.5, 0.8):" << std::endl;
    std::cout << hopping_matrix(VectorType(0.5, 0.8)) << std::endl;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver(hopping_matrix(VectorType(0.5, 0.8)));
    std::cout << "  Eigenvalues at k = (0.5, 0.8):" << std::endl;

    for (size_t i = 0; i < solver.eigenvalues().size(); ++i) {
      std::cout << "    " << i << ": " << solver.eigenvalues()[i] << std::endl;
    }

    std::cout << "  Eigenvectors at k = (0.5, 0.8):" << std::endl;

    for (size_t i = 0; i < solver.eigenvalues().size(); ++i) {
      std::cout << "    " << i << ": " << solver.eigenvectors().col(i).transpose() << std::endl;
    }
  }

  // virtual void Plot(PainterBackend, std::ostream &) const;
  // virtual void PlotBandStructure(std::ostream &) const;
  // virtual void PlotBrillouinZone(std::ostream &) const;

 private:
  GraphType m_graph;
  std::array<VectorType, Dimension> m_basis;
  std::array<VectorType, Dimension> m_reciprocal_basis;
};

template <size_t Orbitals>
class OneDimensionalLattice : public Lattice<1, Orbitals> {
 public:
  using NumberType = double;
  using VectorType = Eigen::Vector<NumberType, 1>;

  explicit OneDimensionalLattice(VectorType basis) : Lattice<1, Orbitals>({basis}) {}

  // void Plot(PainterBackend, std::ostream &) const override;
  // void PlotBandStructure(std::ostream &) const override;
  // void PlotBrillouinZone(std::ostream &) const;
};

template <size_t Orbitals>
class TwoDimensionalLattice : public Lattice<2, Orbitals> {
 public:
  using NumberType = double;
  using VectorType = Eigen::Vector<NumberType, 2>;

  explicit TwoDimensionalLattice(VectorType basis1, VectorType basis2)
      : Lattice<2, Orbitals>({basis1, basis2}) {}

  // void Plot(PainterBackend, std::ostream &) const override;
  // void PlotBandStructure(std::ostream &) const override;
  // void PlotBrillouinZone(PainterBackend, std::ostream &) const;
};

/*
template <size_t Orbitals>
void OneDimensionalLattice<Orbitals>::Plot(PainterBackend backend, std::ostream &out) const {
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
*/

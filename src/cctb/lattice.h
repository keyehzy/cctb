#pragma once

#include <complex>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

#include "Geometry/Point.h"
#include "Geometry/Vector.h"
#include "Painter/Painter.h"
#include "cctb/matrix.h"

template <std::size_t D>
struct Edge {
  int dst;
  std::array<int, D> offset;
  double weight;
};

template <std::size_t D>
struct GraphNode {
  Point<D> position;
  std::vector<Edge<D>> edges;
  int size() const { return edges.size(); }
  const Edge<D> &edge(int i) const { return edges[i]; }
  Edge<D> &edge(int i) { return edges[i]; }
  GraphNode(Point<D> p) : position(p) {}
};

template <std::size_t D>
class Lattice {
 public:
  Lattice() {}
  virtual ~Lattice() {}

  void add_site(GraphNode<D> node) { m_nodes.push_back(node); }

  void add_edge(int src, int dst, std::array<int, D> offset, double weight) {
    m_nodes[src].edges.push_back({dst, offset, weight});
  }

  const std::vector<GraphNode<D>> &sites() const { return m_nodes; }

  const GraphNode<D> &site(int i) const { return m_nodes[i]; }

  GraphNode<D> &site(int i) { return m_nodes[i]; }

  int size() const { return m_nodes.size(); }

  virtual void Plot(PainterBackend, std::ostream &) const = 0;

  Matrix<int> AdjMatrix() const;

  virtual Matrix<std::complex<double>> HoppingMatrix(Vector<D> k) const = 0;

 protected:
  std::vector<GraphNode<D>> m_nodes;
};

template <std::size_t D>
Matrix<int> Lattice<D>::AdjMatrix() const {
  Matrix<int> A(m_nodes.size(), m_nodes.size());
  for (int i = 0; i < m_nodes.size(); ++i) {
    for (auto edge : m_nodes[i].edges) {
      A(i, edge.dst) = 1;
      A(edge.dst, i) = 1;
    }
  }
  return A;
}

class OneDimensionalLattice : public Lattice<1> {
 public:
  OneDimensionalLattice(Vector<1> a1) : m_a1(a1) {}

  void Plot(PainterBackend, std::ostream &) const override;

  Matrix<std::complex<double>> HoppingMatrix(Vector<1> k) const override;

 protected:
  Vector<1> m_a1;
};

class TwoDimensionalLattice : public Lattice<2> {
 public:
  TwoDimensionalLattice(const Vector<2> &a1, const Vector<2> &a2) : m_a1(a1), m_a2(a2) {
    // https://physics.stackexchange.com/questions/340860/reciprocal-lattice-in-2d
    double det = a1[0] * a2[1] - a1[1] * a2[0];
    double pref = 2.0 * M_PI / det;
    m_b1 = pref * Vector<2>(a2[1], -a2[0]);
    m_b2 = pref * Vector<2>(-a1[1], a1[0]);
  }

  Vector<2> a1() const { return m_a1; }
  Vector<2> a2() const { return m_a2; }
  Vector<2> b1() const { return m_b1; }
  Vector<2> b2() const { return m_b2; }

  void Plot(PainterBackend, std::ostream &) const override;

  void PlotBrillouinZone(PainterBackend, std::ostream &) const;

  Matrix<std::complex<double>> HoppingMatrix(Vector<2> k) const override;

 protected:
  Vector<2> m_a1;
  Vector<2> m_a2;
  Vector<2> m_b1;
  Vector<2> m_b2;
};

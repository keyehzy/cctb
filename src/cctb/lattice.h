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

struct Edge {
  std::vector<int> relative_index;
  int from_index;
  int to_index;
  double weight;
  Edge(std::vector<int> idx, int from, int to, double w = 1.0)
      : relative_index(idx), from_index(from), to_index(to), weight(w) {}
};

template <std::size_t D>
class Lattice {
 public:
  Lattice() {}
  virtual ~Lattice() {}

  void AddSite(Point<D> site) { m_sites.push_back(site); }
  void AddEdge(Edge edge) { m_edges.push_back(edge); }

  std::vector<Point<D>> Sites() const { return m_sites; }
  Point<D> SiteAt(int index) const { return m_sites[index]; }

  std::vector<Edge> Edges() const { return m_edges; }
  int Size() const { return m_sites.size(); }

  virtual void Plot(PainterBackend, std::ostream &) const = 0;

  Matrix<int> AdjMatrix() const;

  virtual Matrix<std::complex<double>> HoppingMatrix(Vector<D> k) const = 0;

 protected:
  std::vector<Point<D>> m_sites;
  std::vector<Edge> m_edges;
};

template <std::size_t D>
Matrix<int> Lattice<D>::AdjMatrix() const {
  Matrix<int> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    hamiltonian(edge.from_index, edge.to_index) = 1;
    hamiltonian(edge.to_index, edge.from_index) = 1;
  }
  return hamiltonian;
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

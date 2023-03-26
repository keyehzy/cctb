#pragma once

#include <cfloat>
#include <complex>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

#include "Painter/Painter.h"
#include "cctb/matrix.h"
#include "cctb/vec.h"

struct Site {
  Vec<float> position;
  Site(Vec<float> pos) : position(pos) {}
};

struct Edge {
  std::vector<int> relative_index;
  int from_index;
  int to_index;
  float weight;
  Edge(std::vector<int> idx, int from, int to, float w = 1.0f)
      : relative_index(idx), from_index(from), to_index(to), weight(w) {}
};

class Lattice {
 public:
  Lattice() {}

  void AddSite(Site site) { m_sites.push_back(site); }
  void AddEdge(Edge edge) { m_edges.push_back(edge); }

  std::vector<Site> Sites() const { return m_sites; }
  Site SiteAt(int index) const { return m_sites[index]; }

  std::vector<Edge> Edges() const { return m_edges; }
  int Size() const { return m_sites.size(); }

  virtual void Plot(PainterBackend, std::ostream &) const = 0;

  Matrix<int> AdjMatrix() const;

  virtual Matrix<std::complex<float>> HoppingMatrix(Vec<float> k) const = 0;

 protected:
  std::vector<Site> m_sites;
  std::vector<Edge> m_edges;
};

class OneDimensionalLattice : public Lattice {
 public:
  OneDimensionalLattice(Vec<float> a1) : m_a1(a1) {}

  void Plot(PainterBackend, std::ostream &) const override;

  Matrix<std::complex<float>> HoppingMatrix(Vec<float> k) const override;

 protected:
  Vec<float> m_a1;
};

struct EdgeData {
  float from_position_x;
  float from_position_y;
  float to_position_x;
  float to_position_y;
  float maxflt;
};

class TwoDimensionalLattice : public Lattice {
 public:
  TwoDimensionalLattice(Vec<float> a1, Vec<float> a2) : m_a1(a1), m_a2(a2) {
    // https://physics.stackexchange.com/questions/340860/reciprocal-lattice-in-2d
    float det = a1[0] * a2[1] - a1[1] * a2[0];
    float pref = 2 * M_PI / det;
    m_b1 = Vec<float>(a2[1], -a2[0]) * pref;
    m_b2 = Vec<float>(-a1[1], a1[0]) * pref;
  }

  Vec<float> a1() const { return m_a1; }
  Vec<float> a2() const { return m_a2; }
  Vec<float> b1() const { return m_b1; }
  Vec<float> b2() const { return m_b2; }

  void Plot(PainterBackend, std::ostream &) const override;

  void PlotBrillouinZone(PainterBackend, std::ostream &) const;

  Matrix<std::complex<float>> HoppingMatrix(Vec<float> k) const override;

 protected:
  Vec<float> m_a1;
  Vec<float> m_a2;
  Vec<float> m_b1;
  Vec<float> m_b2;
};

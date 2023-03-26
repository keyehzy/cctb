#pragma once

#include <cfloat>
#include <complex>
#include <fstream>
#include <iostream>
#include <vector>

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

  virtual void Plot() const = 0;

  Matrix<int> AdjMatrix() const;

  virtual Matrix<std::complex<float>> HoppingMatrix(Vec<float> k) const = 0;

 protected:
  std::vector<Site> m_sites;
  std::vector<Edge> m_edges;
};

class OneDimensionalLattice : public Lattice {
 public:
  OneDimensionalLattice(Vec<float> a1) : m_a1(a1) {}

  void Plot() const override;

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
  TwoDimensionalLattice(Vec<float> a1, Vec<float> a2) : m_a1(a1), m_a2(a2) {}

  void Plot() const override;

  Matrix<std::complex<float>> HoppingMatrix(Vec<float> k) const override;

 protected:
  Vec<float> m_a1;
  Vec<float> m_a2;
};

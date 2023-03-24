#pragma once

#include <vector>

#include "cctb/vec.h"

struct Edge {
  std::vector<int> rel_idx;
  int from_site;
  int to_site;
  Edge(std::vector<int> idx, int from, int to)
      : rel_idx(idx), from_site(from), to_site(to) {}
};

class Lattice {
 public:
  Lattice() {}

  void AddSite(Vec2f site) { m_sites.push_back(site); }

  void AddEdge(Edge edge) { m_edges.push_back(edge); }

  std::vector<Vec2f> Sites() { return m_sites; }

  std::vector<Edge> Edges() { return m_edges; }

 private:
  std::vector<Vec2f> m_sites;
  std::vector<Edge> m_edges;
};

class OneDimensionalLattice : public Lattice {
 public:
  OneDimensionalLattice(Vec2f a1, Vec2f a2) : m_a1(a1), m_a2(a2) {}

 private:
  Vec2f m_a1;
  Vec2f m_a2;
};

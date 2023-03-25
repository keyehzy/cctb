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
  Edge(std::vector<int> idx, int from, int to)
      : relative_index(idx), from_index(from), to_index(to) {}
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

  Matrix<int> AdjMatrix() const {
    Matrix<int> hamiltonian(Size(), Size());
    for (auto edge : Edges()) {
      hamiltonian(edge.from_index, edge.to_index) = 1;
      hamiltonian(edge.to_index, edge.from_index) = 1;
    }
    return hamiltonian;
  }

  virtual Matrix<std::complex<float>> HoppingMatrix(float k) const = 0;

 protected:
  std::vector<Site> m_sites;
  std::vector<Edge> m_edges;
};

class OneDimensionalLattice : public Lattice {
 public:
  OneDimensionalLattice(Vec<float> a1) : m_a1(a1) {}

  void Plot() const override {
    for (const Site &site : Sites()) {
      std::printf("%f %f\n", site.position[0], 0.0);
    }

    for (auto edge : Edges()) {
      float r = edge.relative_index[0] * m_a1[0];
      std::printf("%f %f\n", SiteAt(edge.from_index).position[0] + r, 0.0);
      std::printf("%f %f\n", SiteAt(edge.from_index).position[0] - r, 0.0);
    }
  }

  Matrix<std::complex<float>> HoppingMatrix(float k) const override {
    Matrix<std::complex<float>> hamiltonian(Size(), Size());
    for (auto edge : Edges()) {
      float r = edge.relative_index[0] * m_a1[0];
      float dot = k * r;
      std::complex<float> phase = std::exp(std::complex<float>(0, dot));
      hamiltonian(edge.from_index, edge.to_index) = phase;
      hamiltonian(edge.to_index, edge.from_index) = std::conj(phase);
    }
    return hamiltonian;
  }

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

  void Plot() const override {
    std::ofstream sites_file("sites.dat");
    std::ofstream edges_file("edges.dat");
    std::ofstream lattice_file("lattice.dat");

    lattice_file << "#m=0,S=16\n";
    for (const Site &site : Sites()) {
      for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
          Vec<float> offset = m_a1 * i + m_a2 * j;
          Vec<float> p = site.position + offset;
          lattice_file << p[0] << " " << p[1] << "\n";
        }
      }
    }

    lattice_file << "\n#m=1,S=1\n";
    for (const Edge &edge : Edges()) {
      for (int i = -1; i <= 1; i++) {
        for (int j = -1; j <= 1; j++) {
          Vec<float> from_position =
              SiteAt(edge.from_index).position + m_a1 * i + m_a2 * j;
          Vec<float> to_position =
              SiteAt(edge.to_index).position + (m_a1 * i + m_a2 * j) +
              (m_a1 * edge.relative_index[0] + m_a2 * edge.relative_index[1]);
          lattice_file << from_position[0] << " " << from_position[1] << "\n";
          lattice_file << to_position[0] << " " << to_position[1] << "\n\n";
        }
      }
    }
  }

  Matrix<std::complex<float>> HoppingMatrix(float k) const override {
    Matrix<std::complex<float>> hamiltonian(Size(), Size());
    for (auto edge : Edges()) {
      float r1 = SiteAt(edge.from_index).position[0] +
                 edge.relative_index[0] * m_a1[0];
      float r2 =
          SiteAt(edge.to_index).position[0] + edge.relative_index[1] * m_a2[1];
      float dot = k * (r1 + r2);
      std::complex<float> phase = std::exp(std::complex<float>(0, dot));
      hamiltonian(edge.from_index, edge.to_index) = phase;
      hamiltonian(edge.to_index, edge.from_index) = std::conj(phase);
    }
    return hamiltonian;
  }

 protected:
  Vec<float> m_a1;
  Vec<float> m_a2;
};

#include "cctb/lattice.h"

#include <ostream>

Matrix<int> Lattice::AdjMatrix() const {
  Matrix<int> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    hamiltonian(edge.from_index, edge.to_index) = 1;
    hamiltonian(edge.to_index, edge.from_index) = 1;
  }
  return hamiltonian;
}

void OneDimensionalLattice::Plot(PainterBackend backend,
                                 std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);

  plotter->Prepare();

  for (const Site &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      Vec<float> offset = m_a1 * i;
      Vec<float> p = site.position + offset;
      plotter->DrawPoint(p[0], 0.0f);
    }
  }

  for (const Edge &edge : Edges()) {
    for (int i = -1; i <= 1; i++) {
      Vec<float> from_position = SiteAt(edge.from_index).position + m_a1 * i;
      Vec<float> to_position = SiteAt(edge.to_index).position + m_a1 * i +
                               m_a1 * edge.relative_index[0];
      plotter->DrawLine(from_position[0], 0.0f, to_position[0], 0.0f);
    }
  }

  plotter->Finish();
}

Matrix<std::complex<float>> OneDimensionalLattice::HoppingMatrix(
    Vec<float> k) const {
  Matrix<std::complex<float>> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    Vec<float> from_position = SiteAt(edge.from_index).position;
    Vec<float> to_position =
        SiteAt(edge.to_index).position + m_a1 * edge.relative_index[0];
    float dot = k.dot(to_position - from_position);
    std::complex<float> phase = std::complex<float>(0, dot);
    hamiltonian(edge.from_index, edge.to_index) +=
        edge.weight * std::exp(phase);
    hamiltonian(edge.to_index, edge.from_index) +=
        edge.weight * std::exp(-phase);
  }
  return hamiltonian;
}

void TwoDimensionalLattice::Plot(PainterBackend backend,
                                 std::ostream &out) const {
  auto plotter = PainterFactory::create(backend, out);

  plotter->Prepare();

  for (const Site &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Vec<float> offset = m_a1 * i + m_a2 * j;
        Vec<float> p = site.position + offset;
        plotter->DrawPoint(p[0], p[1]);
      }
    }
  }

  for (const Edge &edge : Edges()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Vec<float> from_position =
            SiteAt(edge.from_index).position + m_a1 * i + m_a2 * j;
        Vec<float> to_position =
            SiteAt(edge.to_index).position + (m_a1 * i + m_a2 * j) +
            (m_a1 * edge.relative_index[0] + m_a2 * edge.relative_index[1]);
        plotter->DrawLine(from_position[0], from_position[1], to_position[0],
                          to_position[1]);
      }
    }
  }

  plotter->DrawText(m_a1[0] / 2, m_a1[1] / 2 + 0.25, "$a_1$");
  plotter->DrawArrow(0, 0, m_a1[0], m_a1[1]);
  plotter->DrawText(m_a2[0] / 2, m_a2[1] / 2 - 0.25, "$a_2$");
  plotter->DrawArrow(0, 0, m_a2[0], m_a2[1]);

  plotter->Finish();
}

Matrix<std::complex<float>> TwoDimensionalLattice::HoppingMatrix(
    Vec<float> k) const {
  Matrix<std::complex<float>> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    Vec<float> from_position = SiteAt(edge.from_index).position;
    Vec<float> to_position = SiteAt(edge.to_index).position +
                             m_a1 * edge.relative_index[0] +
                             m_a2 * edge.relative_index[1];
    float dot = k.dot(to_position - from_position);
    std::complex<float> phase = std::complex<float>(0, dot);
    hamiltonian(edge.from_index, edge.to_index) +=
        edge.weight * std::exp(phase);
    hamiltonian(edge.to_index, edge.from_index) +=
        edge.weight * std::exp(-phase);
  }
  return hamiltonian;
}
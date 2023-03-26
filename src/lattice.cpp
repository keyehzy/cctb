#include "cctb/lattice.h"

Matrix<int> Lattice::AdjMatrix() const {
  Matrix<int> hamiltonian(Size(), Size());
  for (auto edge : Edges()) {
    hamiltonian(edge.from_index, edge.to_index) = 1;
    hamiltonian(edge.to_index, edge.from_index) = 1;
  }
  return hamiltonian;
}

void OneDimensionalLattice::Plot() const {
  std::ofstream lattice_file("lattice.dat");

  lattice_file << "#m=0,S=16\n";
  for (const Site &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      Vec<float> offset = m_a1 * i;
      Vec<float> p = site.position + offset;
      lattice_file << p[0] << " " << 0.0f << "\n";
    }
  }

  lattice_file << "\n#m=1,S=1\n";
  for (const Edge &edge : Edges()) {
    for (int i = -1; i <= 1; i++) {
      Vec<float> from_position = SiteAt(edge.from_index).position + m_a1 * i;
      Vec<float> to_position = SiteAt(edge.to_index).position + m_a1 * i +
                               m_a1 * edge.relative_index[0];
      lattice_file << from_position[0] << " " << 0.0f << "\n";
      lattice_file << to_position[0] << " " << 0.0f << "\n\n";
    }
  }
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

void TwoDimensionalLattice::Plot() const {
  std::ofstream lattice_file("lattice.tex");

  // Tikz backend
  lattice_file << "\\documentclass[tikz]{standalone}\n";
  lattice_file << "\\begin{document}\n";
  lattice_file << "\\begin{tikzpicture}\n";

  lattice_file << "\\draw[->] (0,0) -- (" << m_a1[0] << "," << m_a1[1]
               << ");\n";
  lattice_file << "\\draw[->] (0,0) -- (" << m_a2[0] << "," << m_a2[1]
               << ");\n";

  for (const Site &site : Sites()) {
    for (int i = -1; i <= 1; i++) {
      for (int j = -1; j <= 1; j++) {
        Vec<float> offset = m_a1 * i + m_a2 * j;
        Vec<float> p = site.position + offset;
        // lattice_file << p[0] << " " << p[1] << "\n";
        lattice_file << "\\fill[black] (" << p[0] << "," << p[1]
                     << ") circle (0.1);\n";
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
        // lattice_file << from_position[0] << " " << from_position[1] <<
        // "\n"; lattice_file << to_position[0] << " " << to_position[1] <<
        // "\n\n";
        lattice_file << "\\draw (" << from_position[0] << ","
                     << from_position[1] << ") -- (" << to_position[0] << ","
                     << to_position[1] << ");\n";
      }
    }
  }

  lattice_file << "\\end{tikzpicture}\n";
  lattice_file << "\\end{document}\n";

  // lattice_file << "#m=1,S=3\n";
  // lattice_file << 0.0 << " " << 0.0 << "\n";
  // lattice_file << m_a1[0] << " " << m_a1[1] << "\n\n";
  // lattice_file << 0.0 << " " << 0.0 << "\n";
  // lattice_file << m_a2[0] << " " << m_a2[1] << "\n\n";

  // lattice_file << "#m=0,S=16\n";
  // for (const Site &site : Sites()) {
  //   for (int i = -1; i <= 1; i++) {
  //     for (int j = -1; j <= 1; j++) {
  //       Vec<float> offset = m_a1 * i + m_a2 * j;
  //       Vec<float> p = site.position + offset;
  //       lattice_file << p[0] << " " << p[1] << "\n";
  //     }
  //   }
  // }

  // lattice_file << "\n#m=1,S=1\n";
  // for (const Edge &edge : Edges()) {
  //   for (int i = -1; i <= 1; i++) {
  //     for (int j = -1; j <= 1; j++) {
  //       Vec<float> from_position =
  //           SiteAt(edge.from_index).position + m_a1 * i + m_a2 * j;
  //       Vec<float> to_position =
  //           SiteAt(edge.to_index).position + (m_a1 * i + m_a2 * j) +
  //           (m_a1 * edge.relative_index[0] + m_a2 *
  //           edge.relative_index[1]);
  //       lattice_file << from_position[0] << " " << from_position[1] <<
  //       "\n"; lattice_file << to_position[0] << " " << to_position[1] <<
  //       "\n\n";
  //     }
  //   }
  // }
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
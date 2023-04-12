#pragma once

#include <vector>

#include "Geometry/Point.h"

class Mesh {
  using Vertex = Point<2>;

 public:
  Mesh() = default;
  Mesh(const Mesh&) = default;
  explicit Mesh(std::vector<Vertex> vertices) : m_vertices(std::move(vertices)) {}

  const std::vector<Vertex>& vertices() const { return m_vertices; }

  void Add(const Vertex& vertex) { m_vertices.push_back(vertex); }

  void Combine(const Mesh& mesh) {
    m_vertices.reserve(m_vertices.size() + mesh.m_vertices.size());
    m_vertices.insert(m_vertices.end(), mesh.m_vertices.begin(), mesh.m_vertices.end());
  }

  Vertex& operator[](size_t index) { return m_vertices[index]; }

  const Vertex& operator[](size_t index) const { return m_vertices[index]; }

 private:
  std::vector<Vertex> m_vertices;
};

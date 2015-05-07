#ifndef SRC_MATERIAL_H_
#define SRC_MATERIAL_H_

#include <glm/glm.hpp>

class Material {
 public:
  Material(float mp, glm::vec4 c) : melting_point(mp), base_color(c) {}
  // The point at which the vertices begin to move
  float melting_point;
  // The normal color of a vertex
  glm::vec4 base_color;
};

#endif  // SRC_MATERIAL_H_

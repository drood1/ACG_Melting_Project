#ifndef SRC_HEAT_SOURCE_H_
#define SRC_HEAT_SOURCE_H_

#include <glm/glm.hpp>

class HeatSource {
 public:
  HeatSource(float ho, glm::vec3 p) : heat_output(ho), position(p) {}
  // The amount of heat that is emitted. The units are J/s
  float heat_output;
  // The location of the heat source
  glm::vec3 position;
  // By default the color is red, likely not ever going to change it
  glm::vec4 color = glm::vec4(1, 0, 0, 1);
};

#endif  // SRC_HEAT_SOURCE_H_

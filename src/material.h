#ifndef SRC_MATERIAL_H_
#define SRC_MATERIAL_H_

#include <glm/glm.hpp>

class Material {
 public:
  Material(float mp, float sh, float sd, glm::vec4 c) : melting_point(mp), specific_heat(sh), solidity(sd), base_color(c) {}
  // The point at which the vertices begin to move in degrees K
  float melting_point;
  // Resistance to gaining heat, a larger number is higher resistance (ie heat is gained less)
  // Unit is J/gK, as in this amount of energy in Joules will raise
  // one gram of the material by one degree Kelvin
  float specific_heat;
  // Resistance to moving because of heat, 0 indicates it will never move and 100 indicates it will be fully liquid as soon
  // as it melts
  float solidity;
  // The color of the material
  glm::vec4 base_color;
};

#endif  // SRC_MATERIAL_H_

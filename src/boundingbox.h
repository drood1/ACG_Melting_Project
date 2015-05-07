#ifndef SRC_BOUNDINGBOX_H_
#define SRC_BOUNDINGBOX_H_

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cassert>
#include <algorithm>

// needed by Windows
// allows us to use std::min & std::max
#define NOMINMAX

class BoundingBox {
 public:
  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  BoundingBox() {
    Set(glm::vec3(0, 0, 0), glm::vec3(0, 0, 0)); }
  BoundingBox(const glm::vec3 &_minimum, const glm::vec3 &_maximum) {
    Set(_minimum, _maximum); }

  // =========
  // ACCESSORS
  void Get(glm::vec3 &_minimum, glm::vec3 &_maximum) const {
    _minimum = minimum;
    _maximum = maximum; }
  glm::vec3 getMin() const { return minimum; }
  glm::vec3 getMax() const { return maximum; }
  void getCenter(glm::vec3 &c) const {
    c = maximum;
    c -= minimum;
    c *= 0.5f;
    c += minimum;
  }
  glm::vec3 getCenter() const {
    glm::vec3 c = maximum;
    c -= minimum;
    c *= 0.5f;
    c += minimum;
    return c;
  }
  double maxDim() const {
    double x = maximum.x - minimum.x;
    double y = maximum.y - minimum.y;
    double z = maximum.z - minimum.z;
    return std::max(x, std::max(y, z));
  }

  // =========
  // MODIFIERS
  void Set(const BoundingBox &bb) {
    minimum = bb.minimum;
    maximum = bb.maximum; }
  void Set(const glm::vec3 &_minimum, const glm::vec3 &_maximum) {
    assert(minimum.x <= maximum.x &&
      minimum.y <= maximum.y &&
      minimum.z <= maximum.z);
    minimum = _minimum;
    maximum = _maximum; }
  void Extend(const glm::vec3 v) {
    minimum = glm::vec3(std::min(minimum.x, v.x),
                        std::min(minimum.y, v.y),
                        std::min(minimum.z, v.z));
    maximum = glm::vec3(std::max(maximum.x, v.x),
                        std::max(maximum.y, v.y),
                        std::max(maximum.z, v.z));
  }
  void Extend(const BoundingBox &bb) {
    Extend(bb.minimum);
    Extend(bb.maximum);
  }

 private:
  // ==============
  // REPRESENTATION
  glm::vec3 minimum;
  glm::vec3 maximum;
};

#endif  // SRC_BOUNDINGBOX_H_

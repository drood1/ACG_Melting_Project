#ifndef SRC_VERTEX_H_
#define SRC_VERTEX_H_

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdlib>

#include "./material.h"

// ==========================================================
// Stores the vertex position, used by the Mesh class
class Vertex {
 public:
  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const glm::vec3 &pos) : position(pos), original_position(pos), heat(0) { index = i; }

  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x; }
  double y() const { return position.y; }
  double z() const { return position.z; }
  const glm::vec3& getPos() const { return position; }
  const glm::vec3& getOriginalPos() const { return original_position; }
  Material* getMaterial() const { return material; }
  float getHeat() { return heat; }
  glm::vec3 getVelocity() {
    // Velocity is heat times gravity
    glm::vec3 gravity(0.0f, -9.8f, 0.0f);
    if (heat >= material->melting_point) {
      return heat * gravity;
    } else {
      return glm::vec3(0.0f, 0.0f, 0.0f);
    }
  }
  bool isOnFloor() {
    return is_on_floor;
  }
  void setOnFloor() {
    is_on_floor = true;
  }
  bool isMelting() {
    return heat >= material->melting_point;
  }

  glm::vec4 getHeatColor() {
    // Red(ish)
    if (heat >= material->melting_point) {
      return heat * 200.0f * glm::vec4(1.0, 0.0, 0.0, 1.0);
    } else {
      return glm::vec4(0.0, 0.0, 1.0, 1.0);
    }
  }

  glm::vec4 getColor(int mode) {
    if (mode == 0) {
      return material->base_color;
    } else {
      return getHeatColor();
    }
  }

  // =========
  // MODIFIERS
  void setPos(glm::vec3 v) { position = v; }
  void setHeat(float h) { heat = h; }
  void setMaterial(Material* m) { material = m; }
  void loseHeat(float heat_loss, float timestep) { heat -= heat * heat_loss * timestep; }

 private:
  // don't use these constructors
  Vertex() { assert(0); exit(0); }
  Vertex(const Vertex&) { assert(0); exit(0); }
  Vertex& operator=(const Vertex&) { assert(0); exit(0); }

  // ==============
  // REPRESENTATION
  glm::vec3 position;
  glm::vec3 original_position;
  // Some value of heat/energy, Could probably just be temperature
  // but for now I am using this as a fraction of gravity
  // valid values are [0, 1]
  // Velocity of vertex is modified using this
  float heat;
  Material* material;
  bool is_on_floor = false;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure,
  // but we use it for hashing
  int index;

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation, so it has been omitted.
};

#endif  // SRC_VERTEX_H_

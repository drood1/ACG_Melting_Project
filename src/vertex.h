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
  Vertex(int i, float t, const glm::vec3 &pos) : position(pos), original_position(pos), temperature(293) {
    index = i;
    force = glm::vec3(0.0, 0.0, 0.0);
    velocity = glm::vec3(0.0, 0.0, 0.0);
  }

  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x; }
  double y() const { return position.y; }
  double z() const { return position.z; }
  const glm::vec3& getPos() const { return position; }
  const glm::vec3& getOriginalPos() const { return original_position; }
  const glm::vec3& getForce() const { return force; }
  const glm::vec3& getForceUnder() const { return force_under; }
  const glm::vec3& getVelocity() const { return velocity; }
  Material* getMaterial() const { return material; }
  float getMass() { return 10000090.0f; }
  float getTemperature() { return temperature; }
  float getEnergy() { return energy; }
  bool isOnFloor() { return is_on_floor; }
  void setOnFloor() { is_on_floor = true; }
  bool isMelting() { return temperature >= material->melting_point; }
  glm::vec4 getColor(int mode) {
    if (mode == 0) {
      return material->base_color;
    } else if (mode == 1) {
      if (isMelting()) {
        // Maximum temperature is 600K
        return (temperature - material->melting_point) / material->melting_point * glm::vec4(1.0, 0.0, 0.0, 1.0);
      } else {
        return glm::vec4(0.0, 0.0, 1.0, 1.0);
      }
    } else {
      // Force mode
      return glm::length(force) / 100.0f * glm::vec4(0.0, 1.0, 0.0, 1.0);
    }
  }

  // =========
  // MODIFIERS
  void setPos(glm::vec3 v) { position = v; }
  void setTemperature(float t) { temperature = t; }
  void setMaterial(Material* m) { material = m; }
  void setForce(glm::vec3 f) { force = f; }
  void setVelocity(glm::vec3 v) { velocity = v; }
  void setForceUnder(glm::vec3 fu) { force_under = fu; }
  // void loseHeat(float heat_loss, float timestep) { heat -= heat * heat_loss * timestep; }
  void addEnergy(float e) { energy += e; }
  void receiveEnergy() {
    temperature += energy / material->specific_heat;
    energy = 0.0;
  }

 private:
  // don't use these constructors
  Vertex() { assert(0); exit(0); }
  Vertex(const Vertex&) { assert(0); exit(0); }
  Vertex& operator=(const Vertex&) { assert(0); exit(0); }

  // ==============
  // REPRESENTATION
  glm::vec3 position;
  glm::vec3 original_position;
  glm::vec3 force;
  glm::vec3 velocity;
  glm::vec3 force_under;
  // Some value of heat/energy, Could probably just be temperature
  // but for now I am using this as a fraction of gravity
  // valid values are [0, 1]
  // Velocity of vertex is modified using this
  float temperature;
  float energy;  // Extra energy that has not been converted to temperature
  Material* material;
  bool is_on_floor = false;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure,
  // but we use it for hashing
  int index;
};

#endif  // SRC_VERTEX_H_

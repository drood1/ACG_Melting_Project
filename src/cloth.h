#ifndef _CLOTH_H_
#define _CLOTH_H_

#include "argparser.h"
#include "boundingbox.h"
#include "vbo_structs.h"
#include <vector>
#include <iostream>

// =====================================================================================
// Cloth Particles
// =====================================================================================

class ClothParticle {
public:
  // ACCESSORS
  const glm::vec3& getOriginalPosition() const{ return original_position; }
  const glm::vec3& getPosition() const{ return position; }
  const glm::vec3& getVelocity() const{ return velocity; }
  const glm::vec3& getAcceleration() const { return acceleration; }
  glm::vec3 getForce() const { return float(mass)*acceleration; }
  double getMass() const { return mass; }
  bool isFixed() const { return fixed; }
  // MODIFIERS
  void setOriginalPosition(const glm::vec3 &p) { original_position = p; }
  void setPosition(const glm::vec3 &p) { position = p; }
  void setVelocity(const glm::vec3 &v) { velocity = v; }
  void setAcceleration(const glm::vec3 &a) { acceleration = a; }
  void setMass(double m) { mass = m; }
  void setFixed(bool b) { fixed = b; }
private:
  // REPRESENTATION
  glm::vec3 original_position;
  glm::vec3 position;
  glm::vec3 velocity;
  glm::vec3 acceleration;
  double mass;
  bool fixed;
};

// =====================================================================================
// Cloth System
// =====================================================================================

class Cloth {

public:
  Cloth(ArgParser *args);
  ~Cloth() { delete [] particles; cleanupVBOs(); }

  // ACCESSORS
  const BoundingBox& getBoundingBox() const { return box; }

  // PAINTING & ANIMATING
  void Paint() const;
  void Animate();

  void initializeVBOs();
  void setupVBOs();
  void drawVBOs();
  void cleanupVBOs();

private:

  // PRIVATE ACCESSORS
  const ClothParticle& getParticle(int i, int j) const {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }
  ClothParticle& getParticle(int i, int j) {
    assert (i >= 0 && i < nx && j >= 0 && j < ny);
    return particles[i + j*nx]; }

  glm::vec3 computeGouraudNormal(int i, int j) const;

  // HELPER FUNCTION
  void computeBoundingBox();

  // HELPER FUNCTIONS FOR ANIMATION
  void AddWireFrameTriangle(const glm::vec3 &apos, const glm::vec3 &bpos, const glm::vec3 &cpos,
                            const glm::vec3 &anormal, const glm::vec3 &bnormal, const glm::vec3 &cnormal,
                            const glm::vec3 &abcolor, const glm::vec3 &bccolor, const glm::vec3 &cacolor);

  glm::vec3 CalculateSpringForce(const ClothParticle& initial, const ClothParticle& ending, float k) {
    float original_distance = glm::distance(initial.getOriginalPosition(), ending.getOriginalPosition());
    float current_distance = glm::distance(initial.getPosition(), ending.getPosition());
    float deformation = original_distance - current_distance;
    glm::vec3 direction = glm::normalize(initial.getPosition() - ending.getPosition());
    glm::vec3 force = k * deformation * direction;
    // if (k == k_structural && current_distance != original_distance) {
    //   std::cout << "k: " << k << std::endl;
    //   std::cout << "original_distance: " << original_distance << std::endl;
    //   std::cout << "current_distance: " << current_distance << std::endl;
    //   std::cout << "deformation: " << deformation << std::endl;
    //   std::cout << "direction: " << direction[0] << ", " << direction[1] << ", " << direction[2] << std::endl;
    //   std::cout << "force: " << force[0] << ", " << force[1] << ", " << force[2] << std::endl << std::endl;
    // }
    return force;
  }

  void ProvotCorrection(ClothParticle& initial, ClothParticle& ending, float correction) {
    float original_distance = glm::distance(initial.getOriginalPosition(), ending.getOriginalPosition());
    float current_distance = glm::distance(initial.getPosition(), ending.getPosition());
    float deformation = original_distance - current_distance;
    float deformation_rate = deformation / original_distance;
    if (std::abs(deformation_rate) > correction) {
      // Need to get deformation_rate == correction
      glm::vec3 direction = glm::normalize(initial.getPosition() - ending.getPosition());
      float distance = deformation - correction * original_distance;
      // std::cout << distance << std::endl;
      if (ending.isFixed()) {
        glm::vec3 position = initial.getPosition();
        position += direction * distance;
        initial.setPosition(position);
      } else {
        glm::vec3 position = initial.getPosition();
        position += direction * distance / 2.0f;
        initial.setPosition(position);

        position = ending.getPosition();
        position -= direction * distance / 2.0f;
        ending.setPosition(position);
      }
    }
  }


  // REPRESENTATION
  ArgParser *args;
  // grid data structure
  int nx, ny;
  ClothParticle *particles;
  BoundingBox box;
  // simulation parameters
  double damping;
  // spring constants
  double k_structural;
  double k_shear;
  double k_bend;
  // correction thresholds
  double provot_structural_correction;
  double provot_shear_correction;

  // VBOs
  GLuint cloth_verts_VBO;
  GLuint cloth_tri_indices_VBO;
  GLuint cloth_velocity_verts_VBO;
  GLuint cloth_velocity_tri_indices_VBO;
  GLuint cloth_force_verts_VBO;
  GLuint cloth_force_tri_indices_VBO;

  // data for rendering
  std::vector<VBOPosNormalColor> cloth_verts; 
  std::vector<VBOIndexedTri> cloth_tri_indices;
  std::vector<VBOPosNormalColor> cloth_velocity_verts; 
  std::vector<VBOIndexedTri> cloth_velocity_tri_indices;
  std::vector<VBOPosNormalColor> cloth_force_verts; 
  std::vector<VBOIndexedTri> cloth_force_tri_indices;
};

// ========================================================================

#endif

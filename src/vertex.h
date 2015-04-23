#ifndef _VERTEX_H
#define _VERTEX_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <cstdlib>

// ==========================================================
// Stores the vertex position, used by the Mesh class

class Vertex {

public:

  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  Vertex(int i, const glm::vec3 &pos) : position(pos) { index = i; }
  
  // =========
  // ACCESSORS
  int getIndex() const { return index; }
  double x() const { return position.x; }
  double y() const { return position.y; }
  double z() const { return position.z; }
  const glm::vec3& getPos() const { return position; }
  float getHeat() { return heat; }
  glm::vec3 getVelocity() {
    // Velocity is heat times gravity
    glm::vec3 gravity(0.0f, -9.8f, 0.0f);
    return heat * gravity;
  }


  // =========
  // MODIFIERS
  void setPos(glm::vec3 v) { position = v; }
  void setHeat(float h) { heat = h; }

private:

  // don't use these constructors
  Vertex() { assert(0); exit(0); }
  Vertex(const Vertex&) { assert(0); exit(0); }
  Vertex& operator=(const Vertex&) { assert(0); exit(0); }
  
  // ==============
  // REPRESENTATION
  glm::vec3 position;
  // Some value of heat/energy, Could probably just be temperature
  // but for now I am using this as a fraction of gravity
  // valid values are [0, 1]
  // Velocity of vertex is modified using this
  float heat;

  // this is the index from the original .obj file.
  // technically not part of the half-edge data structure, 
  // but we use it for hashing
  int index;  

  // NOTE: the vertices don't know anything about adjacency.  In some
  // versions of this data structure they have a pointer to one of
  // their incoming edges.  However, this data is very complicated to
  // maintain during mesh manipulation, so it has been omitted.

};

// ==========================================================

#endif


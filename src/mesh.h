#ifndef SRC_MESH_H_
#define SRC_MESH_H_

#include <cstdlib>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include "./glCanvas.h"
#include "./hash.h"
#include "./boundingbox.h"
#include "./argparser.h"
#include "./mtrand.h"
#include "./material.h"
#include "./heat_source.h"
#include "./vbo_structs.h"


class Vertex;
class Edge;
class Triangle;

// Stores and renders all the vertices, triangles, and edges for a 3D model
class Mesh {
 public:
  // ========================
  // CONSTRUCTOR & DESTRUCTOR
  explicit Mesh(ArgParser *a) { args = a; }
  ~Mesh();
  void Load(const std::string &input_file);

  // ========
  // VERTICES
  int numVertices() const { return vertices.size(); }
  Vertex* addVertex(const glm::vec3 &pos);
  // look up vertex by index from original .obj file
  Vertex* getVertex(int i) const {
    assert(i >= 0 && i < numVertices());
    Vertex *v = vertices[i];
    assert(v != NULL);
    return v; }

  void removeHeat() {
    // heat_position += glm::vec3(1000.0, 1000.0, 1000.0);
  }

  void replaceHeat() {
    // heat_position -= glm::vec3(1000.0, 1000.0, 1000.0);
  }

  // ==================================================
  // PARENT VERTEX RELATIONSHIPS (used for subdivision)
  // this creates a relationship between 3 vertices (2 parents, 1 child)
  void setParentsChild(Vertex *p1, Vertex *p2, Vertex *child);
  // this accessor will find a child vertex (if it exists) when given
  // two parent vertices
  Vertex* getChildVertex(Vertex *p1, Vertex *p2) const;
  Vertex* getCenter(Edge* edge);

  // =====
  // EDGES
  int numEdges() const { return edges.size(); }
  // this efficiently looks for an edge with the given vertices, using a hash table
  Edge* getMeshEdge(Vertex *a, Vertex *b) const;

  // =========
  // TRIANGLES
  int numTriangles() const { return triangles.size(); }
  void addTriangle(Vertex *a, Vertex *b, Vertex *c);
  void removeTriangle(Triangle *t);

  // ===============
  // OTHER ACCESSORS
  const BoundingBox& getBoundingBox() const { return bbox; }

  // ===+=====
  // RENDERING

  // NEW RENDERING
  void initializeVBOs();
  void setupVBOs();
  void SetupMesh();
  void SetupFloor();
  void SetupHeat();
  void drawVBOs(const glm::mat4 &ProjectionMatrix,
    const glm::mat4 &ViewMatrix, const glm::mat4 &ModelMatrix);
  void cleanupVBOs();
  void animate();
  void animateHeat();

 private:
  // don't use these constructors
  Mesh(const Mesh &/*m*/) { assert(0); exit(0); }
  const Mesh& operator=(const Mesh &/*m*/) { assert(0); exit(0); }

  // ==============
  // REPRESENTATION
  ArgParser *args;
  std::vector<Vertex*> vertices;
  edgeshashtype edges;
  triangleshashtype triangles;
  BoundingBox bbox;
  vphashtype vertex_parents;
  int num_mini_triangles;
  MTRand_closed rand;
  std::set<Edge*> special_edges;
  float floor_y;
  float heat_loss;
  float room_temperature = 293.0f;
  float floor_factor = 0.05;
  float min_x, max_x, min_z, max_z;
  std::map<std::string, Material*> materials;
  std::string material;
  std::vector<HeatSource*> heat_sources;
  std::map<Vertex*, Edge*> edges_by_top;

  GLuint mesh_VAO;
  GLuint mesh_tri_verts_VBO;
  GLuint mesh_tri_indices_VBO;
  GLuint floor_tri_verts_VBO;
  GLuint floor_tri_indices_VBO;
  GLuint heat_vert_VBO;

  std::vector<VBOPosNormalColor> mesh_tri_verts;
  std::vector<VBOIndexedTri> mesh_tri_indices;
  std::vector<VBOPosNormalColor> floor_tri_verts;
  std::vector<VBOIndexedTri> floor_tri_indices;
  std::vector<VBOPosNormalColor> heat_vert;
};

#endif  // SRC_MESH_H_

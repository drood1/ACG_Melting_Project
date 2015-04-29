#include "./mesh.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <set>
#include <map>
#include <queue>
#include <string>
#include <vector>

#include "./glCanvas.h"
#include "./edge.h"
#include "./vertex.h"
#include "./triangle.h"
#include "./mtrand.h"

// to give a unique id number to each triangles
int Triangle::next_triangle_id = 0;

Mesh::~Mesh() {
  cleanupVBOs();

  // delete all the triangles
  std::vector<Triangle*> todo;
  for (triangleshashtype::iterator iter = triangles.begin();
       iter != triangles.end(); iter++) {
    Triangle *t = iter->second;
    todo.push_back(t);
  }
  int num_triangles = todo.size();
  for (int i = 0; i < num_triangles; i++) {
    removeTriangle(todo[i]);
  }
  // delete all the vertices
  int num_vertices = numVertices();
  for (int i = 0; i < num_vertices; i++) {
    delete vertices[i];
  }
}

Vertex* Mesh::addVertex(const glm::vec3 &position) {
  int index = numVertices();
  Vertex *v = new Vertex(index, position);
  vertices.push_back(v);
  if (numVertices() == 1)
    bbox = BoundingBox(position, position);
  else
    bbox.Extend(position);
  return v;
}

void Mesh::addTriangle(Vertex *a, Vertex *b, Vertex *c) {
  // create the triangle
  Triangle *t = new Triangle();
  // create the edges
  Edge *ea = new Edge(a, b, t);
  Edge *eb = new Edge(b, c, t);
  Edge *ec = new Edge(c, a, t);
  // point the triangle to one of its edges
  t->setEdge(ea);
  // connect the edges to each other
  ea->setNext(eb);
  eb->setNext(ec);
  ec->setNext(ea);
  // verify these edges aren't already in the mesh
  // (which would be a bug, or a non-manifold mesh)
  assert(edges.find(std::make_pair(a, b)) == edges.end());
  assert(edges.find(std::make_pair(b, c)) == edges.end());
  assert(edges.find(std::make_pair(c, a)) == edges.end());
  // add the edges to the master list
  edges[std::make_pair(a, b)] = ea;
  edges[std::make_pair(b, c)] = eb;
  edges[std::make_pair(c, a)] = ec;
  // connect up with opposite edges (if they exist)
  edgeshashtype::iterator ea_op = edges.find(std::make_pair(b, a));
  edgeshashtype::iterator eb_op = edges.find(std::make_pair(c, b));
  edgeshashtype::iterator ec_op = edges.find(std::make_pair(a, c));
  if (ea_op != edges.end()) { ea_op->second->setOpposite(ea); }
  if (eb_op != edges.end()) { eb_op->second->setOpposite(eb); }
  if (ec_op != edges.end()) { ec_op->second->setOpposite(ec); }
  // add the triangle to the master list
  assert(triangles.find(t->getID()) == triangles.end());
  triangles[t->getID()] = t;
}

void Mesh::removeTriangle(Triangle *t) {
  Edge *ea = t->getEdge();
  Edge *eb = ea->getNext();
  Edge *ec = eb->getNext();
  Vertex *a = ea->getStartVertex();
  Vertex *b = eb->getStartVertex();
  Vertex *c = ec->getStartVertex();
  // remove these elements from master lists
  edges.erase(std::make_pair(a, b));
  edges.erase(std::make_pair(b, c));
  edges.erase(std::make_pair(c, a));
  triangles.erase(t->getID());
  // clean up memory
  delete ea;
  delete eb;
  delete ec;
  delete t;
}

Edge* Mesh::getMeshEdge(Vertex *a, Vertex *b) const {
  edgeshashtype::const_iterator iter = edges.find(std::make_pair(a, b));
  if (iter == edges.end()) return NULL;
  return iter->second;
}

Vertex* Mesh::getChildVertex(Vertex *p1, Vertex *p2) const {
  vphashtype::const_iterator iter = vertex_parents.find(std::make_pair(p1, p2));
  if (iter == vertex_parents.end()) return NULL;
  return iter->second;
}

void Mesh::setParentsChild(Vertex *p1, Vertex *p2, Vertex *child) {
  assert(vertex_parents.find(std::make_pair(p1, p2)) == vertex_parents.end());
  vertex_parents[std::make_pair(p1, p2)] = child;
}

Vertex* Mesh::getCenter(Edge* edge) {
  Vertex* vertex = getChildVertex(edge->getStartVertex(), edge->getEndVertex());
  if (vertex == NULL) {
    vertex = getChildVertex(edge->getEndVertex(), edge->getStartVertex());
  }
  return vertex;
}

#define MAX_CHAR_PER_LINE 200
void Mesh::Load(const std::string &input_file) {
  std::ifstream istr(input_file.c_str());
  if (!istr) {
    std::cout << "ERROR! CANNOT OPEN: " << input_file << std::endl;
    return;
  }

  char line[MAX_CHAR_PER_LINE];
  std::string token, token2;
  float x, y, z;
  int a, b, c;
  int index = 0;
  int vert_count = 0;
  int vert_index = 1;

  // TODO(austin): maybe this shouldn't be hardcoded
  heat_position = glm::vec3(-0.2, 0.2, 0.1);

  // read in each line of the file
  while (istr.getline(line, MAX_CHAR_PER_LINE)) {
    // put the line into a stringstream for parsing
    std::stringstream ss;
    ss << line;

    // check for blank line
    token = "";
    ss >> token;
    if (token == "") continue;

    if (token == std::string("usemtl") || token == std::string("g")) {
      vert_index = 1;
      index++;
    } else if (token == std::string("v")) {
      vert_count++;
      ss >> x >> y >> z;
      Vertex* newVertex = addVertex(glm::vec3(x, y, z));
      setHeat(newVertex);
    } else if (token == std::string("f")) {
      a = b = c = -1;
      ss >> a >> b;
      // handle faces with > 3 vertices
      // assume the face can be triangulated with a triangle fan
      while (ss >> c) {
        int a_ = a-vert_index;
        int b_ = b-vert_index;
        int c_ = c-vert_index;
        assert(a_ >= 0 && a_ < numVertices());
        assert(b_ >= 0 && b_ < numVertices());
        assert(c_ >= 0 && c_ < numVertices());
        addTriangle(getVertex(a_), getVertex(b_), getVertex(c_));
        b = c;
      }
    } else if (token == std::string("e")) {
      a = b = -1;
      ss >> a >> b >> token2;
      // whoops: inconsistent file format, don't subtract 1
      assert(a >= 0 && a <= numVertices());
      assert(b >= 0 && b <= numVertices());
      if (token2 == std::string("inf")) x = 1000000;  // this is close to infinity...
      x = atof(token2.c_str());
      Vertex *va = getVertex(a);
      Vertex *vb = getVertex(b);
      Edge *ab = getMeshEdge(va, vb);
      Edge *ba = getMeshEdge(vb, va);
      assert(ab != NULL);
      assert(ba != NULL);
      ab->setCrease(x);
      ba->setCrease(x);
    } else if (token == std::string("vt")) {
    } else if (token == std::string("vn")) {
    } else if (token[0] == '#') {
    } else {
      printf("LINE: '%s'", line);
    }
  }

  std::cout << "Loaded " << numTriangles() << " triangles." << std::endl;

  floorY = vertices[0]->getPos()[1];
  for (int i = 0; i < vertices.size(); ++i) {
    Vertex* v = vertices[i];
    float y = v->getPos()[1];
    if (y < floorY) {
      floorY = y;
    }
  }

  assert(numTriangles() > 0);
  num_mini_triangles = 0;
}


// =======================================================================
// DRAWING
// =======================================================================

glm::vec3 ComputeNormal(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
  glm::vec3 v12 = p2;
  v12 -= p1;
  glm::vec3 v23 = p3;
  v23 -= p2;
  glm::vec3 normal = glm::cross(v12, v23);
  normal = glm::normalize(normal);
  return normal;
}


void Mesh::initializeVBOs() {
  HandleGLError("enter initialize VBOs");

  // create a pointer for the vertex & index VBOs
  glGenVertexArrays(1, &mesh_VAO);
  glBindVertexArray(mesh_VAO);
  glGenBuffers(1, &mesh_tri_verts_VBO);
  glGenBuffers(1, &mesh_tri_indices_VBO);
  glGenBuffers(1, &heat_vert_VBO);
  // and the data to pass to the shaders
  GLCanvas::MatrixID = glGetUniformLocation(GLCanvas::programID, "MVP");
  GLCanvas::LightID = glGetUniformLocation(GLCanvas::programID, "LightPosition_worldspace");
  GLCanvas::ViewMatrixID = glGetUniformLocation(GLCanvas::programID, "V");
  GLCanvas::ModelMatrixID = glGetUniformLocation(GLCanvas::programID, "M");

  GLCanvas::wireframeID = glGetUniformLocation(GLCanvas::programID, "wireframe");

  // call this the first time...
  setupVBOs();
  HandleGLError("leaving initializeVBOs");
}


// boundary edges are red, crease edges are yellow
glm::vec3 Mesh::EdgeColor(Edge *e) {
  if (special_edges.size() > 0 && special_edges.count(e) != 0) {
    std::cout << "SEPCI" << std::endl;
    return glm::vec3(0, 1, 0);
  } else if (e->getOpposite() == NULL) {
    return glm::vec3(1, 0, 0);
  } else if (e->getCrease() > 0) {
    return glm::vec3(1, 1, 0);
  } else {
    return glm::vec3(0, 0, 0.0);
  }
  return glm::vec3(0, 0, 0.0);
}


void Mesh::TriVBOHelper(std::vector<glm::vec3> &indexed_verts,
                        std::vector<unsigned int> &mesh_tri_indices,
                        const glm::vec3 &pos_a,
                        const glm::vec3 &pos_b,
                        const glm::vec3 &pos_c,
                        const glm::vec3 &normal_a,
                        const glm::vec3 &normal_b,
                        const glm::vec3 &normal_c,
                        const glm::vec3 &color_ab,
                        const glm::vec3 &color_bc,
                        const glm::vec3 &color_ca) {
  /*
  // To create a wireframe rendering...
  // Each mesh triangle is actually rendered as 3 small triangles
  //           b
  //          /|\
  //         / | \
  //        /  |  \
  //       /   |   \
  //      /    |    \
  //     /    .'.    \
  //    /  .'     '.  \
  //   /.'           '.\
  //  a-----------------c
  //
  */

  // the center is white, the colors of the two vertices depend on
  // whether the edge is a boundary edge (red) or crease edge (yellow)
  glm::vec3 center_color(1, 1, 1);
  // use simple averaging to find centroid & average normal
  glm::vec3 centroid = 1.0f / 3.0f * (pos_a + pos_b + pos_c);
  glm::vec3 normal = normal_a + normal_b + normal_c;
  normal = glm::normalize(normal);
  int i = indexed_verts.size()/3;

  if (args->wireframe) {
    // WIREFRAME
    // make the 3 small triangles
    indexed_verts.push_back(pos_a);
    indexed_verts.push_back(normal_a);
    indexed_verts.push_back(color_ab);
    indexed_verts.push_back(pos_b);
    indexed_verts.push_back(normal_b);
    indexed_verts.push_back(color_ab);
    indexed_verts.push_back(centroid);
    indexed_verts.push_back(normal);
    indexed_verts.push_back(center_color);

    indexed_verts.push_back(pos_b);
    indexed_verts.push_back(normal_b);
    indexed_verts.push_back(color_bc);
    indexed_verts.push_back(pos_c);
    indexed_verts.push_back(normal_c);
    indexed_verts.push_back(color_bc);
    indexed_verts.push_back(centroid);
    indexed_verts.push_back(normal);
    indexed_verts.push_back(center_color);

    indexed_verts.push_back(pos_c);
    indexed_verts.push_back(normal_c);
    indexed_verts.push_back(color_ca);
    indexed_verts.push_back(pos_a);
    indexed_verts.push_back(normal_a);
    indexed_verts.push_back(color_ca);
    indexed_verts.push_back(centroid);
    indexed_verts.push_back(normal);
    indexed_verts.push_back(center_color);

    // add all of the triangle vertices to the indices list
    for (int j = 0; j < 9; j++) {
      mesh_tri_indices.push_back(i+j);
    }
  } else {
    // NON WIREFRAME
    // Note: gouraud shading with the mini triangles looks bad... :(

    // make the 3 small triangles
    indexed_verts.push_back(pos_a);
    indexed_verts.push_back(normal_a);
    indexed_verts.push_back(center_color);
    indexed_verts.push_back(pos_b);
    indexed_verts.push_back(normal_b);
    indexed_verts.push_back(center_color);
    indexed_verts.push_back(pos_c);
    indexed_verts.push_back(normal_c);
    indexed_verts.push_back(center_color);

    // add all of the triangle vertices to the indices list
    for (int j = 0; j < 3; j++) {
      mesh_tri_indices.push_back(i+j);
    }
  }
}

void Mesh::setupVBOs() {
  HandleGLError("enter setupVBOs");

  std::vector<glm::vec3> indexed_verts;
  std::vector<unsigned int> mesh_tri_indices;

  // write the vertex & triangle data
  for (triangleshashtype::iterator iter = triangles.begin();
       iter != triangles.end(); iter++) {
    Triangle *t = iter->second;

    // grab the vertex positions
    glm::vec3 a = (*t)[0]->getPos();
    glm::vec3 b = (*t)[1]->getPos();
    glm::vec3 c = (*t)[2]->getPos();

    // determine edge colors (when wireframe is enabled)
    glm::vec3 edgecolor_ab = EdgeColor(t->getEdge());
    glm::vec3 edgecolor_bc = EdgeColor(t->getEdge()->getNext());
    glm::vec3 edgecolor_ca = EdgeColor(t->getEdge()->getNext()->getNext());

    if (args->gouraud) {
      // =====================================
      // ASSIGNMENT: complete this functionality
      // =====================================

      Edge *starting = t->getEdge();
      Edge *current = starting;
      glm::vec3 normal_a(0, 0, 0);
      float triangles = 0.0f;

      while (triangles == 0 ||
        ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_a += ComputeNormal((*triangle)[0]->getPos(),
          (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_a = normal_a / triangles;

      glm::vec3 normal_b(0, 0, 0);
      triangles = 0.0f;
      starting = t->getEdge()->getNext();
      current = starting;
      while (triangles == 0 ||
        ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_b += ComputeNormal((*triangle)[0]->getPos(),
          (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_b = normal_b / triangles;

      glm::vec3 normal_c(0, 0, 0);
      triangles = 0.0f;
      starting = t->getEdge()->getNext()->getNext();
      current = starting;
      while (triangles == 0 || ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_c += ComputeNormal((*triangle)[0]->getPos(),
          (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_c = normal_c / triangles;
      TriVBOHelper(indexed_verts, mesh_tri_indices,
                   a, b, c,
                   normal_a, normal_b, normal_c,
                   edgecolor_ab, edgecolor_bc, edgecolor_ca);
    } else {
      // for flat shading, use the triangle normal at each vertex
      // use the normal of the triangl
      glm::vec3 normal = ComputeNormal(a, b, c);
      TriVBOHelper(indexed_verts, mesh_tri_indices,
                   a, b, c,
                   normal, normal, normal,
                   edgecolor_ab, edgecolor_bc, edgecolor_ca);
    }
  }

  // the vertex data
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glBufferData(GL_ARRAY_BUFFER, indexed_verts.size() * sizeof(glm::vec3),
    &indexed_verts[0], GL_STATIC_DRAW);
  // the index data (refers to vertex data)
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices.size() * sizeof(unsigned int),
    &mesh_tri_indices[0] , GL_STATIC_DRAW);

  num_mini_triangles = mesh_tri_indices.size();

  heat_vert.clear();
  heat_vert.push_back(VBOPosNormalColor(heat_position, glm::vec3(1, 0, 0), glm::vec4(1, 1, 0, 0)));
  glBindBuffer(GL_ARRAY_BUFFER, heat_vert_VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(VBOPosNormalColor)*1,
    &heat_vert[0], GL_STATIC_DRAW);

  HandleGLError("leaving setupVBOs");
}


void Mesh::drawVBOs(const glm::mat4 &ProjectionMatrix,
  const glm::mat4 &ViewMatrix, const glm::mat4 &ModelMatrix) {
  HandleGLError("enter drawVBOs");

  // prepare data to send to the shaders
  glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

  glm::vec3 lightPos = glm::vec3(4, 4, 4);
  glUniform3f(GLCanvas::LightID, lightPos.x, lightPos.y, lightPos.z);
  glUniformMatrix4fv(GLCanvas::MatrixID, 1, GL_FALSE, &MVP[0][0]);
  glUniformMatrix4fv(GLCanvas::ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
  glUniformMatrix4fv(GLCanvas::ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
  glUniform1i(GLCanvas::wireframeID, args->wireframe);

  // triangle vertex positions
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,  // attribute
      3,                    // size
      GL_FLOAT,             // type
      GL_FALSE,             // normalized?
      3*sizeof(glm::vec3),  // stride
      reinterpret_cast<void*>(0));  // array buffer offset

  // triangle vertex normals
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,    // attribute
      3,                      // size
      GL_FLOAT,               // type
      GL_FALSE,               // normalized?
      3*sizeof(glm::vec3),    // stride
      reinterpret_cast<void*>(sizeof(glm::vec3)));
      // array buffer offset
  // triangle vertex colors
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2,        // attribute
      3,                          // size
      GL_FLOAT,                   // type
      GL_FALSE,                   // normalized?
      3*sizeof(glm::vec3),        // stride
      reinterpret_cast<void*>(sizeof(glm::vec3)*2));
      // array buffer offset
  // triangle indices
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glDrawElements(GL_TRIANGLES,          // mode
                 num_mini_triangles*3,  // count
                 GL_UNSIGNED_INT,       // type
                 reinterpret_cast<void*>(0));
                 // element array buffer offset
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  HandleGLError("enter draw heat");
  glPointSize(25);
  glBindBuffer(GL_ARRAY_BUFFER, heat_vert_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(0));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(sizeof(glm::vec3)));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(sizeof(glm::vec3)*2));
  glDrawArrays(GL_POINTS, 0, heat_vert.size());
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  HandleGLError("enter draw heat");

  HandleGLError("leaving drawVBOs");
}

void Mesh::cleanupVBOs() {
  glDeleteBuffers(1, &mesh_VAO);
  glDeleteBuffers(1, &mesh_tri_verts_VBO);
  glDeleteBuffers(1, &mesh_tri_indices_VBO);
  glDeleteBuffers(1, &heat_vert_VBO);
}

// Animation function
// FINAL PROJECT: a lot of our code is going to be here :)
void Mesh::animate() {
  if (args->animate) {
    // do 10 steps of animation before rendering
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        glm::vec3 position = v->getPos();
        position += args->timestep * v->getVelocity();
        if (position[1] < floorY) {
          position[1] = floorY;
        }
        v->setPos(position);
      }
    }
  }
  setupVBOs();
}

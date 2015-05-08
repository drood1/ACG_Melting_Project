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
#include <cmath>
#include <cstdlib>

#include "./glCanvas.h"
#include "./edge.h"
#include "./vertex.h"
#include "./triangle.h"
#include "./material.h"
#include "./mtrand.h"

void print_vec(glm::vec3 vec) {
  std::cout << vec[0] << ", " << vec[1] << ", " << vec[2] << std::endl;
}

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
  Vertex *v = new Vertex(index, room_temperature, position);
  v->setMaterial(materials[material]);
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
  ea->setOriginalDistance();
  eb->setOriginalDistance();
  ec->setOriginalDistance();
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

  // Default values
  material = "default";
  materials[material] = new Material(300, 1, 50, glm::vec4(0.5, 0.5, 0.5, 1.0)); 
  // heat_position = glm::vec3(-0.2, 0.2, 0.1);
  // heat_loss = 0.1;

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
    } else if (token == std::string("rt")) {
      ss >> x;
      room_temperature = x;      
    } else if (token == std::string("md")) {
      // Material definition
      float mp, sh, sd, r, g, b, a;
      ss >> token2 >> mp >> sh >> sd >> r >> g >> b >> a;
      materials[token2] = new Material(mp, sh, sd, glm::vec4(r, g, b, a));
      material = token2;
    } else if (token == std::string("m")) {
      // Use material
      ss >> token2;
      material = token2;
    } else if (token == std::string("h")) {
      // Heat source
      float ho;
      ss >> ho >> x >> y >> z;
      HeatSource* hs = new HeatSource(ho, glm::vec3(x, y, z));
      heat_sources.push_back(hs);
    } else if (token == std::string("hl")) {
      // Heat loss (cold room=more heat loss)
      ss >> x;
      heat_loss = x;
    } else if (token == std::string("v")) {
      vert_count++;
      ss >> x >> y >> z;
      addVertex(glm::vec3(x, y, z));
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

  floor_y = vertices[0]->getPos()[1];
  for (unsigned int i = 0; i < vertices.size(); ++i) {
    Vertex* v = vertices[i];
    float y = v->getPos()[1];
    if (y < floor_y) {
      floor_y = y;
    }
  }

  assert(numTriangles() > 0);
  num_mini_triangles = 0;

  for (edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    edges_by_top[edge->getTopVertex()] = edge;
  }
}

glm::vec3 ComputeNormal(const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3) {
  glm::vec3 v12 = p2;
  v12 -= p1;
  glm::vec3 v23 = p3;
  v23 -= p2;
  glm::vec3 normal = glm::cross(v12, v23);
  normal = glm::normalize(normal);
  return normal;
}

void Mesh::SetupMesh() {
  mesh_tri_verts.clear();
  mesh_tri_indices.clear();
  for (triangleshashtype::iterator iter = triangles.begin();
       iter != triangles.end(); iter++) {
    Triangle *t = iter->second;
    glm::vec3 a = (*t)[0]->getPos();
    glm::vec3 b = (*t)[1]->getPos();
    glm::vec3 c = (*t)[2]->getPos();
    glm::vec3 na = ComputeNormal(a, b, c);
    glm::vec3 nb = na;
    glm::vec3 nc = na;
    int start = mesh_tri_verts.size();
    mesh_tri_verts.push_back(VBOPosNormalColor(a, na, (*t)[0]->getColor(args->color_mode)));
    mesh_tri_verts.push_back(VBOPosNormalColor(b, nb, (*t)[1]->getColor(args->color_mode)));
    mesh_tri_verts.push_back(VBOPosNormalColor(c, nc, (*t)[2]->getColor(args->color_mode)));
    mesh_tri_indices.push_back(VBOIndexedTri(start, start+1, start+2));
  }
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glBufferData(GL_ARRAY_BUFFER,
         sizeof(VBOPosNormalColor) * mesh_tri_verts.size(),
         &mesh_tri_verts[0],
         GL_STATIC_DRAW);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
         sizeof(VBOIndexedTri) * mesh_tri_indices.size(),
         &mesh_tri_indices[0], GL_STATIC_DRAW);
}

void Mesh::SetupHeat() {
  heat_vert.clear();
  for (unsigned int i = 0; i < heat_sources.size(); ++i) {
    HeatSource* heat_source = heat_sources[i];
    heat_vert.push_back(VBOPosNormalColor(heat_source->position, glm::vec3(1, 1, 1), heat_source->color));
  }
  glBindBuffer(GL_ARRAY_BUFFER, heat_vert_VBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(VBOPosNormalColor)*1,
    &heat_vert[0], GL_STATIC_DRAW);
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

void Mesh::setupVBOs() {
  HandleGLError("enter setupVBOs");
  SetupMesh();
  SetupHeat();
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

  HandleGLError("enter draw mesh");
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(0));
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(sizeof(glm::vec3)));
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE,
    sizeof(VBOPosNormalColor), reinterpret_cast<void*>(sizeof(glm::vec3)*2));
  glDrawElements(GL_TRIANGLES, mesh_tri_indices.size()*3, GL_UNSIGNED_INT, 0);
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);
  HandleGLError("leaving draw mesh");

  HandleGLError("enter draw heat");
  glPointSize(20);
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
    for (unsigned int i = 0; i < 10; i++) {
      // Calculate heat from main heat source(s)
      for (unsigned int h = 0; h < heat_sources.size(); ++h) {
        HeatSource* hs = heat_sources[h];
        std::vector<float> distances;
        float distance_sum, distance_max = 0;
        for (unsigned int j = 0; j < vertices.size(); ++j) {
          Vertex* v = vertices[j];
          float distance = glm::distance(v->getPos(), hs->position);
          float distance_squared = distance*distance;
          distances.push_back(distance_squared);
          distance_sum += distance;
          if (distance > distance_max) {
            distance_max = distance;
          }
        }
        for (unsigned int j = 0; j < vertices.size(); ++j) {
          distances[j] = distance_max / distances[j] / distance_sum;
        }
        for (unsigned int j = 0; j < vertices.size(); ++j) {
          Vertex* v = vertices[j];
          float energy = hs->heat_output * distances[j] * args->timestep;
          v->addEnergy(energy);
        }
      }

      for (unsigned int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        v->receiveEnergy();
        // std::cout << v->getTemperature() << std::endl;
      }

      for (unsigned int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        glm::vec3 force(0.0, 0.0, 0.0);
        if (v->isMelting()) {
          force.y -= (v->getTemperature() - v->getMaterial()->melting_point) * v->getMaterial()->solidity / 100.0;
        }
        v->setForce(force);
      }

      for (edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
        Edge* edge = itr->second;
        if (edge->getTopVertex() == edge->getStartVertex()) {
          Vertex* top_vertex = edge->getTopVertex();
          Vertex* bottom_vertex = edge->getBottomVertex();
          if (top_vertex->isOnFloor()) continue;

          // If the spring is extended
          float current_distance = glm::distance(top_vertex->getPos(), bottom_vertex->getPos());
          float original_distance = edge->getOriginalDistance();
          if (current_distance > original_distance) {
            // Pull top vertex down a little
            glm::vec3 top_force = top_vertex->getForce();
            float k = 10.0;
            float deformation = (current_distance-original_distance)/original_distance;
            if (deformation > 1.0) deformation = 1.0;
            // if (deformation > )
            // std::cout << deformation << std::endl;
            top_force += k * deformation * glm::normalize(bottom_vertex->getPos()-top_vertex->getPos());
            // std::cout << "DOING IT" << std::endl;
            top_vertex->setForce(top_force);
          }
        }
      }

      // for (unsigned int j = 0; j < vertices.size(); ++j) {
      //   Vertex* v = vertices[j];
      //   glm::vec3 force_under(0.0, 0.0, 0.0);
      //   glm::vec3 pos_vxz = v->getPos();
      //   pos_vxz.y = 0;
      //   for (unsigned int k = 0; k < vertices.size(); ++k) {
      //     if (j == k) continue;
      //     Vertex* u = vertices[k];
      //     glm::vec3 pos_uxz = u->getPos();
      //     pos_uxz.y = 0;
      //     float distance = glm::distance(pos_vxz, pos_uxz);
      //     std::cout << distance << std::endl;
      //   }
      // }

      for (edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
        Edge* edge = itr->second;
        if (edge->getTopVertex() == edge->getEndVertex()) continue;
        Vertex* top_vertex = edge->getTopVertex();
        glm::vec3 force_under = edge->getBottomVertex()->getForce();
        Edge* next_edge = edges_by_top[edge->getBottomVertex()];
        int n = 25;
        while (next_edge != NULL && n > 0) {
          force_under += next_edge->getBottomVertex()->getForce();
          next_edge = edges_by_top[edge->getBottomVertex()];
          n--;
        }
        // print_vec(force_under);
        force_under.x = 0;
        force_under.z = 0;
        if (force_under.y > 0) {
          force_under.y = 0;
        }
        top_vertex->setForceUnder(force_under);
      }

      for (edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
        Edge* edge = itr->second;
        if (edge->getTopVertex() == edge->getEndVertex()) continue;
        Vertex* v = edge->getStartVertex();
        Vertex* u = edge->getEndVertex();
        if (v->getPos()[1] == floor_y && u->getPos()[1] == floor_y) {
          float current_distance = glm::distance(v->getPos(), u->getPos());
          float original_distance = edge->getOriginalDistance();
          float k = 5.0;
          float deformation = (current_distance-original_distance)/original_distance;
            if (deformation > 0.25) deformation = 0.25;
          glm::vec3 force_away = k * deformation * glm::normalize(v->getPos()-u->getPos());
          force_away.y = 0;
          v->setForce(v->getForce() + force_away);
          u->setForce(u->getForce() - force_away);
        }
      }

      for (unsigned int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        glm::vec3 force = v->getForce();
        if (v->isOnFloor()) {
        } else if (v->isMelting()) {
          // force += 0.05f * v->getForceUnder();
        } else {
          force += 0.25f * v->getForceUnder();
        }
        glm::vec3 acceleration = force / v->getMass();
        glm::vec3 velocity = v->getVelocity();
        if (v->isOnFloor()) {
          velocity.y = 0;
        }
        velocity += acceleration * args->timestep;
        glm::vec3 position = v->getPos();
        position += velocity * args->timestep;
        v->setVelocity(velocity);
        v->setPos(position);
      }

      for (edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
        Edge* edge = itr->second;
        if (edge->getTopVertex() == edge->getEndVertex()) continue;
        glm::vec3 bottom_position = edge->getBottomVertex()->getPos();
        glm::vec3 top_position = edge->getTopVertex()->getPos();
        if (bottom_position.y < floor_y) {
          bottom_position.y = floor_y + 0.0001;
          edge->getBottomVertex()->setPos(bottom_position);
          edge->getBottomVertex()->setOnFloor();
        }
        if (top_position.y < floor_y) {
          if (edge->getBottomVertex()->isOnFloor()) {
            top_position.y = bottom_position.y + 0.001;
          } else {
            top_position.y = floor_y + 0.0001;
          }
          edge->getTopVertex()->setPos(top_position);
          edge->getTopVertex()->setOnFloor();
        }
      }

      for (unsigned int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        glm::vec3 p = v->getPos();
        float size = 1.0;
        if (p.x < bbox.getMin().x * size) {
          p.x = bbox.getMin().x * size;
        }
        if (p.x > bbox.getMax().x * size) {
          p.x = bbox.getMax().x * size;
        }
        if (p.z < bbox.getMin().z * size) {
          p.z = bbox.getMin().z * size;
        }
        if (p.z > bbox.getMax().z * size) {
          p.z = bbox.getMax().z * size;
        }
        v->setPos(p);
      }
    }
  }
  setupVBOs();
}

void Mesh::animateHeat() {
  // float x = heat_position[0];
  // float y = heat_position[1];
  // float z = heat_position[2];
  // float angle = 0.05;
  // heat_position = glm::vec3(x * cos(angle) - z * sin(angle), y, x * sin(angle) + z * cos(angle));
}

#include "glCanvas.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cfloat>
#include <set>
#include <map>
#include <queue>

#include "mesh.h"
#include "edge.h"
#include "vertex.h"
#include "triangle.h"
#include "mtrand.h"

// to give a unique id number to each triangles
int Triangle::next_triangle_id = 0;

// =======================================================================
// MESH DESTRUCTOR 
// =======================================================================

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

// =======================================================================
// MODIFIERS:   ADD & REMOVE
// =======================================================================

Vertex* Mesh::addVertex(const glm::vec3 &position) {
  int index = numVertices();
  Vertex *v = new Vertex(index, position);
  vertices.push_back(v);
  if (numVertices() == 1)
    bbox = BoundingBox(position,position);
  else 
    bbox.Extend(position);
  return v;
}


void Mesh::addTriangle(Vertex *a, Vertex *b, Vertex *c) {
  // create the triangle
  Triangle *t = new Triangle();
  // create the edges
  Edge *ea = new Edge(a,b,t);
  Edge *eb = new Edge(b,c,t);
  Edge *ec = new Edge(c,a,t);
  // point the triangle to one of its edges
  t->setEdge(ea);
  // connect the edges to each other
  ea->setNext(eb);
  eb->setNext(ec);
  ec->setNext(ea);
  // verify these edges aren't already in the mesh 
  // (which would be a bug, or a non-manifold mesh)
  assert (edges.find(std::make_pair(a,b)) == edges.end());
  assert (edges.find(std::make_pair(b,c)) == edges.end());
  assert (edges.find(std::make_pair(c,a)) == edges.end());
  // add the edges to the master list
  edges[std::make_pair(a,b)] = ea;
  edges[std::make_pair(b,c)] = eb;
  edges[std::make_pair(c,a)] = ec;
  // connect up with opposite edges (if they exist)
  edgeshashtype::iterator ea_op = edges.find(std::make_pair(b,a)); 
  edgeshashtype::iterator eb_op = edges.find(std::make_pair(c,b)); 
  edgeshashtype::iterator ec_op = edges.find(std::make_pair(a,c)); 
  if (ea_op != edges.end()) { ea_op->second->setOpposite(ea); }
  if (eb_op != edges.end()) { eb_op->second->setOpposite(eb); }
  if (ec_op != edges.end()) { ec_op->second->setOpposite(ec); }
  // add the triangle to the master list
  assert (triangles.find(t->getID()) == triangles.end());
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
  edges.erase(std::make_pair(a,b)); 
  edges.erase(std::make_pair(b,c)); 
  edges.erase(std::make_pair(c,a));
  triangles.erase(t->getID());
  // clean up memory
  delete ea;
  delete eb;
  delete ec;
  delete t;
}


// =======================================================================
// Helper functions for accessing data in the hash table
// =======================================================================

Edge* Mesh::getMeshEdge(Vertex *a, Vertex *b) const {
  edgeshashtype::const_iterator iter = edges.find(std::make_pair(a,b));
  if (iter == edges.end()) return NULL;
  return iter->second;
}

Vertex* Mesh::getChildVertex(Vertex *p1, Vertex *p2) const {
  vphashtype::const_iterator iter = vertex_parents.find(std::make_pair(p1,p2)); 
  if (iter == vertex_parents.end()) return NULL;
  return iter->second; 
}

void Mesh::setParentsChild(Vertex *p1, Vertex *p2, Vertex *child) {
  assert (vertex_parents.find(std::make_pair(p1,p2)) == vertex_parents.end());
  vertex_parents[std::make_pair(p1,p2)] = child; 
}

Vertex* Mesh::getCenter(Edge* edge) {
  Vertex* vertex = getChildVertex(edge->getStartVertex(), edge->getEndVertex());
  if (vertex == NULL) {
    vertex = getChildVertex(edge->getEndVertex(), edge->getStartVertex());
  }
  return vertex;
}


// =======================================================================
// the load function parses very simple .obj files
// the basic format has been extended to allow the specification 
// of crease weights on the edges.
// =======================================================================

#define MAX_CHAR_PER_LINE 200

void Mesh::Load(const std::string &input_file) {

  std::ifstream istr(input_file.c_str());
  if (!istr) {
    std::cout << "ERROR! CANNOT OPEN: " << input_file << std::endl;
    return;
  }

  char line[MAX_CHAR_PER_LINE];
  std::string token, token2;
  float x,y,z;
  int a,b,c;
  int index = 0;
  int vert_count = 0;
  int vert_index = 1;

  // TODO: maybe this shouldn't be hardcoded
  heat_position = glm::vec3(0.1, 0.1, 0.1);

  // read in each line of the file
  while (istr.getline(line,MAX_CHAR_PER_LINE)) { 
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
      Vertex* newVertex = addVertex(glm::vec3(x,y,z));
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
        assert (a_ >= 0 && a_ < numVertices());
        assert (b_ >= 0 && b_ < numVertices());
        assert (c_ >= 0 && c_ < numVertices());
        addTriangle(getVertex(a_),getVertex(b_),getVertex(c_));
        b = c;
      }
    } else if (token == std::string("e")) {
      a = b = -1;
      ss >> a >> b >> token2;
      // whoops: inconsistent file format, don't subtract 1
      assert (a >= 0 && a <= numVertices());
      assert (b >= 0 && b <= numVertices());
      if (token2 == std::string("inf")) x = 1000000; // this is close to infinity...
      x = atof(token2.c_str());
      Vertex *va = getVertex(a);
      Vertex *vb = getVertex(b);
      Edge *ab = getMeshEdge(va,vb);
      Edge *ba = getMeshEdge(vb,va);
      assert (ab != NULL);
      assert (ba != NULL);
      ab->setCrease(x);
      ba->setCrease(x);
    } else if (token == std::string("vt")) {
    } else if (token == std::string("vn")) {
    } else if (token[0] == '#') {
    } else {
      printf ("LINE: '%s'",line);
    }
  }

  std::cout << "Loaded " << numTriangles() << " triangles." << std::endl;

  assert (numTriangles() > 0);
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
  glm::vec3 normal = glm::cross(v12,v23);
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
    return glm::vec3(0,1,0); 
  } else if (e->getOpposite() == NULL) {
    return glm::vec3(1,0,0); 
  } else if (e->getCrease() > 0) {
    return glm::vec3(1,1,0);
  } else {
    return glm::vec3(0,0,0.0);
  }
  return glm::vec3(0,0,0.0); 
}


void Mesh::TriVBOHelper( std::vector<glm::vec3> &indexed_verts,
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
  glm::vec3 center_color(1,1,1);
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

      Edge *starting = t->getEdge(); // ab
      Edge *current = starting;
      glm::vec3 normal_a(0, 0, 0);
      float triangles = 0.0f;

      while(triangles == 0 || ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_a += ComputeNormal((*triangle)[0]->getPos(), (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_a = normal_a / triangles;

      glm::vec3 normal_b(0, 0, 0);
      triangles = 0.0f;
      starting = t->getEdge()->getNext();
      current = starting;
      while(triangles == 0 || ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_b += ComputeNormal((*triangle)[0]->getPos(), (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_b = normal_b / triangles;

      glm::vec3 normal_c(0, 0, 0);
      triangles = 0.0f;
      starting = t->getEdge()->getNext()->getNext();
      current = starting;
      while(triangles == 0 || ((current != starting) && current != NULL)) {
        Triangle *triangle = current->getTriangle();
        normal_c += ComputeNormal((*triangle)[0]->getPos(), (*triangle)[1]->getPos(), (*triangle)[2]->getPos());
        current = current->getNext()->getNext()->getOpposite();
        triangles++;
      }
      normal_c = normal_c / triangles;


      TriVBOHelper(indexed_verts,mesh_tri_indices,
                   a,b,c,
                   normal_a,normal_b,normal_c,
                   edgecolor_ab,edgecolor_bc,edgecolor_ca);

      
    } else {
      // for flat shading, use the triangle normal at each vertex
      // use the normal of the triangl
      glm::vec3 normal = ComputeNormal(a,b,c);
      TriVBOHelper(indexed_verts,mesh_tri_indices,
                   a,b,c,
                   normal,normal,normal,
                   edgecolor_ab,edgecolor_bc,edgecolor_ca);

    }
  }
        
  // the vertex data
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glBufferData(GL_ARRAY_BUFFER, indexed_verts.size() * sizeof(glm::vec3), &indexed_verts[0], GL_STATIC_DRAW);
  // the index data (refers to vertex data)
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices.size() * sizeof(unsigned int), &mesh_tri_indices[0] , GL_STATIC_DRAW);

  num_mini_triangles = mesh_tri_indices.size();

  HandleGLError("leaving setupVBOs");
}


void Mesh::drawVBOs(const glm::mat4 &ProjectionMatrix,const glm::mat4 &ViewMatrix,const glm::mat4 &ModelMatrix) {
  HandleGLError("enter drawVBOs");

  // prepare data to send to the shaders
  glm::mat4 MVP = ProjectionMatrix * ViewMatrix * ModelMatrix;

  glm::vec3 lightPos = glm::vec3(4,4,4);
  glUniform3f(GLCanvas::LightID, lightPos.x, lightPos.y, lightPos.z);
  glUniformMatrix4fv(GLCanvas::MatrixID, 1, GL_FALSE, &MVP[0][0]);
  glUniformMatrix4fv(GLCanvas::ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
  glUniformMatrix4fv(GLCanvas::ViewMatrixID, 1, GL_FALSE, &ViewMatrix[0][0]);
  glUniform1i(GLCanvas::wireframeID, args->wireframe);

  // triangle vertex positions
  glBindBuffer(GL_ARRAY_BUFFER, mesh_tri_verts_VBO);
  glEnableVertexAttribArray(0);
  glVertexAttribPointer(0,                  // attribute
			3,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			3*sizeof(glm::vec3),// stride
			(void*)0            // array buffer offset
                        );
  // triangle vertex normals
  glEnableVertexAttribArray(1);
  glVertexAttribPointer(1,                      // attribute
			3,                      // size
			GL_FLOAT,               // type
			GL_FALSE,               // normalized?
			3*sizeof(glm::vec3),    // stride
			(void*)sizeof(glm::vec3)// array buffer offset
                        );
  // triangle vertex colors
  glEnableVertexAttribArray(2);
  glVertexAttribPointer(2,                          // attribute
			3,                          // size
			GL_FLOAT,                   // type
			GL_FALSE,                   // normalized?
			3*sizeof(glm::vec3),        // stride
			(void*)(sizeof(glm::vec3)*2)// array buffer offset
                        );
  // triangle indices
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_tri_indices_VBO);
  glDrawElements(GL_TRIANGLES,         // mode
                 num_mini_triangles*3, // count
                 GL_UNSIGNED_INT,      // type
                 (void*)0              // element array buffer offset
                 );
  glDisableVertexAttribArray(0);
  glDisableVertexAttribArray(1);
  glDisableVertexAttribArray(2);

  // =================================
  // draw the different types of edges
  //  if (args->wireframe) {

  HandleGLError("leaving drawVBOs");
}


void Mesh::cleanupVBOs() {
  glDeleteBuffers(1, &mesh_VAO);
  glDeleteBuffers(1, &mesh_tri_verts_VBO);
  glDeleteBuffers(1, &mesh_tri_indices_VBO);
}

// Animation function
void Mesh::animate() {
  if (args->animate) {
    // do 10 steps of animation before rendering
    for (int i = 0; i < 10; i++) {
      for (int j = 0; j < vertices.size(); ++j) {
        Vertex* v = vertices[j];
        glm::vec3 position = v->getPos();
        position += args->timestep * v->getVelocity();
        v->setPos(position);
      }
    }
  }
  setupVBOs();
}


// =================================================================
// SUBDIVISION
// =================================================================
void Mesh::LoopSubdivision() {
  printf ("Subdivide the mesh!\n");
  // =====================================
  // ASSIGNMENT: complete this functionality
  // =====================================
  std::vector<std::pair<Vertex*, Vertex*>> creased_edges;
  std::set<Vertex*> creased_vertices;
  for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    Vertex* start_vertex = edge->getStartVertex();
    Vertex* end_vertex = edge->getEndVertex();
    if (getChildVertex(start_vertex, end_vertex) == NULL && getChildVertex(end_vertex, start_vertex) == NULL) {
      Vertex* center_vertex;
      if(edge->getOpposite() == NULL || edge->getCrease() > 0) {
        center_vertex = addVertex(start_vertex->getPos() * 0.5f + end_vertex->getPos() * 0.5f);
        creased_vertices.insert(start_vertex);
        creased_vertices.insert(end_vertex);
        creased_vertices.insert(center_vertex);
        creased_edges.push_back(std::make_pair(start_vertex, center_vertex));
        creased_edges.push_back(std::make_pair(center_vertex, end_vertex));
      } else {
        Vertex* top_vertex = edge->getNext()->getEndVertex();
        Vertex* bottom_vertex = edge->getOpposite()->getNext()->getEndVertex();
        assert(start_vertex != end_vertex);
        assert(start_vertex != top_vertex);
        assert(start_vertex != bottom_vertex);
        assert(end_vertex != top_vertex);
        assert(end_vertex != bottom_vertex);
        assert(top_vertex != bottom_vertex);
        center_vertex = addVertex(start_vertex->getPos() * 0.375f + top_vertex->getPos() * 0.125f + end_vertex->getPos() * 0.375f + bottom_vertex->getPos() * 0.125f);
      }
      setParentsChild(start_vertex, end_vertex, center_vertex);
    }
  }

  std::set<Triangle*> triangles_to_remove;
  std::vector<Vertex**> triangles_to_add;
  for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    Triangle* triangle = edge->getTriangle();
    if (triangles_to_remove.count(triangle) == 0) {
      triangles_to_remove.insert(triangle);
      if(edge->getNext() == NULL || edge->getNext()->getNext() == NULL)
        continue;
      // Need to create four triangles
      Vertex** triangle_one = new Vertex*[3];
      triangle_one[0] = getCenter(edge->getNext()->getNext());
      triangle_one[1] = edge->getStartVertex();
      triangle_one[2] = getCenter(edge);
      triangles_to_add.push_back(triangle_one);

      Vertex** triangle_two = new Vertex*[3];
      triangle_two[0] = getCenter(edge);
      triangle_two[1] = getCenter(edge->getNext());
      triangle_two[2] = getCenter(edge->getNext()->getNext());
      triangles_to_add.push_back(triangle_two);

      Vertex** triangle_three = new Vertex*[3];
      triangle_three[0] = getCenter(edge);
      triangle_three[1] = edge->getEndVertex();
      triangle_three[2] = getCenter(edge->getNext());
      triangles_to_add.push_back(triangle_three);

      Vertex** triangle_four = new Vertex*[3];
      triangle_four[0] = getCenter(edge->getNext());
      triangle_four[1] = edge->getNext()->getEndVertex();
      triangle_four[2] = getCenter(edge->getNext()->getNext());
      triangles_to_add.push_back(triangle_four);
    }
  }

  for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    if(edge->getOpposite() == NULL || edge->getCrease() > 0) {
      Vertex* start_vertex = edge->getStartVertex(); // Center vertex
      glm::vec3 position = 0.75F * start_vertex->getPos();
      for (std::vector<Vertex*>::iterator vtr = vertices.begin(); vtr != vertices.end(); ++vtr) {
        Edge* possible_edge = getMeshEdge(*vtr, start_vertex);
        if (possible_edge != NULL && possible_edge->getOpposite() == NULL) {
          position += 0.125F * (*vtr)->getPos();
        }
      }
      position += 0.125F * edge->getEndVertex()->getPos();
      start_vertex->setPos(position);
    }
  }

  // Handle vertices with valence != 6
  std::map<Vertex*, std::vector<Edge*>> vertices_to_edges;
  for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    if (edge->getOpposite() != NULL && edge->getCrease() == 0) {
      Vertex* vertex = edge->getStartVertex();
      vertices_to_edges[vertex].push_back(edge);
    }
  }

  for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
    Edge* edge = itr->second;
    if (edge->getOpposite() == NULL || edge->getCrease() > 0) {
      Vertex* vertex = edge->getStartVertex();
      vertices_to_edges.erase(vertex);
      vertex = edge->getEndVertex();
      vertices_to_edges.erase(vertex);
    }
  }

  for(std::map<Vertex*, std::vector<Edge*>>::iterator itr = vertices_to_edges.begin(); itr != vertices_to_edges.end(); ++itr) {
    if (creased_vertices.count(itr->first) > 0) {
      continue;
    }
    int n = itr->second.size(), m = 0;
    float coeff = 3.0f / (8.0f * n);
    if (n == 3) {
      coeff = 3.0f / 16.0f;
    }
    glm::vec3 position;
    bool can_move = true;
    for(std::vector<Edge*>::iterator jtr = itr->second.begin(); jtr != itr->second.end(); ++jtr) {
      Edge* edge = *jtr;
      std::cout << edge->getCrease() << std::endl;
      if (edge->getEndVertex() != itr->first) {
        if (creased_vertices.count(edge->getEndVertex()) > 0)
          can_move = false;
        assert(edge->getStartVertex() == itr->first);
        position += (coeff * edge->getEndVertex()->getPos());
        m++;
      } else {
        if (creased_vertices.count(edge->getStartVertex()) > 0)
          can_move = false;
        assert(edge->getEndVertex() == itr->first);
        position += (coeff * edge->getStartVertex()->getPos());
        m++;
      }
    }
    assert(m == n);
    position += (1.0f - n * coeff) * itr->first->getPos();
    if (can_move)
      itr->first->setPos(position);
  }

  for(std::set<Triangle*>::iterator itr = triangles_to_remove.begin(); itr != triangles_to_remove.end(); ++itr) {
    removeTriangle(*itr);
  }

  for(std::vector<Vertex**>::iterator jtr = triangles_to_add.begin(); jtr != triangles_to_add.end(); ++jtr) {
    Vertex **vertices = *jtr;
    AttemptAddTriangle(vertices[0],vertices[1],vertices[2]); // 01, 12, 20
  }

  for(std::vector<std::pair<Vertex*, Vertex*>>::iterator ctr = creased_edges.begin(); ctr != creased_edges.end(); ++ctr) {
    Edge* tmp_edge;
    tmp_edge = getMeshEdge(ctr->first, ctr->second);
    if (tmp_edge != NULL) {
      tmp_edge->setCrease(2);
    }
    tmp_edge = getMeshEdge(ctr->second, ctr->first);
    if (tmp_edge != NULL) {
      tmp_edge->setCrease(2);
    }
  }
}


// =================================================================
// SIMPLIFICATION
// =================================================================

std::vector<Vertex*> Mesh::AddRemoveTriangles(Edge* starting, std::vector<Triangle*>& to_remove, std::vector<Vertex**>& to_add, Vertex* new_vertex, Triangle* left_triangle, Triangle* right_triangle, Vertex* left_vertex, Vertex* right_vertex) {
  Edge *current = starting;
  int triangles = 0;
  std::vector<Vertex*> removed_vertices;

  while(triangles == 0 || ((current != starting) && current != NULL)) {
    Triangle *triangle = current->getTriangle();
    if(triangle != left_triangle && triangle != right_triangle) {
      for(int i = 0; i < 3; ++i) {
        if ((*triangle)[i] == left_vertex || (*triangle)[i] == right_vertex) {
          Vertex **new_triangle = new Vertex*[3];
          new_triangle[i] = new_vertex;
          new_triangle[(i + 1) % 3] = (*triangle)[(i + 1) % 3];
          removed_vertices.push_back((*triangle)[(i + 1) % 3]);
          new_triangle[(i + 2) % 3] = (*triangle)[(i + 2) % 3];
          removed_vertices.push_back((*triangle)[(i + 2) % 3]);
          to_add.push_back(new_triangle);
          break;
        }
      }
      to_remove.push_back(triangle);
    }
    triangles++;
    current = current->getNext()->getNext()->getOpposite();
  }
  // std::cout << "FOUND TRIS: " << triangles << std::endl;
  // if (triangles >= 8) return false;
  return removed_vertices;
}

void Mesh::Simplification(int target_tri_count) {
  // clear out any previous relationships between vertices
  vertex_parents.clear();

  printf ("Simplify the mesh! %d -> %d\n", numTriangles(), target_tri_count);

  // =====================================
  // ASSIGNMENT: complete this functionality
  // =====================================

  std::set<Edge*> edges_to_skip;

  while(numTriangles() > target_tri_count) {
    std::cout << numTriangles() << std::endl;
    float shortest_edge = FLT_MAX;
    Edge *edge;
    for(edgeshashtype::iterator itr = edges.begin(); itr != edges.end(); ++itr) {
      if (itr->second->Length() < shortest_edge && itr->second->getOpposite() != NULL && edges_to_skip.count(itr->second) == 0 && itr->second->getOpposite() != NULL) {
        shortest_edge = itr->second->Length();
        edge = itr->second;
      }
    }

    Vertex *left_vertex = edge->getStartVertex();
    Vertex *right_vertex = edge->getEndVertex();
    glm::vec3 position_sums = left_vertex->getPos() + right_vertex->getPos();
    glm::vec3 position_avg = glm::vec3(position_sums[0] / 2.0, position_sums[1] / 2.0, position_sums[2] / 2.0);
    Vertex *new_vertex = addVertex(position_avg);

    Triangle *left_triangle = edge->getTriangle();
    Triangle *right_triangle = edge->getOpposite()->getTriangle();

    std::vector<Triangle*> to_remove;
    to_remove.push_back(left_triangle);
    to_remove.push_back(right_triangle);
    std::vector<Vertex**> to_add;

    std::vector<Vertex*> left_vertices = AddRemoveTriangles(edge, to_remove, to_add, new_vertex, left_triangle, right_triangle, left_vertex, right_vertex);
    std::vector<Vertex*> right_vertices = AddRemoveTriangles(edge->getOpposite(), to_remove, to_add, new_vertex, left_triangle, right_triangle, left_vertex, right_vertex);
    int in_common = 0;
    for(std::vector<Vertex*>::iterator ltr = left_vertices.begin(); ltr != left_vertices.end(); ++ltr) {
      for(std::vector<Vertex*>::iterator rtr = right_vertices.begin(); rtr != right_vertices.end(); ++rtr) {
        if(*ltr == *rtr) {
          in_common++;
        }
      }
    }
    if (in_common > 2) {
      std::cout << "TOO MANY IN COMMON" << std::endl;
      edges_to_skip.insert(edge);
      continue;
    }

    // std::cout << "REMOVING COUNT: " << to_remove.size() << std::endl;
    // std::cout << "ADDING COUNT: " << to_add.size() << std::endl;

    for(std::vector<Triangle*>::iterator jtr = to_remove.begin(); jtr != to_remove.end(); ++jtr) {
      Triangle* triangle = *jtr;
      removeTriangle(triangle);
    }

    for(std::vector<Vertex**>::iterator jtr = to_add.begin(); jtr != to_add.end(); ++jtr) {
      Vertex **vertices = *jtr;
      AttemptAddTriangle(vertices[1],vertices[0],vertices[2]); // 10, 02, 21
      AttemptAddTriangle(vertices[0],vertices[1],vertices[2]); // 01, 12, 20
    }
  }
}

void Mesh::AttemptAddTriangle(Vertex* a, Vertex* b, Vertex* c) {
  if (a != NULL && b != NULL && c != NULL && getMeshEdge(a,b) == NULL
    && getMeshEdge(b,c) == NULL && getMeshEdge(c,a) == NULL) {
      addTriangle(a,b,c);
    } else {
      // std::cout << "CAN'T ADD" << std::endl;
    }
}


// =================================================================

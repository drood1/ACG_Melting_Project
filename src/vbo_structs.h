#ifndef SRC_VBO_STRUCTS_H_
#define SRC_VBO_STRUCTS_H_


struct VBOPosNormalColor {
  VBOPosNormalColor(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) {
    x = p.x; y = p.y; z = p.z;
    nx = n.x; ny = n.y; nz = n.z;
    r  =  c.x;  g =  c.y;  b =  c.z; a = c.z;
  }

  VBOPosNormalColor(const glm::vec3 &p, const glm::vec3 &n,
    const glm::vec4 &c, const glm::vec4 &wc, float s_, float t_) {
    x = p.x; y = p.y; z = p.z;
    nx = n.x; ny = n.y; nz = n.z;
    r  =  c.x;  g =  c.y;  b =  c.z;  a = c.a;
  }

  float x, y, z;         // position
  float nx, ny, nz;      // normal
  float r, g, b, a;      // color
};

struct VBOIndexedTri {
  VBOIndexedTri() {}
  VBOIndexedTri(unsigned int a, unsigned int b, unsigned int c) {
    verts[0] = a;
    verts[1] = b;
    verts[2] = c;
  }
  unsigned int verts[3];
};

#endif  // SRC_VBO_STRUCTS_H_

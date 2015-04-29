#ifndef SRC_GLCANVAS_H_
#define SRC_GLCANVAS_H_

// Graphics Library Includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <string>

class ArgParser;
class Mesh;
class Camera;

// ====================================================================
// NOTE:  All the methods and variables of this class are static
// ====================================================================

// ======================================================================
// helper structures for VBOs, for rendering (note, the data stored in
// each of these is application specific, adjust as needed!)
struct VBOPosNormalColor {
  VBOPosNormalColor(const glm::vec3 &p, const glm::vec3 &n, const glm::vec4 &c) {
    x = p.x; y = p.y; z = p.z;
    nx = n.x; ny = n.y; nz = n.z;
    r  =  c.x;  g =  c.y;  b =  c.z; a = c.z;
//    wr = 1; wg = 1; wb = 1; wa = 1;
//    s = 0;
//    t = 0;
  }

  VBOPosNormalColor(const glm::vec3 &p, const glm::vec3 &n,
    const glm::vec4 &c, const glm::vec4 &wc, float s_, float t_) {
    x = p.x; y = p.y; z = p.z;
    nx = n.x; ny = n.y; nz = n.z;
    r  =  c.x;  g =  c.y;  b =  c.z;  a = c.a;
//    wr = wc.x; wg = wc.y; wb = wc.z; wa = wc.a;
//    s = s_;
//    t = t_;
  }

  float x, y, z;         // position
  float nx, ny, nz;      // normal
  float r, g, b, a;      // color
//  float wr, wg, wb, wa;  // wireframe color
//  float s, t;            // texture coordinates
};

class GLCanvas {
 public:
  // various static variables
  static ArgParser *args;
  static Mesh *mesh;
  static Camera* camera;
  static GLFWwindow* window;

  static GLuint ViewMatrixID;
  static GLuint ModelMatrixID;
  static GLuint LightID;
  static GLuint MatrixID;
  static GLuint programID;
  static GLuint wireframeID;

  // mouse position
  static int mouseX;
  static int mouseY;
  // which mouse button
  static bool leftMousePressed;
  static bool middleMousePressed;
  static bool rightMousePressed;
  // current state of modifier keys
  static bool shiftKeyPressed;
  static bool controlKeyPressed;
  static bool altKeyPressed;
  static bool superKeyPressed;

  static void initialize(ArgParser *_args, Mesh *_mesh);

  // Callback functions for mouse and keyboard events
  static void mousebuttonCB(GLFWwindow *window, int which_button, int action, int mods);
  static void mousemotionCB(GLFWwindow *window, double x, double y);
  static void keyboardCB(GLFWwindow *window, int key, int scancode, int action, int mods);
  static void error_callback(int error, const char* description);
};

// ====================================================================

// helper functions
GLuint LoadShaders(const std::string &vertex_file_path, const std::string &fragment_file_path);
std::string WhichGLError(GLenum &error);
int HandleGLError(const std::string &message = "", bool ignore = false);

#endif  // SRC_GLCANVAS_H_

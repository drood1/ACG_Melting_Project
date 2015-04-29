// Graphics Library Includes
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <cstdlib>

// for sleep
#if defined(_WIN32)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include "./argparser.h"
#include "./mesh.h"
#include "./glCanvas.h"
#include "./camera.h"

int main(int argc, char *argv[]) {
  // parse the command line arguments
  ArgParser args(argc, argv);
  // create & initialize the mesh object
  Mesh mesh(&args);
  GLCanvas::initialize(&args, &mesh);

  glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glDisable(GL_CULL_FACE);

  // transform the object as necessary to fit in the
  // (-1,-1,-1)->(1, 1, 1) box
  const BoundingBox &bb = mesh.getBoundingBox();
  glm::vec3 center;
  bb.getCenter(center);
  glm::mat4 myTranslateMatrix = glm::translate(-center);
  double maxDim = bb.maxDim();
  float scaleFactor = 2.0 / static_cast<float>(maxDim);
  glm::mat4 myScalingMatrix = glm::scale(glm::vec3(scaleFactor, scaleFactor, scaleFactor));
  glm::mat4 ModelMatrix = myScalingMatrix*myTranslateMatrix;

  while (!glfwWindowShouldClose(GLCanvas::window)) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glUseProgram(GLCanvas::programID);

    GLCanvas::camera->glPlaceCamera();

    // Build the matrix to position the camera based on keyboard and mouse input
    glm::mat4 ProjectionMatrix = GLCanvas::camera->getProjectionMatrix();
    glm::mat4 ViewMatrix = GLCanvas::camera->getViewMatrix();

    mesh.drawVBOs(ProjectionMatrix, ViewMatrix, ModelMatrix);
    mesh.animate();

    // Swap buffers
    glfwSwapBuffers(GLCanvas::window);
    glfwPollEvents();

    #if defined(_WIN32)
      Sleep(100);
    #else
      usleep(100);
    #endif
  }

  mesh.cleanupVBOs();
  glDeleteProgram(GLCanvas::programID);

  // Close OpenGL window and terminate GLFW
  glfwDestroyWindow(GLCanvas::window);
  glfwTerminate();
  exit(EXIT_SUCCESS);
}

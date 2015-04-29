#ifndef SRC_ARGPARSER_H_
#define SRC_ARGPARSER_H_

#include <string>
#include <cassert>
#include "./mtrand.h"


inline void separatePathAndFile(const std::string &input, std::string &path, std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;
  while (1) {
    int next = input.find('/', last+1);
    if (next != static_cast<int>(std::string::npos)) {
      last = next;
      continue;
    }
    next = input.find('\\', last+1);
    if (next != static_cast<int>(std::string::npos)) {
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last+1, input.size()-last-1);
    path = input.substr(0, last);
  }
}

class ArgParser {
 public:
  ArgParser() { DefaultValues(); }

  ArgParser(int argc, char *argv[]) {
    DefaultValues();
    for (int i = 1; i < argc; i++) {
      // TODO(austin): add timestep here
      if (std::string(argv[i]) == std::string("-input") ||
          std::string(argv[i]) == std::string("-i")) {
        i++; assert(i < argc);
        separatePathAndFile(argv[i], path, input_file);
      } else if (argv[i] == std::string("-size")) {
        i++; assert(i < argc);
        width = height = atoi(argv[i]);
      } else if (argv[i] == std::string("-wireframe")) {
        wireframe = 1;
      } else if (argv[i] == std::string("-gouraud")) {
        gouraud = true;
      } else {
        printf("whoops error with command line argument %d: '%s'\n", i, argv[i]);
        assert(0);
      }
    }
  }

  void DefaultValues() {
    width = 500;
    height = 500;
    wireframe = 0;
    gouraud = false;
    animate = false;
    timestep = 0.001;
  }

  // ==============
  // REPRESENTATION
  // all public! (no accessors)
  std::string input_file;
  std::string path;
  int width;
  int height;
  GLint wireframe;
  bool gouraud;
  MTRand mtrand;
  bool animate;
  float timestep;
};

#endif  // SRC_ARGPARSER_H_

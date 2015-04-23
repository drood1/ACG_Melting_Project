#include "glCanvas.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "cloth.h"
#include "argparser.h"
#include "utils.h"

// ================================================================================
// ================================================================================

Cloth::Cloth(ArgParser *_args) {
  args =_args;

  // open the file
  std::ifstream istr(std::string(args->path+'/'+args->cloth_file).c_str());
  assert (istr.good());
  std::string token;

  // read in the simulation parameters
  istr >> token >> k_structural; assert (token == "k_structural");  // (units == N/m)  (N = kg*m/s^2)
  istr >> token >> k_shear; assert (token == "k_shear");
  istr >> token >> k_bend; assert (token == "k_bend");
  istr >> token >> damping; assert (token == "damping");
  // NOTE: correction factor == .1, means springs shouldn't stretch more than 10%
  //       correction factor == 100, means don't do any correction
  istr >> token >> provot_structural_correction; assert (token == "provot_structural_correction");
  istr >> token >> provot_shear_correction; assert (token == "provot_shear_correction");

  // the cloth dimensions
  istr >> token >> nx >> ny; 
  assert (token == "m");
  assert (nx >= 2 && ny >= 2);

  // the corners of the cloth
  // (units == meters)
  glm::vec3 a,b,c,d;
  istr >> token >> a.x >> a.y >> a.z; assert (token == "p");
  istr >> token >> b.x >> b.y >> b.z; assert (token == "p");
  istr >> token >> c.x >> c.y >> c.z; assert (token == "p");
  istr >> token >> d.x >> d.y >> d.z; assert (token == "p");

  // fabric weight  (units == kg/m^2)
  // denim ~300 g/m^2
  // silk ~70 g/m^2
  double fabric_weight;
  istr >> token >> fabric_weight; assert (token == "fabric_weight");
  double area = AreaOfTriangle(a,b,c) + AreaOfTriangle(a,c,d);

  // create the particles
  particles = new ClothParticle[nx*ny];
  double mass = area*fabric_weight / double(nx*ny);
  for (int i = 0; i < nx; i++) {
    double x = i/double(nx-1);
    glm::vec3 ab = float(1-x)*a + float(x)*b;
    glm::vec3 dc = float(1-x)*d + float(x)*c;
    for (int j = 0; j < ny; j++) {
      double y = j/double(ny-1);
      ClothParticle &p = getParticle(i,j);
      glm::vec3 abdc = float(1-y)*ab + float(y)*dc;
      p.setOriginalPosition(abdc);
      p.setPosition(abdc);
      p.setVelocity(glm::vec3(0,0,0));
      p.setMass(mass);
      p.setFixed(false);
    }
  }

  // the fixed particles
  while (istr >> token) {
    assert (token == "f");
    int i,j;
    double x,y,z;
    istr >> i >> j >> x >> y >> z;
    ClothParticle &p = getParticle(i,j);
    p.setPosition(glm::vec3(x,y,z));
    p.setFixed(true);
  }

  computeBoundingBox();
}

// ================================================================================

void Cloth::computeBoundingBox() {
  box = BoundingBox(getParticle(0,0).getPosition());
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      box.Extend(getParticle(i,j).getPosition());
      box.Extend(getParticle(i,j).getOriginalPosition());
    }
  }
}

// ================================================================================



void Cloth::Animate() {


  // *********************************************************************  
  // ASSIGNMENT:
  //
  // Compute the forces on each particle, and update the state
  // (position & velocity) of each particle.
  //
  // Also, this is where you'll put the Provot correction for super-elasticity
  //
  // *********************************************************************    

  // std::cout << "Animate()" << std::endl;

  for(int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      ClothParticle& particle = getParticle(i, j);
      if (particle.isFixed()) continue;
      // glm::vec3 gravity_force = args->gravity * (float) particle.getMass();
      glm::vec3 spring_force;

      // Structural spring
      if (i > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i-1, j), (float) k_structural);
      }
      if (i+1 < nx) {
        spring_force += CalculateSpringForce(particle, getParticle(i+1, j), (float) k_structural);
      }
      if (j > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i, j-1), (float) k_structural);
      }
      if (j+1 < ny) {
        spring_force += CalculateSpringForce(particle, getParticle(i, j+1), (float) k_structural);
      }

      // Shear spring
      if (i > 0 && j > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i-1, j-1), (float) k_shear);
      }
      if (i+1 < nx && j > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i+1, j-1), (float) k_shear);
      }
      if (i > 0 && j+1 < ny) {
        spring_force += CalculateSpringForce(particle, getParticle(i-1, j+1), (float) k_shear);
      }
      if (i+1 < nx && j+1 < ny) {
        spring_force += CalculateSpringForce(particle, getParticle(i+1, j+1), (float) k_shear);
      }

      // Bend spring
      if (i-1 > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i-2, j), (float) k_bend);
      }
      if (i+2 < nx) {
        spring_force += CalculateSpringForce(particle, getParticle(i+2, j), (float) k_bend);
      }
      if (j-1 > 0) {
        spring_force += CalculateSpringForce(particle, getParticle(i, j-2), (float) k_bend);
      }
      if (j+2 < ny) {
        spring_force += CalculateSpringForce(particle, getParticle(i, j+2), (float) k_bend);
      }

      glm::vec3 acceleration = args->gravity + (spring_force / (float) particle.getMass());
      particle.setAcceleration(acceleration);

      glm::vec3 velocity = particle.getVelocity();
      velocity = (float) damping * 10.0f * (velocity + acceleration * (float) args->timestep);
      particle.setVelocity(velocity);
    }
  }

  for(int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      ClothParticle& particle = getParticle(i, j);
      if (particle.isFixed()) continue;

      glm::vec3 velocity = particle.getVelocity();
      glm::vec3 position = particle.getPosition();
      position += velocity * (float) args->timestep;
      particle.setPosition(position);

      // Structural spring
      if (i > 0) {
        ProvotCorrection(particle, getParticle(i-1, j), (float) provot_structural_correction);
      }
      if (i+1 < nx) {
        ProvotCorrection(particle, getParticle(i+1, j), (float) provot_structural_correction);
      }
      if (j > 0) {
        ProvotCorrection(particle, getParticle(i, j-1), (float) provot_structural_correction);
      }
      if (j+1 < ny) {
        ProvotCorrection(particle, getParticle(i, j+1), (float) provot_structural_correction);
      }

      // Shear spring
      if (i > 0 && j > 0) {
        ProvotCorrection(particle, getParticle(i-1, j-1), (float) provot_shear_correction);
      }
      if (i+1 < nx && j > 0) {
        ProvotCorrection(particle, getParticle(i+1, j-1), (float) provot_shear_correction);
      }
      if (i > 0 && j+1 < ny) {
        ProvotCorrection(particle, getParticle(i-1, j+1), (float) provot_shear_correction);
      }
      if (i+1 < nx && j+1 < ny) {
        ProvotCorrection(particle, getParticle(i+1, j+1), (float) provot_shear_correction);
      }
    }
  }

  // float max_accel = 0;
  double threshold = 100;
  float mean_acceleration = 0;
  // Check for instability
  for(int i = 0; i < nx; ++i) {
    for(int j = 0; j < ny; ++j) {
      ClothParticle& particle = getParticle(i, j);
      if (particle.isFixed()) continue;
      glm::vec3 acceleration = particle.getAcceleration();
      mean_acceleration += glm::length(acceleration);
    }
  }
  mean_acceleration /= nx * ny;

  if (mean_acceleration > threshold) {
    for(int i = 0; i < nx; ++i) {
      for(int j = 0; j < ny; ++j) {
        ClothParticle& particle = getParticle(i, j);
        if (particle.isFixed()) continue;
        particle.setPosition(particle.getOriginalPosition());
        glm::vec3 velocity;
        particle.setVelocity(velocity);
        glm::vec3 acceleration;
        particle.setAcceleration(acceleration);
      }
    }
    // Half timestep
    args->timestep /= 2.0;
    // std::cout << "Halfing, timestep=" << args->timestep << std::endl;
    // std::cout << "mean a=" << mean_acceleration << std::endl;
    return;
  }

  // commented out because an animated bounding box can be weird
  // computeBoundingBox();
}


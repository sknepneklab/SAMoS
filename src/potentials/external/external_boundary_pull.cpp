/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file external_boundary_pull.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Sep-2016
 * \brief Implementation of the self propulsion for active particles.
 */ 

#include "external_boundary_pull.hpp"

/*! Apply pulling force onto each boundary particle. */
void ExternalBoundaryPull::compute()
{
  int N = m_system->size();
  double alpha = m_alpha;

  // First we compute centre of mass to be albe to determine outward direction
  double xcm = 0.0, ycm = 0.0, zcm = 0.0;
  int cnt = 0;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.in_tissue)
    {
      xcm += pi.x;
      ycm += pi.y;
      zcm += pi.z;
      cnt++;
    }
  }
  xcm /= cnt;  ycm /= cnt;  zcm /= cnt;
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.boundary)
    {
      double force_sign = 1.0;
      double X = pi.x - xcm,  Y = pi.y - ycm,  Z = pi.z - zcm;
      double len_R = sqrt(X*X + Y*Y + Z*Z);
      X /= len_R;  Y /= len_R;  Y /= len_R;
      Particle& pj = m_system->get_particle(pi.boundary_neigh[0]);
      Particle& pk = m_system->get_particle(pi.boundary_neigh[1]);
      double xji = pj.x - pi.x, yji = pj.y - pi.y, zji = pj.z - pi.z;
      double len_ji = sqrt(xji*xji + yji*yji + zji*zji);
      xji /= len_ji;  yji /= len_ji;  zji /= len_ji; 
      double xki = pk.x - pi.x, yki = pk.y - pi.y, zki = pk.z - pi.z;
      double len_ki = sqrt(xki*xki + yki*yki + zki*zki);
      xki /= len_ki;  yki /= len_ki;  zki /= len_ki;
      double x = -(xji+xki), y = -(yji+yki), z = -(zji+zki);
      double len_r = sqrt(x*x + y*y + z*z);
      x /= len_r;  y /= len_r;  z /= len_r;
      // compute dot product with radius the vector connecting center of mass and pi
      if ((x*X + y*Y + z*Z) < 0.0)
        force_sign = -1.0;
      double factor = force_sign*alpha;
      pi.fx += factor*x;
      pi.fy += factor*y;
      pi.fz += factor*z;
    }
  }
}

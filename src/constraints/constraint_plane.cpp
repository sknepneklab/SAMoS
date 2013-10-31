/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_plane.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Implementation of the planar constraint
 */ 

/*! Force all particles to be confined to the surface of the xy plane 
 *  by setting z = 0, vz = 0, fz = 0
 */
void ConstraintPlane::enforce()
{
  bool periodic = m_system->get_periodic();
  int N = m_system->size();
  double xlo = -0.5*m_lx, xhi = 0.5*m_lx;
  double ylo = -0.5*m_ly, yhi = 0.5*m_ly;
  double zlo = -0.5*m_lz, zhi = 0.5*m_lz;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    p.z = 0.0;
    p.vz = 0.0;
    p.fz = 0.0;
    // Check periodic boundary conditions 
    if (periodic)
    {
      if (p.x < xlo) p.x += m_lx;
      else if (p.x > xhi) p.x -= m_lx;
      if (p.y < ylo) p.y += m_ly;
      else if (p.y > yhi) p.y -= m_ly;
      if (p.z < zlo) p.z += m_lz;
      else if (p.z > zhi) p.z -= m_lz;
    }
  }
}
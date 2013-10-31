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
 * \file constraint_sphere.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Implementation of the spherical constraint
 */ 

/*! Force all particles to be confined to the surface of the sphere and
 *  all velocities to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down using simple radial projection.
 */
void ConstraintSphere::enforce()
{
  int N = m_system->size();
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    double x = p.x, y = p.y, z = p.z;
    double R = sqrt(x*x + y*y + z*z);
    double s = m_r/R;
    // Scale back to the surface
    p.x *= s; p.y *= s; p.z *= s;
    // Compute unit normal
    double nx = p.x/m_r, ny = p.y/m_r, nz = p.z/m_r;
    // Compute tangent component of the velocity
    p.vx -= p.vx*nx; p.vy -= p.vy*ny; p.vz -= p.vz*nz;
  }
}
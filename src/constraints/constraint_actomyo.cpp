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
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file constraint_actomyo.cpp
 * \author Amit Das, dosamit@gmail.com
 * \date 20-Nov-2013
 * \brief Implementation of the actomyo constraint
 */ 

#include "constraint_actomyo.hpp"

/*! Force all particles of actin to be confined on a plane to a surface parallel to the xy direction and all particles of
 *   myosin to move at another plane of same size above the former containing the actins 
 *  by setting z = 0, vz = 0, fz = 0 for type 1 particles (actin)
 *  by setting z = z0, vz = 0, fz = 0 for type 2 particles (myosin backbone)
 *  by setting z = z0/2, for type 3 particles (mysoin sidechains); z0 ~ 2 x Morse minimum
 *  \param p particle to project onto the sphere
 */
void ConstraintActomyo::enforce(Particle& p)
{
  bool periodic = m_system->get_periodic();
  double xlo = -0.5*m_lx, xhi = 0.5*m_lx;
  double ylo = -0.5*m_ly, yhi = 0.5*m_ly;
  
  if (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end())
  {
    // NOTE: Currently, types are hard-coded. This breaks generality of the code. Needs to be rethought.
    if (p.get_type() == 1)  
    {  
      p.z = 0.0;
      p.vz = 0.0;
      p.fz = 0.0;
    }  

    if (p.get_type() == 2)
    {  
      p.z = 0.5*m_lz;
      p.vz = 0.0;
      p.fz = 0.0;
    }  

    if (p.get_type() == 3)
    {  
      p.z = 0.25*m_lz;
    }  

    // Check periodic boundary conditions 
    if (periodic)
    {
      if (p.x <= xlo) p.x += m_lx;
      else if (p.x >= xhi) p.x -= m_lx;
      if (p.y <= ylo) p.y += m_ly;
      else if (p.y >= yhi) p.y -= m_ly;
    }
    else // reflective boundary conditions
    {
      if (p.x < xlo) 
      {
        p.x = xlo;
        p.vx = -p.vx;
      }
      else if (p.x > xhi)
      {
        p.x = xhi;
        p.vx = -p.vx;
      }
      if (p.y < ylo) 
      {
        p.y = ylo;
        p.vy = -p.vy;
      }
      else if (p.y > yhi)
      {
        p.y = yhi;
        p.vy = -p.vy;
      }
    }
  }
}

/*! Rotate director of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintActomyo::rotate_director(Particle& p, double phi)
{
  if (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end())
  {
    // Sine and cosine of the rotation angle
    double c = cos(phi), s = sin(phi);
    // Rotation matrix around z axis
    double Rxx = c, Rxy = -s;
    double Ryx = s, Ryy = c;
    // Apply rotation matrix
    double nx = Rxx*p.nx + Rxy*p.ny;
    double ny = Ryx*p.nx + Ryy*p.ny;
    double len = sqrt(nx*nx + ny*ny);
    // Update particle director
    p.nx = nx/len;
    p.ny = ny/len;
  }
}

/*! Rotate velocity of a particle around the normal vector (z axis)
 *  \note This function assumes that the particle has already been
 *  projected onto the plane and that its director is also in plane
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintActomyo::rotate_velocity(Particle& p, double phi)
{
  if (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end())
  {
    // Sine and cosine of the rotation angle
    double c = cos(phi), s = sin(phi);
    // Rotation matrix around z axis
    double Rxx = c, Rxy = -s;
    double Ryx = s, Ryy = c;
    // Apply rotation matrix
    double vx = Rxx*p.vx + Rxy*p.vy;
    double vy = Ryx*p.vx + Ryy*p.vy;
    // Update particle director
    p.vx = vx;
    p.vy = vy;
  }
}

/*! Project particle torque onto the normal vector (z axis). The assumption here is that 
 *  the particle's director and velocity are all already in the xy plane
 *  and that it is constrained to the plane.
 *  \param p Particle whose torque to project
*/ 
double ConstraintActomyo::project_torque(Particle& p)
{
  if (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end())
    return p.tau_z;  
  else
    return 0.0;
}

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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file constraint_hourglass.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Feb-2015
 * \brief Implementation of the hourglass constraint
 */ 

#include "constraint_hourglass.hpp"

static void compute_gradient(double x, double y, double z, double R, double A, double n, double L, double& gx, double& gy, double& gz)
{
  gx = 2.0*x;
  gy = 2.0*y; 
  gz = (-4*A*n*M_PI*cos((2*n*M_PI*z)/L)*(R + A*sin((2*n*M_PI*z)/L)))/L;
}

/*! Compute normal to the surface at p
*/
void ConstraintHourglass::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lz = box->Lz;
  // Compute unit gradient of the implicit function (normal will be normalized gradient)
  compute_gradient(x,y,z,m_R,m_A,m_n,lz,Nx,Ny,Nz);
  // Normalize N
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}


/*! Force particle to be confined to the surface of a hourglass and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface of the hourglass
 *  \param p Particle which is to be projected onto the hourglass
 */
void ConstraintHourglass::enforce(Particle& p)
{
  // Project particle back onto the surface
  double x = p.x, y = p.y, z = p.z;
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double lz = box->Lz;
  // Compute reference gradient (SHAKE method)
  double ref_gx, ref_gy, ref_gz;
  compute_gradient(x,y,z,m_R,m_A,m_n,lz,ref_gx, ref_gy, ref_gz);
    
  bool satisfied;
  int iter = 0;
  while (iter++ < m_max_iter)
  {
    satisfied = true;
    x = p.x; y = p.y; z = p.z;
    // Compute unit gradient of the implicit function
    double gx, gy, gz;
    compute_gradient(x,y,z,m_R,m_A,m_n,lz,gx,gy,gz);
    double s = gx*ref_gx + gy*ref_gy + gz*ref_gz;
    double zfact = m_R + m_A*sin(2.0*M_PI/lz*m_n*z);
    double g = x*x + y*y - zfact*zfact;
    double lambda = g/s;
    if (fabs(g) > m_tol)
      satisfied = false;
    if (satisfied) break;
    p.x -= lambda*ref_gx;
    p.y -= lambda*ref_gy;
    p.z -= lambda*ref_gz;
  }
    
  double Nx, Ny, Nz;
  this->compute_normal(p,Nx,Ny,Nz);
  // compute v.N
  double v_dot_N = p.vx*Nx + p.vy*Ny + p.vz*Nz;
  // compute n.N
  double n_dot_N = p.nx*Nx + p.ny*Ny + p.nz*Nz;
  // Project velocity onto tangent plane
  p.vx -= v_dot_N*Nx; p.vy -= v_dot_N*Ny; p.vz -= v_dot_N*Nz;
  // Project director onto tangent plane
  p.nx -= n_dot_N*Nx; p.ny -= n_dot_N*Ny; p.nz -= n_dot_N*Nz;
  // normalize director
  double inv_len = 1.0/sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
  p.nx *= inv_len;  p.ny *= inv_len;  p.nz *= inv_len;
  if (periodic)
  {
    if (p.z > box->zhi) p.z -= box->Lz;
    else if (p.z < box->zlo) p.z += box->Lz;
  }
}

/*! Rotate director of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto hourglass and that its director is laying in 
 *  the tangent plane.
 *  \param p Particle whose director to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintHourglass::rotate_director(Particle& p, double phi)
{
  double U, V, W;
  this->compute_normal(p,U,V,W);
  // Compute angle sins and cosines
  double c = cos(phi), s = sin(phi);
  // Compute new velocity coordinates
  double nx = U*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nx*c + (-W*p.ny+V*p.nz)*s;
  double ny = V*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.ny*c + ( W*p.nx-U*p.nz)*s;
  double nz = W*(U*p.nx+V*p.ny+W*p.nz)*(1.0-c) + p.nz*c + (-V*p.nx+U*p.ny)*s;
  double len = sqrt(nx*nx + ny*ny + nz*nz);
  // Update particle director (normalize it along the way to collect for any numerical drift that may have occurred)
  p.nx = nx/len;
  p.ny = ny/len;
  p.nz = nz/len;  
}

/*! Rotate velocity of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto hourglass and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintHourglass::rotate_velocity(Particle& p, double phi)
{
  double U, V, W;
  this->compute_normal(p,U,V,W);
  // Compute angle sins and cosines
  double c = cos(phi), s = sin(phi);
  // Compute new velocity coordinates
  double vx = U*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vx*c + (-W*p.vy+V*p.vz)*s;
  double vy = V*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vy*c + ( W*p.vx-U*p.vz)*s;
  double vz = W*(U*p.vx+V*p.vy+W*p.vz)*(1.0-c) + p.vz*c + (-V*p.vx+U*p.vy)*s;
  // Update particle velocity
  p.vx = vx;
  p.vy = vy;
  p.vz = vz;  
}


/*! Project particle torque onto the normal vector. The assumption here is that 
 *  the particle's director and velocity are all already in the tangent plane 
 *  and that it is constrained to the hourglass.
 *  \param p Particle whose torque to project
*/ 
double ConstraintHourglass::project_torque(Particle& p)
{
  double Nx, Ny, Nz;
  this->compute_normal(p,Nx,Ny,Nz);
  return (p.tau_x*Nx + p.tau_y*Ny + p.tau_z*Nz);  
}
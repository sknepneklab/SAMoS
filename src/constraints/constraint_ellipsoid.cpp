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
 * \file constraint_ellipsoid.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2014
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_ellipsoid.hpp"

static void compute_gradient(double x, double y, double z, double a, double b, double c, double& gx, double& gy, double& gz)
{
  gx = 2.0*x/a/a;
  gy = 2.0*y/b/b;
  gz = 2.0*z/c/c;
}

/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = \frac{2x}{a^2}, \f$
 *  \f$ g_y = \frac{2y}{b^2}, \f$
 *  \f$ g_z = \frac{2z}{c^2}. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintEllipsoid::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  double x = p.x, y = p.y, z = p.z;
  Nx = 2.0*x/m_a/m_a;
  Ny = 2.0*y/m_b/m_b;
  Nz = 2.0*z/m_c/m_c;
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}


/*! Force particle to be confined to the surface of a ellipsoid and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface of the ellipsoid.
 *  
 *  Ellipsoid is defined by an implicit equation
 *  \f$ g(x,y,z) = \frac{x^2}{a^2} + \frac{y^2}{b^2} + \frac{z^2}{c^2} = 1\f$, where \f$ a \f$, \f$ b \f$, \f$ c \f$ are 
 *  ellipsoid parameters.
 * 
 *  \param p Particle which is to be projected onto the ellipsoid
 */
void ConstraintEllipsoid::enforce(Particle& p)
{
  // Project particle back onto the surface
  double x = p.x, y = p.y, z = p.z;
  // Compute reference gradient (SHAKE method)
  double ref_gx, ref_gy, ref_gz;
  compute_gradient(x,y,z,m_a,m_b,m_c,ref_gx, ref_gy, ref_gz);
  
  bool satisfied;
  int iter = 0;
  while (iter++ < m_max_iter)
  {
    satisfied = true;
    x = p.x; y = p.y; z = p.z;
    // Compute unit gradient of the implicit function
    double gx, gy, gz;
    compute_gradient(x,y,z,m_a,m_b,m_c,gx,gy,gz);
    double s = gx*ref_gx + gy*ref_gy + gz*ref_gz;
    double g = x*x/m_a/m_a + y*y/m_b/m_b + z*z/m_c/m_c - 1.0;  // implicit equation for ellipsoid
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
}

/*! Rotate director of a particle around the normal vector
 *  \note This function assumes that the particle has already been
 *  projected onto ellipsoid and that its director is laying in 
 *  the tangent plane.
 *  \param p Particle whose director to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintEllipsoid::rotate_director(Particle& p, double phi)
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
 *  projected onto ellipsoid and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintEllipsoid::rotate_velocity(Particle& p, double phi)
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
 *  and that it is constrained to the ellipsoid.
 *  \param p Particle whose torque to project
*/ 
double ConstraintEllipsoid::project_torque(Particle& p)
{
  double Nx, Ny, Nz;
  this->compute_normal(p,Nx,Ny,Nz);
  return (p.tau_x*Nx + p.tau_y*Ny + p.tau_z*Nz);  
}
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
 * \file constraint_gyroid.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 15-Dec-2014
 * \brief Implementation of the toroidal constraint
 */ 

#include "constraint_gyroid.hpp"

static void compute_gradient(double x, double y, double z, double lx, double ly, double lz, double& gx, double& gy, double& gz)
{
  gx = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*y)/ly)))/lx;
  gy = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*y)/ly) - sin((2*M_PI*y)/ly)*sin((2*M_PI*z)/lz)))/ly;
  gz = (2*M_PI*(cos((2*M_PI*y)/ly)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*z)/lz)))/lz;
}

/*! Compute normal to the surface at p
 *  
 *  Normal vector is computed as the normalized gradient vector at point p.
 *  Coordinates of gradient vector are given as:
 * 
 *  \f$ g_x = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  x}{\text{lx}}\right) \cos \left(\frac{2 \pi  z}{\text{lz}}\right)-\sin \left(\frac{2 \pi  x}{\text{lx}}\right) \sin \left(\frac{2 \pi 
   y}{\text{ly}}\right)\right)}{\text{lx}}, \f$
 *  \f$ g_y = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  x}{\text{lx}}\right) \cos \left(\frac{2 \pi  y}{\text{ly}}\right)-\sin \left(\frac{2 \pi  y}{\text{ly}}\right) \sin \left(\frac{2 \pi 
   z}{\text{lz}}\right)\right)}{\text{ly}}, \f$
 *  \f$ g_z = \frac{2 \pi  \left(\cos \left(\frac{2 \pi  y}{\text{ly}}\right) \cos \left(\frac{2 \pi  z}{\text{lz}}\right)-\sin \left(\frac{2 \pi  x}{\text{lx}}\right) \sin \left(\frac{2 \pi 
   z}{\text{lz}}\right)\right)}{\text{lz}}. \f$  
 *  
 *  \param p Reference to the particle object
 *  \param Nx x coordinate of the normal
 *  \param Ny y coordinate of the normal
 *  \param Nz z coordinate of the normal
*/
void ConstraintGyroid::compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
{
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lx = box->Lx, ly = box->Ly, lz = box->Lz;
  Nx = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*y)/ly)))/lx;
  Ny = (2*M_PI*(cos((2*M_PI*x)/lx)*cos((2*M_PI*y)/ly) - sin((2*M_PI*y)/ly)*sin((2*M_PI*z)/lz)))/ly;
  Nz = (2*M_PI*(cos((2*M_PI*y)/ly)*cos((2*M_PI*z)/lz) - sin((2*M_PI*x)/lx)*sin((2*M_PI*z)/lz)))/lz;
  double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
  Nx /= len_N; Ny /= len_N; Nz /= len_N;
}


/*! Force particle to be confined to the surface of a gyroid and
 *  its velocity to be tangent to it. We assume that during the integration 
 *  step particles do not move to much away from the surface and we project 
 *  them down to the surface of the gyroid.
 *  
 *  Gyroid is defined by an implicit equation
 *  \f$ g(x,y,z) = \cos \left(\frac{2 \pi  x}{\text{lx}}\right) \sin \left(\frac{2 \pi y}{\text{ly}}\right)+\sin \left(\frac{2 \pi  x}{\text{lx}}\right) \cos
   \left(\frac{2 \pi  z}{\text{lz}}\right)+\cos \left(\frac{2 \pi y}{\text{ly}}\right) \sin \left(\frac{2 \pi  z}{\text{lz}}\right)  = 0 \f$.
 * 
 *  \param p Particle which is to be projected onto the gyroid
 */
void ConstraintGyroid::enforce(Particle& p)
{
  // Project particle back onto the surface
  double x = p.x, y = p.y, z = p.z;
  BoxPtr box = m_system->get_box();
  double lx = box->Lx, ly = box->Ly, lz = box->Lz;
  bool periodic = m_system->get_periodic();
  // Compute reference gradient (SHAKE method)
  double ref_gx, ref_gy, ref_gz;
  compute_gradient(x,y,z,lx,ly,lz,ref_gx, ref_gy, ref_gz);
  
  bool satisfied;
  int iter = 0;
  while (iter++ < m_max_iter)
  {
    satisfied = true;
    x = p.x; y = p.y; z = p.z;
    // Compute unit gradient of the implicit function
    double gx, gy, gz;
    compute_gradient(x,y,z,lx,ly,lz,gx,gy,gz);
    double s = gx*ref_gx + gy*ref_gy + gz*ref_gz;
    double g = cos((2*M_PI*z)/lz)*sin((2*M_PI*x)/lx) + cos((2*M_PI*x)/lx)*sin((2*M_PI*y)/ly) + cos((2*M_PI*y)/ly)*sin((2*M_PI*z)/lz);  // implicit equation for gyroid
    double lambda = g/s;
    if (fabs(g) > m_tol)
      satisfied = false;
    if (satisfied) break;
    p.x -= lambda*ref_gx;
    p.y -= lambda*ref_gy;
    p.z -= lambda*ref_gz;
    // Check periodic boundary conditions 
    if (periodic)
    {
      if (p.x <= box->xlo) p.x += lx;
      else if (p.x >= box->xhi) p.x -= lx;
      if (p.y <= box->ylo) p.y += ly;
      else if (p.y >= box->yhi) p.y -= ly;
      if (p.z <= box->zlo) p.z += lz;
      else if (p.z >= box->zhi) p.z -= lz;
    }
    else // reflective boundary conditions
    {
      if (p.x < box->xlo) 
      {
        p.x = box->xlo;
        p.vx = -p.vx;
      }
      else if (p.x > box->xhi)
      {
        p.x = box->xhi;
        p.vx = -p.vx;
      }
      if (p.y < box->ylo) 
      {
        p.y = box->ylo;
        p.vy = -p.vy;
      }
      else if (p.y > box->yhi)
      {
        p.y = box->yhi;
        p.vy = -p.vy;
      }
      if (p.z < box->zlo) 
      {
        p.z = box->zlo;
        p.vz = -p.vz;
      }
      else if (p.z > box->zhi)
      {
        p.z = box->zhi;
        p.vz = -p.vz;
      }
    }
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
 *  projected onto gyroid and that its director is laying in 
 *  the tangent plane.
 *  \param p Particle whose director to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintGyroid::rotate_director(Particle& p, double phi)
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
 *  projected onto gyroid and that its velocity is laying in 
 *  the tangent plane.
 *  \param p Particle whose velocity to rotate
 *  \param phi angle by which to rotate it
*/
void ConstraintGyroid::rotate_velocity(Particle& p, double phi)
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
 *  and that it is constrained to the gyroid.
 *  \param p Particle whose torque to project
*/ 
double ConstraintGyroid::project_torque(Particle& p)
{
  double Nx, Ny, Nz;
  this->compute_normal(p,Nx,Ny,Nz);
  return (p.tau_x*Nx + p.tau_y*Ny + p.tau_z*Nz);  
}
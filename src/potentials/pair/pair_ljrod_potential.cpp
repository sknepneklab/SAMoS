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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_ljrod_potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Mar-2015
 * \brief Implementation of PairLJRodPotential class
 */ 

#include "pair_ljrod_potential.hpp"

//! \param dt time step sent by the integrator 
//! Algorithm for determining closest between two rods follows closely 
//! "A fast algorithm to evaluate the shortest distance between rods",
//!  C. Vega, S. Lago, Computers & Chemistry, Volume 18, 55–59 (1994)
void PairLJRodPotential::compute(double dt)
{
  int N = m_system->size();
  double sigma = m_sigma;
  double eps = m_eps;
  double rcut = m_rcut;
  double sigma_sq = sigma*sigma, rcut_sq = rcut*rcut;
  double li2, lj2;
  double lambda, mu;
  double dx, dy, dz;
  double fx, fy, fz;
  double force_factor;
  double potential_energy;
  double alpha = 1.0;  // phase in factor
  double k;
  double inv_core_sq, inv_core_6, lj_core_sq;
  bool parallel;
  double A, B, C, D, overlap;
  double ifact, jfact;
    
  if (m_system->compute_per_particle_energy())
  {
    for  (int i = 0; i < N; i++)
    {
      Particle& p = m_system->get_particle(i);
      p.set_pot_energy("ljrod",0.0);
    }
  }
  
  m_potential_energy = 0.0;
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_phase_in)
      alpha = m_val->get_val(static_cast<int>(pi.age/dt));
    li2 = 0.5*pi.get_length();
    double ni_x = pi.nx, ni_y = pi.ny, ni_z = pi.nz; 
    vector<int>& neigh = m_nlist->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      if (m_has_pair_params)
      {
        rcut = m_pair_params[pi.get_type()-1][pj.get_type()-1].rcut;
        rcut_sq = rcut*rcut;
      }
      lj2 = 0.5*pj.get_length();
      double nj_x = pj.nx, nj_y = pj.ny, nj_z = pj.nz;
      double dx_cm = pj.x - pi.x, dy_cm = pj.y - pi.y, dz_cm = pj.z - pi.z;
      m_system->apply_periodic(dx_cm,dy_cm,dz_cm);
      
      double ni_dot_nj = ni_x*nj_x + ni_y*nj_y + ni_z*nj_z; 
      double drcm_dot_ni = dx_cm*ni_x + dy_cm*ni_y + dz_cm*ni_z;
      double drcm_dot_nj = dx_cm*nj_x + dy_cm*nj_y + dz_cm*nj_z;
      double cc = 1.0 - ni_dot_nj*ni_dot_nj;
      
      parallel = false;
      
      if (cc <= 1e-6) // rods are nearly parallel
      {
        parallel = true;
      }
      else
      {
        lambda = (drcm_dot_ni-ni_dot_nj*drcm_dot_nj)/cc;
        mu = (-drcm_dot_nj+ni_dot_nj*drcm_dot_ni)/cc;
        if (fabs(lambda) > li2 || fabs(mu) > lj2)
        {
          if (fabs(lambda) - li2 > fabs(mu) - lj2)
          {
            lambda = copysign(li2,lambda);
            mu = lambda*ni_dot_nj - drcm_dot_nj;
            if (fabs(mu) > lj2) mu = copysign(lj2,mu);
          }
          else
          {
            mu =  copysign(lj2,mu);
            lambda = mu*ni_dot_nj + drcm_dot_ni;
            if (fabs(lambda) > li2) lambda = copysign(li2,lambda);      
          }
        }
      }
      
      if (parallel) // if rods are parallel, distance vector is \delta\vec r = \delta\vec r_cm - (\delta \vec r_cm \dot \vec n_i)\vec n_i
      {
        dx = dx_cm - drcm_dot_ni*ni_x;
        dy = dy_cm - drcm_dot_ni*ni_y;
        dz = dz_cm - drcm_dot_ni*ni_z;
        A = fmax(-li2,-lj2+drcm_dot_ni);
        B = fmin(li2,lj2+drcm_dot_ni);
        C = fmax(-lj2,-li2-drcm_dot_nj);
        D = fmin(lj2,li2-drcm_dot_nj);
        overlap = fmax(0.0,B-A);
      }
      else
      {
        double x1 = pi.x + lambda*ni_x, y1 = pi.y + lambda*ni_y, z1 = pi.z + lambda*ni_z;
        double x2 = pj.x + mu*nj_x, y2 = pj.y + mu*nj_y, z2 = pj.z + mu*nj_z;
        m_system->apply_periodic(x1,y1,z1);
        m_system->apply_periodic(x2,y2,z2);
        dx = x2 - x1; dy = y2 - y1; dz = z2 - z1;
        //dx = dx_cm + mu*nj_x - lambda*ni_x;
        //dy = dy_cm + mu*nj_y - lambda*ni_y;
        //dz = dz_cm + mu*nj_z - lambda*ni_z;
      }
      m_system->apply_periodic(dx,dy,dz);
      
      double r_sq = dx*dx + dy*dy + dz*dz;
      
      if (r_sq <= rcut_sq)
      { 
        if (m_has_pair_params)
        {
          int pi_t = pi.get_type() - 1, pj_t = pj.get_type() - 1;
          sigma = m_pair_params[pi_t][pj_t].sigma;
          eps = m_pair_params[pi_t][pj_t].eps;
          sigma_sq = sigma*sigma;
        }
        if (r_sq <= LJ_HARD_CORE_DISTANCE*LJ_HARD_CORE_DISTANCE*sigma_sq)
        {
          lj_core_sq = LJ_HARD_CORE_DISTANCE*LJ_HARD_CORE_DISTANCE;
          inv_core_sq = sigma_sq/lj_core_sq;
          inv_core_6 = inv_core_sq*inv_core_sq*inv_core_sq;
          k = 6.0*eps/(lj_core_sq*lj_core_sq)*inv_core_6*(2.0*inv_core_6 - 1.0);
          potential_energy = 4.0*eps*alpha*inv_core_6*(inv_core_6 - 1.0);
          if (r_sq <= SMALL_NUMBER)
          {
            dx = dx_cm; dy = dy_cm; dz = dz_cm;
            r_sq = dx_cm*dx_cm + dy_cm*dy_cm + dz_cm*dz_cm;
          }
          force_factor = 4.0*k*r_sq*sqrt(r_sq);
          //cout << "Too close. F = " << force_factor << endl;
        }
        else
        {
          double inv_r_sq = sigma_sq/r_sq;
          double inv_r_6  = inv_r_sq*inv_r_sq*inv_r_sq;
          // Handle potential 
          potential_energy = 4.0*eps*alpha*inv_r_6*(inv_r_6 - 1.0);
          if (m_shifted)
          {
            double inv_r_cut_sq = sigma_sq/rcut_sq;
            double inv_r_cut_6 = inv_r_cut_sq*inv_r_cut_sq*inv_r_cut_sq;
            potential_energy -= 4.0 * eps * alpha * inv_r_cut_6 * (inv_r_cut_6 - 1.0);
          }
          force_factor = 48.0*eps*alpha*inv_r_6*(inv_r_6 - 0.5)*inv_r_sq;
          //cout << "LJ part with F = " << force_factor << endl;
        }
        m_potential_energy += potential_energy;
        // Handle force
        fx = force_factor*dx;  fy = force_factor*dy;  fz = force_factor*dz;
        if (parallel)
        {
          fx *= overlap;  fy *= overlap;  fz *= overlap;
        }
        pi.fx -= fx;
        pi.fy -= fy;
        pi.fz -= fz;
        // Use 3d Newton's law
        pj.fx += fx;
        pj.fy += fy;
        pj.fz += fz;
        // handle torques
        if (parallel)
        {
          ifact = 0.5*(B*B - A*A);
          jfact = 0.5*(D*D - C*C);
        }
        else
        {
          ifact = lambda;
          jfact = mu;
        }
        pi.tau_x -= ifact*(ni_y*fz - ni_z*fy);
        pi.tau_y -= ifact*(ni_z*fx - ni_x*fz);
        pi.tau_z -= ifact*(ni_x*fy - ni_y*fx);
        pj.tau_x += jfact*(nj_y*fz - nj_z*fy);
        pj.tau_y += jfact*(nj_z*fx - nj_x*fz);
        pj.tau_z += jfact*(nj_x*fy - nj_y*fx);
        if (m_system->compute_per_particle_energy())
        {
          pi.add_pot_energy("ljrod",potential_energy);
          pj.add_pot_energy("ljrod",potential_energy);
        }
      }
    }
  }
}
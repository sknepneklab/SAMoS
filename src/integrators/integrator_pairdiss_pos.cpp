/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

/*!
 * \file integrator_brownian_pos.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \author Silke Henkes shenkes@lorentz.leidenuniv.nl
 * \date 8 April 2024
 * \brief Implementation of the Pair dissipation dynamics integrator for particle position.
 */ 

#include "integrator_pairdiss_pos.hpp"

/*! Integrates equation of motion in the over-damped limit using a first order 
 *  scheme. Computes the contact matrix first, and then uses Eigen Sparse matrix to invert it
 *  Equation of motion:
 *  \zeta \dot{x} = F + \xipair C \dot{x}, where C is the mattrix of contacts, i.e.
 *  C_iajb = \delta_ab (-1 if i!=j and z_i if i = j)
 *  (\zeta I + \xipair C) \dot{x} = F
 * \dot{x} = (\zeta I + \xipair C)^{-1} F = M^{-1} F
 *  with M_iajb = \delta_ab (-xipair if i!=j and zeta + z_i \xipair if i = j)
 * This is a Hermitian D. P. matrix, and a close cousin to a Dynamical matrix in granular. Always invertible as long as zeta>0, xipair>0.
 *  \note This integrator applies only to the particle position and does not implement activity.
 *  In order to use activity, you should define, e.g. external self propulsion.  
**/

// Following the same steps as Yann Keta's vertex model implementation
// while I cam factor out the \delta_ab part completely, leading to a matrix that's only NxN, not 3N x 3N
// We never actually construct the full inverse matrix. So need to keep the 3Nx3N approach here

// 25.11.24: This has a bug: if we run this for any subgroup of a system whih is not the whole system, there are interactions between group and non-group
// particles which can't formulated like this. Hence they crash. Need to always compute with everything, though we can throw away the updates for non-group members
void IntegratorPairdissPos::integrate()
{

    // number of particles sets size
    //int N = m_system->get_group(m_group_name)->get_size();
    int N = m_system->size();

    // Eigen linear algebra components
    // triplet i,j,value of nonzero sparse entries
    vector<Eigen::Triplet<double>> modContactTrip(0);
    // actual sparse matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> modContact(3*N,3*N);
    // velocities and forces as full vectors
    Eigen::VectorXd generalForces(3*N);
    Eigen::VectorXd generalVelocities(3*N);
    // velocities will be class variable

    // Fill the matrix
    double alpha_i = 1.0;  // phase in factor for particle i
    double alpha_j = 1.0;  // phase in factor for particle j
    double alpha = 1.0; // phase in factor for pair interaction (see below)
    //vector<int> particles = m_system->get_group('all')->get_particles();
    vector<double> zpart;
    for (int i = 0; i < N; i++) zpart.push_back(0.0);
    for (int i = 0; i < N; i++)
    {
        //int plab = particles[i];

        //Particle& pi = m_system->get_particle(plab);
        Particle& pi = m_system->get_particle(i);


        if (m_phase_in) {
            alpha_i = 0.5*(1.0 + std::min(1.0,pi.age/m_tphase));
        }

        // locate neighbours
        vector<int>& neigh = m_nlist->get_neighbours(i);
        // now determine which ones are actually interacting ...
        for (unsigned int j = 0; j < neigh.size(); j++)
        {
            Particle& pj = m_system->get_particle(neigh[j]);
            if (m_phase_in)
            {
                alpha_j = 0.5*(1.0 + std::min(1.0,pi.age/m_tphase));
                // Determine global phase in factor: particles start at 0.5 strength (both daugthers of a division replace the mother)
                // Except for the interaction between daugthers which starts at 0
                if (alpha_i < 1.0 && alpha_j < 1.0)
	                alpha = alpha_i + alpha_j - 1.0;
	            else 
	                alpha = alpha_i*alpha_j;
            }
      
            double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
            m_system->apply_periodic(dx,dy,dz);
            double r_sq = dx*dx + dy*dy + dz*dz;
            double r = sqrt(r_sq);

            double ai_p_aj;
            if (!m_use_particle_radii)
                ai_p_aj = 2.0;
            else
                ai_p_aj = pi.get_radius() + pj.get_radius();
            // interaction radius
            double rintij = m_r_int*ai_p_aj;

            if (r < rintij) // we are in contact
            {
                for (int dim=0; dim<3; dim++) 
                {
                    modContactTrip.push_back(Eigen::Triplet<double>(3*i+dim,3*neigh[j]+dim,-alpha*m_xipair));
                    // and the Neighbour list is one-sided ...
                    modContactTrip.push_back(Eigen::Triplet<double>(3*neigh[j]+dim,3*i+dim,-alpha*m_xipair));
                }
                // Sum up friction contribution in both directions
                zpart[i] += alpha;
                zpart[neigh[j]] += alpha;
            }
        }
    }
    // And now that's all done we can do the diagonal properly
    double zav=0;
    for (int i = 0; i < N; i++)
    {
        for (int dim=0; dim<3; dim++) 
        {
            modContactTrip.push_back(Eigen::Triplet<double>(3*i+dim,3*i+dim,m_zeta + m_xipair*zpart[i]));
        }
        zav+=zpart[i];
    }
    zav = zav/N;

    // Which now will allow us to solve the linear system for the velocities
    // set sparse matrix to zero
    modContact.setZero();
    // and fill with the triplets we just constructed
    modContact.setFromTriplets(modContactTrip.begin(),modContactTrip.end());

    // Set the forces vector
    for (int i = 0; i < N; i++)
    {
        //int pi = particles[i];
        //Particle& p = m_system->get_particle(pi);
        Particle& p = m_system->get_particle(i);
        generalForces[3*i]=p.fx;
        generalForces[3*i+1]=p.fy;
        generalForces[3*i+2]=p.fz;
    }

    // Eigen does the actual solution
    // compute decomposition
    sol.compute(modContact);
    // check that it's working
    assert(sol.info()==Eigen::Success);
    // invert matrix by linear algebra to compute the velocities
    generalVelocities = sol.solve(generalForces);
    // and check again that it's working
    assert(sol.info()==Eigen::Success);

    // now we can do the actual update of the velocities and positions


    // FDT following Doi + Edwards page 50-55 and a bit of work
    // The noise needs to follow the correlation function
    // <fn fm> = 2 L_nm k_b T delta(t-t'), where L_nm is the mobility matrix, i.e. M^{-1}
    // To generate that kind of noise is nontrivial. Method:
    // Find the eigenvalues L_mu of the mobility matrix. Compute the eigenvalues D_mu = sqrt{2 L_mu k_b T} 
    // Generate the matrix D = \Sigma^{-1} diag(sqrt{2 L_mu k_b T}) \Sigma, where Sigma is the eigenvector matrix of L
    // then the appropirate noise fm = D \eta, where \eta is Gaussian white
    // Uncertain how to do this in Eigen
    // Stopgap measure: approximate mobility by mu = 1/(zeta + z xipair/2)
    // TO BE IMPLEMENTD PROPERLY
    
    // general integrator things
    double T = m_temp->get_val(m_system->get_run_step());
    double sqrt_dt = sqrt(m_dt);
    double fr_x = 0.0, fr_y = 0.0, fr_z = 0.0;  // Random part of the force
    double muMF = 1.0/(m_zeta+0.5*zav*m_xipair);
    double B = sqrt(2.0*muMF*T);
    if (T > 0.0)
    {
        std::cout<<"Warning: Non-FDT compliant thermal fluctuations. Do not use as thermal integrator in current form."<<std::endl;
    }

    // iterate over all particles of only the correct group
    // Finally check for group and update only if it's the correct one
    vector<int> particles = m_system->get_group(m_group_name)->get_particles();
    int Ngroup = m_system->get_group(m_group_name)->get_size(); 
    for (int i = 0; i < Ngroup; i++)
    {
        int pi = particles[i];

        Particle& p = m_system->get_particle(pi);
        // Update velocity
        p.vx = generalVelocities[3*i];
        p.vy = generalVelocities[3*i+1];
        p.vz = generalVelocities[3*i+2];
        // Update particle position 
        p.x += m_dt*p.vx;
        p.y += m_dt*p.vy;
        p.z += m_dt*p.vz;
        // Check is non-zero T and if non-zero add stochastic part
        // TO BE IMPLEMENTED PROPERLY
        if (T > 0.0)
        {
            fr_x = B*m_rng->gauss_rng(1.0);
            fr_y = B*m_rng->gauss_rng(1.0);
            fr_z = B*m_rng->gauss_rng(1.0);
            p.vx += fr_x; 
            p.vy += fr_y;
            p.vz += fr_z;  
            p.x += sqrt_dt*fr_x;
            p.y += sqrt_dt*fr_y;
            p.z += sqrt_dt*fr_z;
        }
        // Project everything back to the manifold
        m_constrainer->enforce(p);
        p.age += m_dt;
    }
    // Update vertex mesh
    m_system->update_mesh();

}


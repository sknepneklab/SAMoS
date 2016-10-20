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
 * \file population_actomyosin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 08-Jun-2016
 * \brief Implementation of PopulationActomyosin class.
 */ 

#include "population_actomyosin_molecule.hpp"


/*! This function controls detachment of the entire myosin molecule. This is a
 *  achieved via changing types of myosin head groups from attached (type "A") to detached (type "D") a
 *  and vice versa. 
 *
 *  \note In order to avoid head groups, e.g. changing  their type from "A" to "D" and back to "A" in the 
 *  single time step, we implement both processes in the same function using lists of indices.
 *  \param t current time step
 *  
*/
void PopulationActomyosinMolecule::divide(int t)
{
  if ((m_freq > 0) && (t % m_freq == 0) && (m_detach_rate > 0.0))  // Attempt D to A transition only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before actomyosin molecule detachment: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    int Nmol = m_system->number_of_molecules();
    // Probability of detachment for a given particle is rate per particle multiplied with time,
    // where time is equal to m_freq*integrator_time_step.   
    double detach_prob = m_detach_rate*m_freq*m_system->get_integrator_step();
    // loop over all molecules in the system
    for (int mol = 0; mol < Nmol; mol++)
    {
      vector<int>& particles = m_system->get_mol_particles(mol);  // get all particle ids that belong to this molecule
      bool detach = false;
      for (unsigned int i = 0; i < particles.size(); i++)   // loop over all particles in the molecule
      {
        Particle& pi = m_system->get_particle(particles[i]);
        if (pi.get_type() == m_type_a)    // if at least one of the beads is of type "attached", mark for detachment
        {
          detach = true;
          break;
        }
      }
      // if marked for detachment and random number below detachment probability, make transition attached --> detached on all beads
      if (detach && m_rng->drnd() < detach_prob)
        for (unsigned int i = 0; i < particles.size(); i++)   // loop over all particles in the molecule
        {
          Particle& pi = m_system->get_particle(particles[i]);
          if (pi.get_type() == m_type_a)    
            pi.set_type(m_type_d);          
        }
    }
  }
}

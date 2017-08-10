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
 * \file potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Oct-2013
 * \brief Definition of Potential class members.
*/ 

#include "potential.hpp"

/*! Iterate over all pair and external potential and compute 
 *  potential energies and forces
 *  \param dt step size (used to phase in particles)
 */
void Potential::compute(double dt)
{
  //m_system->reset_forces();
  PairPotType::iterator it_pair;
  ExternPotType::iterator it_ext;
  BondPotType::iterator it_bond;
  AnglePotType::iterator it_angle;
  
  for(it_pair = m_pair_interactions.begin(); it_pair != m_pair_interactions.end(); it_pair++)
    (*it_pair).second->compute(dt);
  for(it_ext = m_external_potentials.begin(); it_ext != m_external_potentials.end(); it_ext++)
    (*it_ext).second->compute();
  for(it_bond = m_bond.begin(); it_bond != m_bond.end(); it_bond++)
    (*it_bond).second->compute();
  for(it_angle = m_angle.begin(); it_angle != m_angle.end(); it_angle++)
    (*it_angle).second->compute();
}

/*! Iterate over all pair, external, bond and angle force computes  
 *  and calculate total potential energy of the system.
 */
double Potential::compute_potential_energy()
{
  double pot_eng = 0.0;
  PairPotType::iterator it_pair;
  ExternPotType::iterator it_ext;
  BondPotType::iterator it_bond;
  AnglePotType::iterator it_angle;
  
  for(it_pair = m_pair_interactions.begin(); it_pair != m_pair_interactions.end(); it_pair++)
    pot_eng += (*it_pair).second->get_potential_energy();
  for(it_ext = m_external_potentials.begin(); it_ext != m_external_potentials.end(); it_ext++)
    pot_eng += (*it_ext).second->get_potential_energy();
  for(it_bond = m_bond.begin(); it_bond != m_bond.end(); it_bond++)
    pot_eng += (*it_bond).second->get_potential_energy();
  for(it_angle = m_angle.begin(); it_angle != m_angle.end(); it_angle++)
    pot_eng += (*it_angle).second->get_potential_energy();

  return pot_eng;
}

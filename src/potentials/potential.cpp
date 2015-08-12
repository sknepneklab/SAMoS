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
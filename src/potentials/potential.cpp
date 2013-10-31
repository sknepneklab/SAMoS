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
 * \file potential.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Oct-2013
 * \brief Definition of Potential class members.
*/ 


/*! Iterate over all pair and external potential and compute 
 *  potential energies and forces
 */
void Potential::compute()
{
  m_system->reset_forces();
  PairPotType::iterator it_pair;
  ExternPotType::iterator it_ext;
  
  for(it_pair = m_pair_interactions.begin(); it_pair != m_pair_interactions.end(); it_pair++)
    (*it_pair).second->compute();
  for(it_ext = m_external_potentials.begin(); it_ext != m_external_potentials.end(); it_ext++)
    (*it_ext).second->compute();
}
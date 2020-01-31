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
 * \file potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Oct-2013
 * \brief Declaration of Potential class.
 */ 

#ifndef __POTENTIAL_HPP__
#define __POTENTIAL_HPP__

#include <map>
#include <string>
#include <vector>

using std::map;
using std::string;

#include "messenger.hpp"
#include "system.hpp"

#include "pair_potential.hpp"
#include "external_potential.hpp"
#include "bond_potential.hpp"
#include "angle_potential.hpp"

typedef map<string,PairPotentialPtr> PairPotType;
typedef map<string,ExternalPotentialPtr> ExternPotType;
typedef map<string,BondPotentialPtr> BondPotType;
typedef map<string,AnglePotentialPtr> AnglePotType;

/*! Potential class handles all potentials (pair interactions, external forces, bonds and angles) present 
 *  in the system. All potentials are stored in two STL maps which both have strings as keys
 *  while items are pointers (std shared_ptr) to the PairPotential and ExternalPotential abstract classes,
 *  respectively. The integrator part of the code calls member functions of this class
 *  to perform actual force and potential computations.
*/
class Potential
{
public:
  
  //! Construct Potential object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  Potential(SystemPtr sys, const MessengerPtr msg) : m_system(sys), m_msg(msg), m_need_nlist(false) { }
  
  //! Destructor
  ~Potential()
  {
    m_pair_interactions.clear();
    m_external_potentials.clear();
    m_bond.clear();
    m_angle.clear();
  }
  
  //! Add pair potential to the list of all pair interactions
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the pair potential object
  void add_pair_potential(const string& name, PairPotentialPtr pot)
  {
    m_pair_interactions[name] = pot;
    m_need_nlist = m_need_nlist | pot->need_nlist();
    m_msg->msg(Messenger::INFO,"Added pair potential : " + name + " to the list of pair interactions.");
    if (pot->need_nlist())
      m_msg->msg(Messenger::INFO,"Pair potential " + name + " has neighbour list. Neighbour list updates will be performed during the simulation.");
  }
  
  //! Add external potential to the list of all external potentials
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the external potential object
  void add_external_potential(const string& name, ExternalPotentialPtr pot)
  {
    m_external_potentials[name] = pot;
    m_msg->msg(Messenger::INFO,"Added external potential : " + name + " to the list of external forces.");
  }
  
  //! Add bond potential to the list of all bond interactions
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the bond potential object
  void add_bond_potential(const string& name, BondPotentialPtr pot)
  {
    m_bond[name] = pot;
    m_msg->msg(Messenger::INFO,"Added bond potential : " + name + " to the list of bond interactions.");
  }
  
  //! Add angle potential to the list of all angle interactions
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the angle potential object
  void add_angle_potential(const string& name, AnglePotentialPtr pot)
  {
    m_angle[name] = pot;
    m_msg->msg(Messenger::INFO,"Added angle potential : " + name + " to the list of angle interactions.");
  }
  
  //! Add pair potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_pair_potential_parameters(const string& name, pairs_type& params)
  {
    m_pair_interactions[name]->set_pair_parameters(params);
  }
  
  //! Add pair type parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new type parameters
  void add_pair_type_parameters(const string& name, pairs_type& params)
  {
    m_pair_interactions[name]->set_type_parameters(params);
  }
  
  //! Add external potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_external_potential_parameters(const string& name, pairs_type& params)
  {
    m_external_potentials[name]->set_parameters(params);
  }
  
  //! Add bond potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_bond_potential_parameters(const string& name, pairs_type& params)
  {
    m_bond[name]->set_bond_parameters(params);
  }
  
  //! Add angle potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_angle_potential_parameters(const string& name, pairs_type& params)
  {
    m_angle[name]->set_angle_parameters(params);
  }
  
  //! Compute total pair potential energy of a given type (for measurement)
  //! \param type pair potential type
  double compute_pair_potential_energy_of_type(const string& type)
  {
    if (m_pair_interactions.find(type) == m_pair_interactions.end())
    {
      m_msg->msg(Messenger::ERROR,"Trying to compute pair potential of type " + type + " that is not defined for this system.");
      throw runtime_error("Pair potential of type " + type + " not defined.");
    }
    return m_pair_interactions[type]->get_potential_energy();
  }
  
  //! Compute total external potential energy of a given type (for measurement)
  //! \param type pair potential type
  double compute_external_potential_energy_of_type(const string& type)
  {
    if (m_external_potentials.find(type) == m_external_potentials.end())
    {
      m_msg->msg(Messenger::ERROR,"Trying to compute external potential of type " + type + " that is not defined for this system.");
      throw runtime_error("External potential of type " + type + " not defined.");
    }
    return m_external_potentials[type]->get_potential_energy();
  }
  
  //! Compute total bond potential energy of a given type (for measurement)
  //! \param type bond potential type
  double compute_bond_potential_energy_of_type(const string& type)
  {
    if (m_bond.find(type) == m_bond.end())
    {
      m_msg->msg(Messenger::ERROR,"Trying to compute bond potential of type " + type + " that is not defined for this system.");
      throw runtime_error("Bond potential of type " + type + " not defined.");
    }
    return m_bond[type]->get_potential_energy();
  }
  
  //! Compute total angle potential energy of a given type (for measurement)
  //! \param type angle potential type
  double compute_angle_potential_energy_of_type(const string& type)
  {
    if (m_angle.find(type) == m_angle.end())
    {
      m_msg->msg(Messenger::ERROR,"Trying to compute angle potential of type " + type + " that is not defined for this system.");
      throw runtime_error("Angle potential of type " + type + " not defined.");
    }
    return m_angle[type]->get_potential_energy();
  }
  
  //! Returns true if any of the potentials need neighbour list
  bool need_nlist() { return m_need_nlist; } 
  
  // Compute total potential energy
  double compute_potential_energy();

  //! Compute all forces and potentials in the system
  void compute(double);
  
private:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  
  PairPotType m_pair_interactions;      //!< Contains information about all pair interactions
  ExternPotType m_external_potentials;  //!< Contains information about all external forces
  BondPotType m_bond;                   //!< Contains information about all bonds
  AnglePotType m_angle;                 //!< Contains information about all angles
  
  bool m_need_nlist;                  //!< If true, there are potentials that need neighbour list
   
};

typedef shared_ptr<Potential> PotentialPtr;

#endif

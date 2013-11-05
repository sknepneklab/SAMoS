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
 * \file potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Oct-2013
 * \brief Declaration of Potential class.
 */ 

#ifndef __POTENTIAL_HPP__
#define __POTENTIAL_HPP__

#include <map>
#include <string>

using std::map;
using std::string;

#include "messenger.hpp"
#include "system.hpp"

#include "pair_potential.hpp"
#include "external_potential.hpp"

typedef map<string,PairPotentialPtr> PairPotType;
typedef map<string,ExternalPotentialPtr> ExternPotType;

/*! Potential class handles all potentials (pair interactions and external forces) present 
 *  in the system. All potentials are stored in two STL maps which both have strings as keys
 *  while items are pointers (boost shared_ptr) to the PairPotential and ExternalPotential abstract classes,
 *  respectively. The integrator part of the code calls member functions of this class
 *  to perform actual force and potential computations.
*/
class Potential
{
public:
  
  //! Construct Potential object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  Potential(SystemPtr sys, const MessengerPtr msg) : m_system(sys), m_msg(msg) { }
  
  //! Destructor
  ~Potential()
  {
    m_pair_interactions.clear();
    m_external_potentials.clear();
  }
  
  //! Add pair potential to the list of all pair interactions
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the pair potential object
  void add_pair_potential(const string& name, PairPotentialPtr pot)
  {
    m_pair_interactions[name] = pot;
    m_msg->msg(Messenger::INFO,"Added pair potential : " + name + " to the list of pair interactions.");
  }
  
  //! Add external potential to the list of all external potentials
  //! \param name Unique name of the potential 
  //! \param pot Pointer to the external potential object
  void add_external_potential(const string& name, ExternalPotentialPtr pot)
  {
    m_external_potentials[name] = pot;
    m_msg->msg(Messenger::INFO,"Added external potential : " + name + " to the list of external forces.");
  }
  
  //! Add pair potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_pair_potential_parameters(const string& name, pairs_type& params)
  {
    m_pair_interactions[name].set_pair_parameters(params);
  }
  
  //! Add external potential parameters
  //! \param name Unique name of the potential 
  //! \param params maps with new potential parameters
  void add_external_potential_parameters(const string& name, pairs_type& params)
  {
    m_external_potentials[name].set_parameters(params);
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
    return m_pair_interactions[type].second->get_potential_energy();
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
    return m_external_potentials[type].second->get_potential_energy();
  }
  
  //! Compute all forces and potentials in the system
  void compute();
  
private:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  
  PairPotType m_pair_interactions;    //!< Contains information about all pair interactions
  ExternPotType m_external_potentials;  //!< Contains information about all external forces
  
};

typedef shared_ptr<Potential> PotentialPtr;

#endif
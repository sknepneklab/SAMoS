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
 * \file bond_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Declaration of BondPotential class
 */ 

#ifndef __BOND_POTENTIAL_HPP__
#define __BOND_POTENTIAL_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;


/*! BondPotential is the base class (abstract) that handles
 *  all calls to the bond potential evaluations. Children 
 *  of this class will implement actual force fields.
 */
class BondPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  BondPotential(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), 
                                                                      m_msg(msg),
                                                                      m_has_bond_params(false),
                                                                      m_use_particle_radii(false)
                                                                      { }
                                                                                                       
  //! Destructor 
  virtual ~BondPotential() { }
  
                                                                                                   
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set bond parameters   
  virtual void set_bond_parameters(pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  bool m_has_bond_params;         //!< Flag that controls if bond parameters are set
  double m_potential_energy;       //!< Total potential energy
  bool m_use_particle_radii;       //!< If true, base native bond length on particle radii
  
};


typedef shared_ptr<BondPotential> BondPotentialPtr;


#endif

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
 * \file angle_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Declaration of AnglePotential class
 */ 

#ifndef __ANGLE_POTENTIAL_HPP__
#define __ANGLE_POTENTIAL_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;


/*! AnglePotential is the base class (abstract) that handles
 *  all calls to the angle potential evaluations. Children 
 *  of this class will implement actual force fields.
 */
class AnglePotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  AnglePotential(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), 
                                                                       m_msg(msg),
                                                                       m_has_angle_params(false)
                                                                       { }
                                                                                                       
  //! Destructor 
  virtual ~AnglePotential() { }
  
                                                                                                   
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set angle parameters   
  virtual void set_angle_parameters(pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  bool m_has_angle_params;         //!< Flag that controls if angle parameters are set
  double m_potential_energy;       //!< Total potential energy
  
};


typedef shared_ptr<AnglePotential> AnglePotentialPtr;


#endif

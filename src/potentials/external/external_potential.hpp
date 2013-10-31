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
 * \file external_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Oct-2013
 * \brief Declaration of ExternalPotential class
 */ 

#ifndef __EXTERNAL_POTENTIAL_HPP__
#define __EXTERNAL_POTENTIAL_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;

typedef map<int, map<string,double> > ExternData;

/*! ExternalPotential is the base class (abstract) that handles
 *  all calls to the external potential evaluations. Children 
 *  of this class will implement actual external potentials.
 */
class ExternalPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalPotential(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg),
                                                                          m_param(param),
                                                                          m_has_params(false)
                                                                          { }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  virtual void set_parameters(int, pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  pairs_type m_param;              //!< Handles global potential parameters
  ExternData m_type_params;        //!< Handles specific parameters for a given particle type
  bool m_has_params;               //!< Flag that controls if pair parameters are set
  double m_potential_energy;       //!< Total potential energy
  
};


typedef shared_ptr<ExternalPotential> ExternalPotentialPtr;


#endif

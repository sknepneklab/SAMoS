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
                                                                          m_has_params(false),
                                                                          m_exclude_boundary(false)
  {
    if (param.find("exclude_boundary") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"External potential/force will not be applied onto boundary particles.");
      m_exclude_boundary = true;
    }
   }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  virtual void set_parameters(pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  pairs_type m_param;              //!< Handles global potential parameters
  ExternData m_type_params;        //!< Handles specific parameters for a given particle type
  bool m_has_params;               //!< Flag that controls if pair parameters are set
  double m_potential_energy;       //!< Total potential energy
  bool m_exclude_boundary;         //!< If true, do not apply external potentail to praticle with boundary flag set to true
  
};


typedef shared_ptr<ExternalPotential> ExternalPotentialPtr;


#endif

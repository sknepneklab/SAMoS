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
 * \file external_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2014
 * \brief Declaration of ExternalAlign class
 */ 

#ifndef __EXTERNAL_ALIGN_HPP__
#define __EXTERNAL_ALIGN_HPP__

#include <map>
#include <vector>
#include <string>

#include "system.hpp"

#include "parse_parameters.hpp"

using std::map;
using std::pair;
using std::string;
using std::vector;

typedef map<int, map<string,double> > ExternAlignData;

/*! ExternalAlign is the base class (abstract) that handles
 *  all calls to the external alignment evaluations. Children 
 *  of this class will implement actual external alignment.
 */
class ExternalAlign
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg),
                                                                      m_param(param),
                                                                      m_has_params(false)
                                                                      { }
    
  //! Set pair parameters data for external alignment 
  virtual void set_parameters(pairs_type&) = 0;
  
  //! Computes potentials and forces for all particles
  virtual void compute() = 0;
  
  
protected:
       
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  pairs_type m_param;              //!< Handles global alignment parameters
  ExternAlignData m_type_params;   //!< Handles specific parameters for a given particle type
  bool m_has_params;               //!< Flag that controls if pair parameters are set
};


typedef shared_ptr<ExternalAlign> ExternalAlignPtr;


#endif

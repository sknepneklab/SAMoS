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
 * \file integrator.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Declaration of Integrator class
 */ 

#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

#include <boost/make_shared.hpp>

#include "system.hpp"
#include "messenger.hpp"
#include "potential.hpp"
#include "constraint.hpp"
#include "aligner.hpp" 
#include "parse_parameters.hpp"

using boost::make_shared;

/*! Integrator class is the base class for handling different numerical 
 *  integrators for integrating equations of motion (e.g., NVE).
 *  This is an abstract class and its children will implement 
 *  actual integrators.
*/
class Integrator
{
public:
  
  //! Constructor
  //! \param sys Pointer to a System object containing all particles
  //! \param msg Internal message handler
  //! \param pot Pairwise and external interaction handler
  //! \param align Pairwise and external alignment handler
  //! \param nlist Neighbour list object
  //! \param cons Enforces constraints to the manifold surface
  //! \param param Contains information about all parameters 
  Integrator(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist, ConstraintPtr cons, pairs_type& param) : m_system(sys),
                                                                                                                                                   m_msg(msg),
                                                                                                                                                   m_potential(pot),
                                                                                                                                                   m_align(align),
                                                                                                                                                   m_nlist(nlist),
                                                                                                                                                   m_constraint(cons)
  { 
    if (param.find("dt") == param.end())
    {
      m_msg->msg(Messenger::ERROR,"Time step for the integrator has not been set.");
      throw runtime_error("Integrator time step not specified.");
    }
    m_dt = lexical_cast<double>(param["dt"]);
  }
  
  //! Propagate system for a time step
  virtual void integrate() = 0;
  
protected:
  
  SystemPtr m_system;          //!< Pointer to the System object
  MessengerPtr m_msg;          //!< Pointer to the messenger object
  PotentialPtr m_potential;    //!< Pointer to the interaction handler 
  AlignerPtr m_align;          //!< Pointer to alignment handler
  NeighbourListPtr m_nlist;    //!< Pointer to the neighbour list object
  ConstraintPtr m_constraint;  //!< Pointer to the handler for constraints
  double m_dt;                 //!< time step
  
};

typedef shared_ptr<Integrator> IntegratorPtr;

#endif

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
 * \file integrator_nve.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Feb-2014
 * \brief Declaration of IntegratorNVE class
 */ 

#ifndef __INTEGRATOR_NVE_H__
#define __INTEGRATOR_NVE_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorNVE class handles NVE dynamics. 
 *  \note No activity. Pure MD.
*/
class IntegratorNVE : public Integrator
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
  IntegratorNVE(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstraintPtr cons, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, param)
  { 
    if (param.find("limit") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"NVE integrator has limit on the maximum particle displacement. It is set to "+param["limit"]+".");
      m_limit = lexical_cast<double>(param["limit"]);;
      m_has_limit = true;
    }
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:

  double  m_limit;        //!< If set, maximum amount particle is allowed to be displaced
  bool    m_has_limit;    //!< Flag that determines if maximum particle displacement has been set
  
};

typedef shared_ptr<IntegratorNVE> IntegratorNVEPtr;

#endif

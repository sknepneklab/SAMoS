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
 * \file integrator_brownian.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 01-Nov-2013
 * \brief Declaration of IntegratorBrownian class
 */ 

#ifndef __INTEGRATOR_BROWNIAN_H__
#define __INTEGRATOR_BROWNIAN_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorBrownian class handles Brownian dynamics.
*/
class IntegratorBrownian : public Integrator
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
  IntegratorBrownian(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstraintPtr cons, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, param)
  { 
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamic integrator. Active velocity v0 not specified. Using default value 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamic integrator. Setting magnitude of active velocity to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
    }
    if (param.find("nu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamic integrator. Rotational diffusion rate not set. Using default value 1.");
      m_nu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamic integrator. Setting rotational diffusion rate to "+param["nu"]+".");
      m_nu = lexical_cast<double>(param["nu"]);
    }
    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamic integrator. Mobility not set. Using default value 1.");
      m_mu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamic integrator. Setting mobility to "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
    }
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamic integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamic integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
    }
    m_stoch_coeff = sqrt(m_nu*m_dt);
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_v0;           //!< Magnitude of the active velocity 
  double  m_nu;           //!< Rotational diffusion 
  double  m_mu;           //!< Mobility 
  double  m_stoch_coeff;  //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  
};

typedef shared_ptr<IntegratorBrownian> IntegratorBrownianPtr;

#endif

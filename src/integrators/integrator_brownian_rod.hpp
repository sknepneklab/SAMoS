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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file integrator_brownian_rod.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 27-Mar-2015
 * \brief Declaration of IntegratorBrownianRod class
 */ 

#ifndef __INTEGRATOR_BROWNIAN_ROD_H__
#define __INTEGRATOR_BROWNIAN_ROD_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;
using std::log;

/*! IntegratorBrownianRod class handles Brownian dynamics of rods.
*/
class IntegratorBrownianRod : public Integrator
{
public:
  
  //! Constructor
  //! \param sys Pointer to a System object containing all particles
  //! \param msg Internal message handler
  //! \param pot Pairwise and external interaction handler
  //! \param align Pairwise and external alignment handler
  //! \param nlist Neighbour list object
  //! \param cons Enforces constraints to the manifold surface
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  IntegratorBrownianRod(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstraintPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. Active velocity v0 not specified. Using default value 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Setting magnitude of active velocity to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
    }
    if (param.find("D0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. Bare diffusion constant D0 not set. Using default value 1.");
      m_D0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Setting bare diffusion constant to "+param["D0"]+".");
      m_D0 = lexical_cast<double>(param["D0"]);
    }
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
    }
    if (param.find("nematic") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. Assuming polar order parameter.");
      m_nematic = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Assuming nematic order parameter.");
      m_nematic = true;
      if (param.find("tau") == param.end())
      {
        m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. Nematic systems No flip rate given. Assuming default 1.");
        m_tau = m_dt/1.0;
      }
      else
      {
        m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Nematic system. Setting flip rate to "+param["tau"]+".");
        m_tau = m_dt/lexical_cast<double>(param["tau"]);
      }
    }
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_v0;           //!< Magnitude of the active velocity 
  double  m_D0;           //!< Bare diffusion constant 
  bool    m_nematic;      //!< If true; assume that the system is nematic, and the velocity will switch direction randomly
  double  m_tau;          //!< Time scale for the direction flip for nematic systems (flip with probability dt/tau)
  
};

typedef shared_ptr<IntegratorBrownianRod> IntegratorBrownianRodPtr;

#endif

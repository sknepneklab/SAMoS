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
 * \file integrator_langevin.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Declaration of IntegratorLangevin class
 */ 

#ifndef __INTEGRATOR_LANGEVIN_H__
#define __INTEGRATOR_LANGEVIN_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;
using std::exp;

/*! IntegratorLangevin class implements Langevin dynamics for the particle position.
    Particle director will not be integrated. 
*/
class IntegratorLangevin : public Integrator
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
  IntegratorLangevin(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator. Active velocity v0 not specified. Using default value 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator. Setting magnitude of active velocity to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
    }
    m_msg->write_config("integrator.langevin.v0",lexical_cast<string>(m_v0));
    if (param.find("gamma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator for particle position. Friction constant not set. Using default value 1.");
      m_gamma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator for particle position. Setting friction constant to "+param["gamma"]+".");
      m_gamma = lexical_cast<double>(param["gamma"]);
    }
    m_msg->write_config("integrator.langevin.gamma",lexical_cast<string>(m_gamma));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator for particle position. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.langevin.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator for particle position. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.langevin.seed",param["seed"]);
    }
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_v0;           //!< Magnitude of the active velocity 
  double  m_gamma;        //!< Friction constant 
  
};

typedef shared_ptr<IntegratorLangevin> IntegratorLangevinPtr;

#endif

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
 * \file integrator_brownian_pos.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Declaration of IntegratorBrownianPos class
 */ 

#ifndef __INTEGRATOR_BROWNIAN_POS_H__
#define __INTEGRATOR_BROWNIAN_POS_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorBrownian class implements Brownian dynamics for the particle position.
    Particle director will not be integrated. 
*/
class IntegratorBrownianPos : public Integrator
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
  IntegratorBrownianPos(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("mu");
    m_known_params.push_back("seed");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for brownian_pos integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in brownian_pos integrator.");
    }
    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for particle position. Mobility not set. Using default value 1.");
      m_mu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for particle position. Setting mobility to "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
    }
    m_msg->write_config("integrator.brownian_pos.mu",lexical_cast<string>(m_mu));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for particle position. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.brownian_pos.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for particle position. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.brownian_pos.seed",param["seed"]);
    }
  }
  
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_mu;           //!< Mobility 
  
};

typedef shared_ptr<IntegratorBrownianPos> IntegratorBrownianPosPtr;

#endif

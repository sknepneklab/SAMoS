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
 * \file integrator_actomyo.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 06-Nov-2014
 * \brief Declaration of IntegratorActomyo class
 */ 

#ifndef __INTEGRATOR_ACTOMYO_H__
#define __INTEGRATOR_ACTOMYO_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorActomyo class handles dynamics of actomyosin system.
*/
class IntegratorActomyo : public Integrator
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
  IntegratorActomyo(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstraintPtr cons, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, param)
  { 
    m_msg->msg(Messenger::WARNING,"Using Actomyo dynamics integrator. All particle groups will be ignored. Working only with group \"all\".");
    if (param.find("zeta") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. Zeta not set. Using default value 1.");
      m_zeta = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Setting zeta to "+param["zeta"]+".");
      m_zeta = lexical_cast<double>(param["zeta"]);
    }
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
    }
    if (param.find("f") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. Active force on myosin not set. Using default value 1.");
      m_f_active = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Setting active force on myosin to "+param["f"]+".");
      m_f_active = lexical_cast<double>(param["f"]);
    }
    if (param.find("T") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. Temperature is not set. Using default value 1.");
      m_D = 1.0/m_zeta;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Setting temperature to "+param["T"]+".");
      m_D = lexical_cast<double>(param["T"])/m_zeta;
    }
    m_stoch_coeff = sqrt(4.0*m_D*m_dt); // Check the factor under square root. I think it is 2d, where d is dimension of the system
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_zeta;         //!< "Friction", prefactor for velocity
  double  m_stoch_coeff;  //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  double  m_f_active;     //!< strength of the active force acting on actin filament
  double  m_D;            //!< Diffusion coefficient
  
};

typedef shared_ptr<IntegratorActomyo> IntegratorActomyoPtr;

#endif

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
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  IntegratorActomyo(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
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
    m_msg->write_config("integrator.actomyo.zeta",lexical_cast<string>(m_zeta));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.actomyo.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.actomyo.seed",param["seed"]);
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
    m_msg->write_config("integrator.actomyo.f",lexical_cast<string>(m_f_active));
    if (param.find("actin_type") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. Actin particle type not set. Using default value 1.");
      m_actin_type = 1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Actin particle types set to "+param["actin_type"]+".");
      m_actin_type = lexical_cast<int>(param["actin_type"]);
    }
    m_msg->write_config("integrator.actomyo.actin_type",lexical_cast<string>(m_actin_type));
    if (param.find("myosin_site") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Actomyo dynamics integrator. Active myosin site particle type not set. Using default value 3.");
      m_mysoin_site_type = 3;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Actomyo dynamics integrator. Active myosin site particle types set to "+param["myosin_site"]+".");
      m_mysoin_site_type = lexical_cast<int>(param["myosin_site"]);
    }
    m_msg->write_config("integrator.actomyo.myosin_site",lexical_cast<string>(m_actin_type));
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;               //!< Random number generator 
  double  m_zeta;              //!< "Friction", prefactor for velocity
  double  m_stoch_coeff;       //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  double  m_f_active;          //!< strength of the active force acting on actin filament
  double  m_D;                 //!< Diffusion coefficient
  int     m_actin_type;        //!< Type of the particles representing actin
  int     m_mysoin_site_type;  //!< Type of the particles representing active site on myosin
  
};

typedef shared_ptr<IntegratorActomyo> IntegratorActomyoPtr;

#endif

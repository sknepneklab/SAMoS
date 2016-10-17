/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  IntegratorBrownian(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("v0");
    m_known_params.push_back("nu");
    m_known_params.push_back("mu");
    m_known_params.push_back("mur");
    m_known_params.push_back("seed");
    m_known_params.push_back("nematic");
    m_known_params.push_back("tau");
    m_known_params.push_back("velocity_align");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for brownian integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in brownian integrator.");
    }
    m_msg->msg(Messenger::WARNING,"Brownian dynamics is a legacy integrator keep for backwards compatibility. It will be removed in a later release. Please consider using brawnian_pos and brownian_align instead.");
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Active velocity v0 not specified. Using default value 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Setting magnitude of active velocity to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
    }
    m_msg->write_config("integrator.brownian.v0",lexical_cast<string>(m_v0));
    if (param.find("nu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Rotational diffusion rate not set. Using default value 1.");
      m_nu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Setting rotational diffusion rate to "+param["nu"]+".");
      m_nu = lexical_cast<double>(param["nu"]);
    }
    m_msg->write_config("integrator.brownian.nu",lexical_cast<string>(m_nu));
    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Mobility not set. Using default value 1.");
      m_mu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Setting mobility to "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
    }
    m_msg->write_config("integrator.brownian.mu",lexical_cast<string>(m_mu));
    if (param.find("mur") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Rotational mobility not set. Using default value 1.");
      m_mur = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Setting rotational mobility to "+param["mur"]+".");
      m_mur = lexical_cast<double>(param["mur"]);
    }
    m_msg->write_config("integrator.brownian.mur",lexical_cast<string>(m_mur));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.brownian.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.brownian.seed",param["seed"]);
    }
    if (param.find("nematic") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Assuming polar order parameter.");
      m_nematic = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Assuming nematic order parameter.");
      m_nematic = true;
      m_msg->write_config("integrator.brownian.nematic","true");
      if (param.find("tau") == param.end())
      {
        m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Nematic systems No flip rate given. Assuming default 1.");
        m_tau = m_dt/1.0;
      }
      else
      {
        m_msg->msg(Messenger::INFO,"Brownian dynamics integrator. Nematic system. Setting flip rate to "+param["tau"]+".");
        m_tau = m_dt/lexical_cast<double>(param["tau"]);
      }
      m_msg->write_config("integrator.brownian.tau",lexical_cast<string>(m_tau));
    }
    if (param.find("velocity_align") == param.end())
    {
      m_velocity = false;
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator. Assuming velocity alignment.");
      m_velocity = true;
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
  double  m_mur;          //!< Rotational mobility
  double  m_stoch_coeff;  //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  bool    m_nematic;      //!< If true; assume that the system is nematic, and the velocity will switch direction randomly
  double  m_tau;          //!< Time scale for the direction flip for nematic systems (flip with probability dt/tau)
  bool    m_velocity;     //!< If true, apply torque to velocity (this is used in simulations with velocity alignmant)
  
};

typedef shared_ptr<IntegratorBrownian> IntegratorBrownianPtr;

#endif

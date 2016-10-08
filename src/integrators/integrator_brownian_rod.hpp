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
  IntegratorBrownianRod(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param), m_pos_noise(false)
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
    m_msg->write_config("integrator.brownian_rod.v0",lexical_cast<string>(m_v0));
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
    m_msg->write_config("integrator.brownian_rod.D0",lexical_cast<string>(m_D0));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.brownian_rod.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.brownian_rod.seed",param["seed"]);
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
      m_msg->write_config("integrator.brownian_rod.nematic","true");
      if (param.find("tau") == param.end())
      {
        m_msg->msg(Messenger::WARNING,"Brownian rod dynamics integrator. Nematic systems No flip rate given. Assuming default 1.");
        m_tau = m_dt/1.0;
        m_msg->write_config("integrator.brownian_rod.tau",lexical_cast<string>(1));
      }
      else
      {
        m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Nematic system. Setting flip rate to "+param["tau"]+".");
        m_tau = m_dt/lexical_cast<double>(param["tau"]);
      }
       m_msg->write_config("integrator.brownian_rod.tau",param["tau"]);
    }
    if (param.find("positional_noise") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Brownian rod dynamics integrator. Adding positional random noise.");
      m_pos_noise = true;
      m_msg->write_config("integrator.brownian_rod.positional_noise","true");
    }
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_v0;           //!< Magnitude of the active velocity 
  double  m_D0;           //!< Bare diffusion constant 
  bool    m_nematic;      //!< If true; assume that the system is nematic, and the velocity will switch direction randomly
  bool    m_pos_noise;    //!< If true, add noise to the rod's position
  double  m_tau;          //!< Time scale for the direction flip for nematic systems (flip with probability dt/tau)
  
};

typedef shared_ptr<IntegratorBrownianRod> IntegratorBrownianRodPtr;

#endif

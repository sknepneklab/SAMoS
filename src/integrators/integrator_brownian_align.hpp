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
 * \file integrator_brownian_align.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Declaration of IntegratorBrownianAlign class
 */ 

#ifndef __INTEGRATOR_BROWNIAN_ALIGN_H__
#define __INTEGRATOR_BROWNIAN_ALIGN_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorBrownianAlign class handles Brownian dynamics.
*/
class IntegratorBrownianAlign : public Integrator
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
  IntegratorBrownianAlign(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("nu");
    m_known_params.push_back("mur");
    m_known_params.push_back("seed");
    m_known_params.push_back("nematic");
    m_known_params.push_back("tau");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for brownian_align integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in brownian_align integrator.");
    }
    if (param.find("nu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for alignment. Rotational diffusion rate not set. Using default value 1.");
      m_nu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for alignment. Setting rotational diffusion rate to "+param["nu"]+".");
      m_nu = lexical_cast<double>(param["nu"]);
    }
    m_msg->write_config("integrator.brownian_align.nu",lexical_cast<string>(m_nu));
    if (param.find("mur") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for alignment. Rotational mobility not set. Using default value 1.");
      m_mur = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for alignment. Setting rotational mobility to "+param["mur"]+".");
      m_mur = lexical_cast<double>(param["mur"]);
    }
    m_msg->write_config("integrator.brownian_align.mur",lexical_cast<string>(m_mur));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for alignment. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.brownian_align.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for alignment. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.brownian_align.seed",param["seed"]);
    }
    if (param.find("nematic") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for alignment. Assuming polar order parameter.");
      m_nematic = false;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for alignment. Assuming nematic order parameter.");
      m_nematic = true;
      m_msg->write_config("integrator.brownian_align.nematic","true");
      if (param.find("tau") == param.end())
      {
        m_msg->msg(Messenger::WARNING,"Brownian dynamics integrator for alignment. Nematic systems No flip rate given. Assuming default 1.");
        m_tau = m_dt/1.0;
      }
      else
      {
        m_msg->msg(Messenger::INFO,"Brownian dynamics integrator for alignment. Nematic system. Setting flip rate to "+param["tau"]+".");
        m_tau = m_dt/lexical_cast<double>(param["tau"]);
      }
      m_msg->write_config("integrator.brownian_align.tau",lexical_cast<string>(m_tau));
    }
    m_stoch_coeff = sqrt(m_nu*m_dt);
  }
  
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_nu;           //!< Rotational diffusion 
  double  m_mur;          //!< Rotational mobility
  double  m_stoch_coeff;  //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  bool    m_nematic;      //!< If true; assume that the system is nematic, and the velocity will switch direction randomly
  double  m_tau;          //!< Time scale for the direction flip for nematic systems (flip with probability dt/tau)
  
};

typedef shared_ptr<IntegratorBrownianAlign> IntegratorBrownianAlignPtr;

#endif

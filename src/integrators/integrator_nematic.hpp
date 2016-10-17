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
 * \file integrator_nematic.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 01-Aug-2014
 * \brief Declaration of IntegratorNematic class
 */ 

#ifndef __INTEGRATOR_NEMATIC_H__
#define __INTEGRATOR_NEMATIC_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorNematic class handles dynamics the active nematic (spatial part) as defined 
 *  in Eq. (2) of E. Bertin, et al. arXiv:1305.0772
 *  Alignment part follows the standard XY-model-like dynamics
*/
class IntegratorNematic : public Integrator
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
  IntegratorNematic(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    if (param.find("d0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Nematic dynamic integrator. Elementary displacement d0 not set. Using default value 0.1.");
      m_d0 = 0.1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Nematic dynamic integrator. Setting elementary displacement d0 to "+param["d0"]+".");
      m_d0 = lexical_cast<double>(param["d0"]);
    }
    m_msg->write_config("integrator.nematic.d0",lexical_cast<string>(m_d0));
    if (param.find("nu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Nematic dynamic integrator. Rotational diffusion rate not set. Using default value 1.");
      m_nu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Nematic dynamic integrator. Setting rotational diffusion rate to "+param["nu"]+".");
      m_nu = lexical_cast<double>(param["nu"]);
    }
    m_msg->write_config("integrator.nematic.nu",lexical_cast<string>(m_nu));
    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Nematic dynamic integrator. Mobility not set. Using default value 1.");
      m_mu = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Nematic dynamic integrator. Setting mobility to "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
    }
    m_msg->write_config("integrator.nematic.mu",lexical_cast<string>(m_mu));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Nematic dynamic integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.nematic.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Nematic dynamic integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.nematic.seed",param["seed"]);
    }
    m_stoch_coeff = sqrt(m_nu*m_dt);
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_d0;           //!< Elementary displacement
  double  m_nu;           //!< Rotational diffusion
  double  m_mu;           //!< Mobility
  double  m_stoch_coeff;  //!< Factor for the stochastic part of the equation of motion (\f$ = \nu \sqrt{dt} \f$)
  
};

typedef shared_ptr<IntegratorNematic> IntegratorNematicPtr;

#endif

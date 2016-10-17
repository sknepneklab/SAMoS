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
 * \file integrator_vicsek.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 10-Dec-2013
 * \brief Declaration of IntegratorVicsek class
 */ 

#ifndef __INTEGRATOR_VICSEK_H__
#define __INTEGRATOR_VICSEK_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;

/*! IntegratorVicsek class handles dynamics of Vicsek model.
*/
class IntegratorVicsek : public Integrator
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
  IntegratorVicsek(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    if (param.find("eta") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Vicsek dynamic integrator. Width of the noise distribution not set. Using default value 0.1.");
      m_eta = 0.1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Vicsek dynamic integrator. Setting random noise distribution width to "+param["eta"]+".");
      m_eta = lexical_cast<double>(param["eta"]);
    }
    m_msg->write_config("integrator.vicsek.eta",lexical_cast<string>(m_eta));
    if (param.find("mu") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Vicsek dynamic integrator. Mobility not set. Using default value 0.");
      m_mu = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Vicsek dynamic integrator. Setting mobility to "+param["mu"]+".");
      m_mu = lexical_cast<double>(param["mu"]);
      if (m_mu != 0.0)
        m_msg->msg(Messenger::WARNING,"Original Vicsek dynamics has no pairwise interaction term beyond the alignment. Using modified approach");
    }
    m_msg->write_config("integrator.vicsek.mu",lexical_cast<string>(m_mu));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Vicsek dynamic integrator. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.vicsek.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Vicsek dynamic integrator. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.vicsek.seed",param["seed"]);
    }
    if (param.find("v0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No velocity magnitude (v0) specified for Vicsek integrator. Setting it to 1.");
      m_v0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Velocity magnitude (v0) for Vicsek integrator is set to "+param["v0"]+".");
      m_v0 = lexical_cast<double>(param["v0"]);
    }
    m_msg->write_config("integrator.vicsek.v0",lexical_cast<string>(m_v0));
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator 
  double  m_eta;          //!< Random noise distribution width
  double  m_mu;           //!< Mobility 
  double  m_v0;       //!< velocity magnitude
  
};

typedef shared_ptr<IntegratorVicsek> IntegratorVicsekPtr;

#endif

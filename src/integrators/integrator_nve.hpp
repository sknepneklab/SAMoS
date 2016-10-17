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
 * \file integrator_nve.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 20-Feb-2014
 * \brief Declaration of IntegratorNVE class
 */ 

#ifndef __INTEGRATOR_NVE_H__
#define __INTEGRATOR_NVE_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

#define SIGN(X) ( (X) >= 0.0 ? 1.0 : -1.0 )

using std::sqrt;
using std::fabs;

/*! IntegratorNVE class handles NVE dynamics. 
 *  \note No activity. Pure MD.
*/
class IntegratorNVE : public Integrator
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
  IntegratorNVE(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_msg->write_config("integrator.nve","");
    if (param.find("limit") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"NVE integrator has limit on the maximum particle displacement. It is set to "+param["limit"]+".");
      m_limit = lexical_cast<double>(param["limit"]);;
      m_has_limit = true;
      m_msg->write_config("integrator.nve.limit",lexical_cast<string>(m_limit));
    }
    if (param.find("angle_limit") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"NVE integrator has limit on the maximum rotation of the director. It is set to "+param["angle_limit"]+".");
      m_theta_limit = lexical_cast<double>(param["angle_limit"]);;
      m_has_theta_limit = true;
      m_msg->write_config("integrator.nve.angle_limit",lexical_cast<string>(m_theta_limit));
    }
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:

  double  m_limit;              //!< If set, maximum amount particle is allowed to be displaced
  double  m_theta_limit;        //!< If set, maximum angular displacement of the director
  bool    m_has_limit;          //!< Flag that determines if maximum particle displacement has been set
  bool    m_has_theta_limit;    //!< Flag that determines if maximum angular displacement of the particle has been set
  
};

typedef shared_ptr<IntegratorNVE> IntegratorNVEPtr;

#endif

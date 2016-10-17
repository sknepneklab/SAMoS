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
 * \file log_vec_velocity.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Mar-2014
 * \brief Declaration of LogVecVelocity class
 */ 

#ifndef __LOG_VEC_VELOCITY_H__
#define __LOG_VEC_VELOCITY_H__

#include <cmath>

#include "log.hpp"

using std::sqrt;

//! LogVecVelocity class
/*! Logs magnitude of total velocity calculated as
 *  \f$ v =  \frac{1}{N}\left|\sum_i \vec v_i\right| \f$. 
 *  It is related to the Vicsek order parameters by dividing this 
 *  quantity by \f$ v_0 \f$.
 */
class LogVecVelocity : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogVecVelocity(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogVecVelocity() { }
  
   //! \return current time step
  string operator()()
  {
    double vx = 0.0, vy = 0.0, vz = 0.0;
    for (int i = 0; i < m_system->size(); i++)
    {
      Particle& p = m_system->get_particle(i);
      vx += p.vx;
      vy += p.vy;
      vz += p.vz;
    }
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    return str(format("%12.6e ") % (v/m_system->size()));
  }
  
};

typedef shared_ptr<LogVecVelocity> LogVecVelocityPtr;

#endif

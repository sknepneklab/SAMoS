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
 * \file log_kinetic_eng.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Jul-2016
 * \brief Declaration of LogKineticEng class
 */ 

#ifndef __LOG_KINETIC_ENG_H__
#define __LOG_KINETIC_ENG_H__

#include <cmath>

#include "log.hpp"

using std::sqrt;

//! LogKineticEng class
/*! Logs kinetic energy per particle for a give time step. 
 */
class LogKineticEng : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogKineticEng(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogKineticEng() { }
  
   //! \return current time step
  string operator()()
  {
    double K = 0.0;
    for (int i = 0; i < m_system->size(); i++)
    {
      Particle& p = m_system->get_particle(i);
      double v2 = p.vx*p.vx + p.vy*p.vy + p.vz*p.vz;
      K += 0.5*v2/p.mass;
    }
    return str(format("%12.6e ") % (K/m_system->size()));
  }
  
};

typedef shared_ptr<LogKineticEng> LogKineticEngPtr;

#endif

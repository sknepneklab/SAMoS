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
 * \file log_perim.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Dec-2015
 * \brief Declaration of LogPerim class
 */ 

#ifndef __LOG_PERIM_H__
#define __LOG_PERIM_H__

#include "log.hpp"


//! LogPerim class
/*! Logs average perimeter for a cell system
 */
class LogPerim : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  LogPerim(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align) : Log(sys, msg, pot, align)  { }
  
  virtual ~LogPerim() { }
  
   //! \return current time step
  string operator()()
  {
    return str(format("%10.8f ") % m_system->compute_average_perimeter());
  }
  
};

typedef shared_ptr<LogPerim> LogPerimPtr;

#endif

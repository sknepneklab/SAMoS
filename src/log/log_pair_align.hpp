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
 * \file log_pair_align.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2014
 * \brief Declaration of LogPairAlign class
 */ 

#ifndef __LOG_PAIR_ALIGN_H__
#define __LOG_PAIR_ALIGN_H__


#include "log.hpp"

//! LogPairAlign class
/*! Logs different types of pair energies
 */
class LogPairAlign : public Log
{
public:
  
  //! Construct Log object
  //! \param sys pointer to a system object
  //! \param compute pointer to a compute object
  //! \param pot Interaction handler object
  //! \param align Alignment handler object
  //! \param type Type of the pair alignment to log (should match the pair_align type)
  LogPairAlign(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, const string& type) : Log(sys, msg, pot, align), m_type(type)  { }
  
  virtual ~LogPairAlign() { }
  
  //! \return current time step
  string operator()()
  {
    return str(format("%12.6e ") % m_aligner->compute_pair_alignment_energy_of_type(m_type));
  }
  
private:
  
  string m_type;    //!< Pair energy type to log
  
};

typedef shared_ptr<LogPairAlign> LogPairAlignPtr;

#endif

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
 * \file parse_timestep.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Jan-2016
 * \brief Grammar for the parsing command that sets the global integrator time step size
 */ 

#ifndef __PARSE_TIMESTEP_HPP__
#define __PARSE_TIMESTEP_HPP__

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct TimeStepData
{
  double dt;    //!< contains integrator step size
};

/*! This is a parser for parsing command that contain number of 
 *  global step size for all integrators
 * 
 *  For example:
 * 
 *  timestep 0.001
 * 
 * This parser will return the size of the integrator step (0.001 in this case)
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class timestep_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  timestep_grammar(TimeStepData& timestep_data) : timestep_grammar::base_type(timestep)
  {
    timestep =    qi::double_[phx::bind(&TimeStepData::dt, phx::ref(timestep_data)) = qi::_1 ]
              >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> timestep;  //!< Rule for parsing timestep lines.
  
};

#endif

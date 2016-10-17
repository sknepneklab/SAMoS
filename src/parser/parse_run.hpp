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
 * \file parse_run.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing command that runs the simulation for a given number of steps
 */ 

#ifndef __PARSE_RUN_HPP__
#define __PARSE_RUN_HPP__

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct RunData
{
  int steps;    //!< contains number of steps to run simulation for
};

/*! This is a parser for parsing command that contain number of 
 *  simulation steps to run the simulation for.
 * 
 *  For example:
 * 
 *  run 100000
 * 
 * This parser will return number of steps (100000 in this case)
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class run_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  run_grammar(RunData& run_data) : run_grammar::base_type(run)
  {
    run =    qi::int_[phx::bind(&RunData::steps, phx::ref(run_data)) = qi::_1 ]
          >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> run;  //!< Rule for parsing run lines.
  
};

#endif

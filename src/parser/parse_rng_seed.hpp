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
 * \file parse_rng_seed.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Grammar for the parsing random number generator seed
 */ 

#ifndef __PARSE_RNG_SEED_HPP__
#define __PARSE_RNG_SEED_HPP__

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct RngSeedData
{
  int seed;    //!< contains random number generator seed
};

/*! This is a parser for parsing command that contain random number 
 *  generator seed.
 * 
 *  For example:
 * 
 *  seed 11
 * 
 * This parser will return the seed (11 in this case)
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class rng_seed_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  rng_seed_grammar(RngSeedData& rng_seed_data) : potential_grammar::base_type(rng_seed)
  {
    rng_seed =   qi::int_[phx::bind(&RngSeedData::seed, phx::ref(rng_seed_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> rng_seed;  //!< Rule for parsing integrator lines.
  
};

#endif

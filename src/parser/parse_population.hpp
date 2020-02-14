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
 * \file parse_population.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Grammar for the parsing populations
 */ 

#ifndef __PARSE_POPULATION_HPP__
#define __PARSE_POPULATION_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct PopulationData
{
  std::string type;      //!< population type (such as random)
  std::string params;    //!< global population parameters 
};

/*! This is a parser for parsing command that contain population control.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the population line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  population random { group=all; division_rate=2.0;  death_rate = 0.5; freq = 100 }
 * 
 * This parser will extract the population type (in this case "random")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class population_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  population_grammar(PopulationData& population_data) : population_grammar::base_type(population)
  {
    population = (
                    qi::as_string[keyword["random"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]      /*! Handles random population control */
                  | qi::as_string[keyword["density"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]     /*! Handles density population control */
                  | qi::as_string[keyword["grow"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]        /*! Handles particle growth */
                  | qi::as_string[keyword["elongate"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]    /*! Handles rod elongation population control */
                  | qi::as_string[keyword["cell"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]        /*! Handles cell growth/division/death */
                  | qi::as_string[keyword["actomyosin"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]  /*! Handles actomyosin attachment/detachment */
                  | qi::as_string[keyword["actomyosin_poisson"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]  /*! Handles actomyosin Poisson attachment/detachment */
                  | qi::as_string[keyword["actomyosin_molecule"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]  /*! Handles actomyosin molecule detachment */
                  | qi::as_string[keyword["actomyosin_head"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]  /*! Handles actomyosin head detachment */
                  | qi::as_string[keyword["region"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]      /*! Handles region population control */
                  /* to add new population: | qi::as_string[keyword["newpopulation"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&PopulationData::params, phx::ref(population_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> population;  //!< Rule for parsing pair population lines.
  
};





#endif

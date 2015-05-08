/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

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
                    qi::as_string[keyword["random"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]         /*! Handles random population control */
                  | qi::as_string[keyword["density"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ]        /*! Handles density population control */
                  /* to add new population: | qi::as_string[keyword["newpopulation"]][phx::bind(&PopulationData::type, phx::ref(population_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&PopulationData::params, phx::ref(population_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> population;  //!< Rule for parsing pair population lines.
  
};





#endif
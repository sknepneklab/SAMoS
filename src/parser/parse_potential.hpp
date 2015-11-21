/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file parse_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Oct-2013
 * \brief Grammar for the parsing pair potentials
 */ 

#ifndef __PARSE_POTENTIAL_HPP__
#define __PARSE_POTENTIAL_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct PotentialData
{
  std::string type;      //!< potential type (such as lj)
  std::string params;    //!< global potential parameters 
};

/*! This is a parser for parsing command that contain pair potentials.
 *  Structurally, it is very similar to the command parsed (\see parse_command.hpp)
 *  This parser extracts the potential line type and returns it together with 
 *  the list of parameters that will be parsed separately.
 * 
 *  for example:
 * 
 *  pair_potential lj { eps = 1.0; sigma = 1.0; r_cut = 2.5; }
 * 
 * This parser will extract the potential type (in this case "lj")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class potential_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  potential_grammar(PotentialData& potential_data) : potential_grammar::base_type(potential)
  {
    potential = (
                  qi::as_string[keyword["lj"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]         /*! Handles Lennard-Jones potential */
                  | qi::as_string[keyword["coulomb"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]  /*! Handles Coulomb potential */
                  | qi::as_string[keyword["debye"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Debye-Hueckel potential */
                  | qi::as_string[keyword["morse"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Morse potential */
                  | qi::as_string[keyword["soft"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]     /*! Handles soft-core potential */
                  | qi::as_string[keyword["gaussian"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] /*! Handles Gaussian potential */
                  | qi::as_string[keyword["morse"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Morse potential */
                  | qi::as_string[keyword["active"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]   /*! Handles active potential */
                  | qi::as_string[keyword["rod"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]      /*! Handles soft rod potential */
                  | qi::as_string[keyword["ljrod"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles Lennard-Jones rod potential */
                  | qi::as_string[keyword["soft_attractive"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]    /*! Handles soft attractive potential */
                  | qi::as_string[keyword["vp"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ]       /*! Handles vertex-particle potential */
                  /* to add new potential: | qi::as_string[keyword["newpotential"]][phx::bind(&PotentialData::type, phx::ref(potential_data)) = qi::_1 ] */
                 )
                 >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&PotentialData::params, phx::ref(potential_data)) = qi::_1 ]
                 >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> potential;  //!< Rule for parsing pair potential lines.
  
};





#endif
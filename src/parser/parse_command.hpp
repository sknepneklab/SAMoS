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
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file parse_command.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 21-Oct-2013
 * \brief Grammar for the command parser
 */ 

#ifndef __PARSE_COMMAND_HPP__
#define __PARSE_COMMAND_HPP__

#include <string>

#include "parse_aux.hpp"

/*! Control data structure for passing around parsed information. */
struct CommandData
{
  std::string command;                //!< contain command
  std::string attrib_param_complex;   //!< attributes and parameters (i.e., the rest of the line)
};

/*! This is a simple parser for parsing commands in the control 
 *  file of the code. The parser is based on Boost Spirit Qi library.
 *  This parser only extracts key word command from a single line in the
 *  control file. 
 *  Syntax of the control file is:
 * 
 *  command [attribute] { [parameters] }
 * 
 *  for example:
 * 
 *  pair_potential lj { eps = 1.0; sigma = 1.0; r_cut = 2.5; }
 * 
 * This parser will only extract the keyword (in this case "pair_potential")
 * and return the rest of the line for post-processing.
 * 
 * \note Based on http://www.boost.org/doc/libs/1_53_0/libs/spirit/repository/test/qi/distinct.cpp
*/
class command_grammar : public qi::grammar<std::string::iterator, qi::space_type>
{
public:
  command_grammar(CommandData& command_data) : command_grammar::base_type(command)
  {
    command = (
                 qi::as_string[keyword["pair_potential"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]    /*! Handles pair potential, such as Lennard-Jones or Coulomb */
                 | qi::as_string[keyword["pair_param"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]      /*! Handles parameters for the pair potentials */
                 | qi::as_string[keyword["align_param"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]     /*! Handles parameters for the alignment */
                 | qi::as_string[keyword["external"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]        /*! Handles external potentials */
                 | qi::as_string[keyword["external_param"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]  /*! Handles parameters for external potentials */
                 | qi::as_string[keyword["integrator"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]      /*! Handles integrators (NVE, NTV, etc.) */
                 | qi::as_string[keyword["input"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]           /*! Handles input file. */
                 | qi::as_string[keyword["constraint"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]      /*! Handles constraints. */
                 | qi::as_string[keyword["run"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]             /*! Execute simulation. */
                 | qi::as_string[keyword["log"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]             /*! Handles various measurement logs. */
                 | qi::as_string[keyword["dump"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]            /*! Handles coordinate dumps. */
                 | qi::as_string[keyword["messages"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]        /*! Handles internal messages. */
                 | qi::as_string[keyword["box"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]             /*! Handles simulation box. */
                 | qi::as_string[keyword["nlist"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]           /*! Handles neighbour list. */
                 | qi::as_string[keyword["pair_align"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]      /*! Handles pairwise alignment. */
                 | qi::as_string[keyword["external_align"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]  /*! Handles external alignment. */
                 | qi::as_string[keyword["external_align_param"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]  /*! Handles external alignment parameters. */
                 | qi::as_string[keyword["group"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]           /*! Handles particle groups. */
                 | qi::as_string[keyword["disable"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]         /*! Handles disabling integrators. */
                 | qi::as_string[keyword["population"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ]      /*! Handles population control. */
                 /* to add new command: | qi::as_string[keyword["newcommand"]][ phoenix::bind(&CommandData::command, phoenix::ref(command_data)) = qi::_1 ] */
               )
               >> qi::as_string[qi::no_skip[+qi::char_]][phoenix::bind(&CommandData::attrib_param_complex, phoenix::ref(command_data)) = qi::_1 ]
               >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> command;  //!< Rule for parsing command lines.
  
};





#endif
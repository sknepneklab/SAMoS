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
                 qi::as_string[keyword["pair_potential"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]    /*! Handles pair potential, such as Lennard-Jones or Coulomb */
                 | qi::as_string[keyword["pair_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles parameters for the pair potentials */
                 | qi::as_string[keyword["align_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]     /*! Handles parameters for the alignment */
                 | qi::as_string[keyword["external"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]        /*! Handles external potentials */
                 | qi::as_string[keyword["external_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]  /*! Handles parameters for external potentials */
                 | qi::as_string[keyword["integrator"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles integrators (NVE, NTV, etc.) */
                 | qi::as_string[keyword["input"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]           /*! Handles input file. */
                 | qi::as_string[keyword["constraint"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles constraints. */
                 | qi::as_string[keyword["run"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]             /*! Execute simulation. */
                 | qi::as_string[keyword["log"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]             /*! Handles various measurement logs. */
                 | qi::as_string[keyword["dump"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]            /*! Handles coordinate dumps. */
                 | qi::as_string[keyword["messages"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]        /*! Handles internal messages. */
                 | qi::as_string[keyword["box"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]             /*! Handles simulation box. */
                 | qi::as_string[keyword["nlist"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]           /*! Handles neighbour list. */
                 | qi::as_string[keyword["pair_align"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles pairwise alignment. */
                 | qi::as_string[keyword["external_align"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]  /*! Handles external alignment. */
                 | qi::as_string[keyword["external_align_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]  /*! Handles external alignment parameters. */
                 | qi::as_string[keyword["group"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]           /*! Handles particle groups. */
                 | qi::as_string[keyword["disable"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]         /*! Handles disabling integrators. */
                 | qi::as_string[keyword["population"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles population control. */
                 | qi::as_string[keyword["zero_momentum"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]   /*! Handles zeroing momentum of a group of particles. */
                 | qi::as_string[keyword["bond"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]            /*! Handles bond potentials. */
                 | qi::as_string[keyword["bond_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles bond potential parameters. */
                 | qi::as_string[keyword["angle"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]           /*! Handles angle potentials. */
                 | qi::as_string[keyword["angle_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]     /*! Handles angle potential parameters. */
                 | qi::as_string[keyword["read_bonds"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]      /*! Handles bonds file. */
                 | qi::as_string[keyword["read_angles"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]     /*! Handles angles file. */
                 | qi::as_string[keyword["read_cell_boundary"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]    /*! Handles cell boundary connectivity file. */
                 | qi::as_string[keyword["ntypes"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]          /*! Handles number of types. */
                 | qi::as_string[keyword["config"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]          /*! Handles configuration file. */
                 | qi::as_string[keyword["pair_type_param"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ] /*! Handles particle type parameters for the pair potentials */
                 | qi::as_string[keyword["timestep"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]        /*! Hangles global integrator step. */
                 | qi::as_string[keyword["disable_population"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ]        /*! Hangles global disable population command. */
                 /* to add new command: | qi::as_string[keyword["newcommand"]][ phx::bind(&CommandData::command, phx::ref(command_data)) = qi::_1 ] */
               )
               >> qi::as_string[qi::no_skip[+qi::char_]][phx::bind(&CommandData::attrib_param_complex, phx::ref(command_data)) = qi::_1 ]
               >> (qi::eol || qi::eoi);
  }

private:
    
  qi::rule<std::string::iterator, qi::space_type> command;  //!< Rule for parsing command lines.
  
};





#endif

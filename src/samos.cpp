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
 * \file samos.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 01-Nov-2013
 * \brief Main program
*/

#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <ctime>

#include <boost/functional/factory.hpp>
#include <boost/function.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/make_shared.hpp>

#include "samos.hpp"

#include "factory_types.hpp"
#include "register.hpp"

using std::map;
using std::string;


int main(int argc, char* argv[])
{
  // Parser data
  CommandData             command_data;
  BoxData                 box_data;
  InputData               input_data;
  ExternalData            external_data;
  LogDumpData             log_dump_data;
  PotentialData           potential_data;
  IntegratorData          integrator_data;
  ConstraintlData         constraint_data;
  RunData                 run_data;
  AlignData               pair_align_data;
  ExternalAlignData       external_align_data;
  GroupData               group_data;
  DisableData             disable_data;
  PopulationData          population_data;
  DisablePopulationData   population_disable_data;
  BondData                bond_data;
  AngleData               angle_data;
  TimeStepData            timestep_data;
  pairs_type              parameter_data;   // All parameters for different commands
  
  // Parser grammars
  command_grammar                 command_parser(command_data);
  box_grammar                     box_parser(box_data);
  input_grammar                   input_parser(input_data);
  external_grammar                external_parser(external_data); 
  log_dump_grammar                log_dump_parser(log_dump_data);
  potential_grammar               potential_parser(potential_data);
  integrator_grammar              integrator_parser(integrator_data);  
  constraint_grammar              constraint_parser(constraint_data);
  run_grammar                     run_parser(run_data);
  align_grammar                   align_parser(pair_align_data);
  external_align_grammar          external_align_parser(external_align_data);
  group_grammar                   group_parser(group_data);
  disable_grammar                 disable_parser(disable_data);
  population_grammar              population_parser(population_data);
  disable_population_grammar      disable_population_parser(population_disable_data);
  bond_grammar                    bond_parser(bond_data);
  angle_grammar                   angle_parser(angle_data);
  timestep_grammar                timestep_parser(timestep_data);
  key_value_sequence              param_parser;
  
  // Class factories 
  ConstraintMap constraints;
  PairPotentialMap pair_potentials;
  ExternalPotentialMap external_potentials;
  IntegratorMap integrators;
  PairAlignerMap pair_aligners;
  ExternalAlignerMap external_aligners;
  PopulationMap populations;
  BondPotentialMap bond_potentials;
  AnglePotentialMap angle_potentials;
  ValueMap values;
  
  std::ifstream command_file;    // File with simulation parameters and controls
  std::string command_line;      // Line from the command file
  
  MessengerPtr msg;                                // Messenger object
  BoxPtr box;                                      // Handles simulation box
  SystemPtr sys;                                   // System object
  PotentialPtr pot;                                // Handles all potentials
  ConstrainerPtr constraint;                       // Handles constrints to the manifold
  NeighbourListPtr nlist;                          // Handles global neighbour list
  map<string,IntegratorPtr> integrator;            // Handles the integrator
  AlignerPtr    aligner;                           // Handles all aligners
  ExternalAlignPtr external_aligner;               // Handles all external aligners
  vector<DumpPtr> dump;                            // Handles all different dumps
  vector<LoggerPtr> log;                           // Handles all different logs
  map<string,PopulationPtr>  population;           // Handles all population methods
  
  bool periodic = false;        // If true, use periodic boundary conditions
  bool has_potential = false;   // If false, potential handling object (Potential class) has not be initialized yet
  bool has_constraints = false; // If false, constraint handling object (Constainer class) has not be initialized yet
  bool has_aligner = false;     // If false, pairwise aligners are not defined
  bool has_population = false;  // If false, population is not defined
  int time_step = 0;            // Counts current time step
  
  // flags that ensure that system the simulation is not in a bad state
  std::map<std::string, bool>  defined;
  defined["messages"] = false;         // If false, system has no messenger (set it to default "messages.msg")
  defined["box"] = false;              // If false, no simulation box defined (cannot run simulation)
  defined["input"] = false;            // If false, no system object defines (cannot run the simulation)
  defined["pair_potential"] = false;   // If false, no pair potentials have been specified (needs at least one external potential)
  defined["external"] = false;         // If false, no external potentials have been specified (needs at least one pair potential)
  defined["constraint"] = false;       // If false, no constraints have been defined (a full 3d simulation with a rather slow neighbour list building algorithm)
  defined["nlist"] = false;            // If false, system has no neighbour list 
  defined["integrator"] = false;       // If false, no integrator has been defined (cannot run simulation)
  defined["pair_aligner"] = false;     // If false, no pairwise alignment has been defined
  defined["external_aligner"] = false; // If false, no external alignment has been defined
  defined["bond"] = false;             // If false, no bonds have been added to the system
  defined["angle"] = false;            // If false, no angles have been defined in the system
  defined["read_bonds"] = false;       // If false, no bond potentials have been defined
  defined["read_angles"] = false;      // If false, no angle potentials have been defined
  defined["population"] = false;       // If false, no populations have been defined
    
  register_constraints(constraints);                     // Register all constraints
  register_pair_potentials(pair_potentials);             // Register all pair potentials 
  register_external_potentials(external_potentials);     // Register all external potentials 
  register_integrators(integrators);                     // Register all integrators
  register_pair_aligners(pair_aligners);                 // Register all pair aligners 
  register_external_aligners(external_aligners);         // Register all external aligners
  register_populations(populations);                     // Register all populations
  register_bond_potentials(bond_potentials);             // Register all bond potentials
  register_angle_potentials(angle_potentials);           // Register all angle potentials
  register_values(values);                               // Register all values
  
  
  if (argc < 2)
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << "    samos <config file name>" << std::endl;
    return -1;
  }
  
  // Measure execution time
  std::time_t start_time = std::time(NULL);
   
  
  // Start parsing the configuration file
  int current_line = 1;
  command_file.open(argv[1],std::ifstream::in);
  if (command_file)
  {
    while (std::getline(command_file, command_line))
    {
      boost::algorithm::trim(command_line);            // get rid of all leading and trailing white spaces 
      //boost::algorithm::to_lower(command_line);        // transform the entire line to lower case (command script is case insensitive)
      // skip empty and comment lines and start parsing real commands
      if (command_line.size() > 0 && command_line[0] != '#')
      {
        parameter_data.clear();
        if (qi::phrase_parse(command_line.begin(), command_line.end(), command_parser, qi::space))
        {
          if (defined.find(command_data.command) != defined.end())
            defined[command_data.command] = true;
          if (command_data.command == "box")               // if command is simulation box, parse it and define it
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), box_parser, qi::space))
            {
              if (box_data.type == "periodic")    periodic = true;
              else if (box_data.type == "fixed")  periodic = false;
              else
              {
                if (defined["messages"])
                  msg->msg(Messenger::ERROR,"Box type "+box_data.type+" is not known.");
                else
                  std::cerr << "Box type " << box_data.type << " is not known." << std::endl;
                throw std::runtime_error("Unknown box type.");
              }
              // Box size
              double lx = DEFAULT_LX, ly = DEFAULT_LY, lz = DEFAULT_LZ;
              // Finally parse box parameters 
              if (qi::phrase_parse(box_data.params.begin(), box_data.params.end(), param_parser, qi::space, parameter_data))
              {
                if (parameter_data.find("lx") != parameter_data.end()) lx = lexical_cast<double>(parameter_data["lx"]);
                if (parameter_data.find("ly") != parameter_data.end()) ly = lexical_cast<double>(parameter_data["ly"]);
                if (parameter_data.find("lz") != parameter_data.end()) lz = lexical_cast<double>(parameter_data["lz"]);
              }
              else
              {
                if (defined["messages"]) msg->msg(Messenger::ERROR,"Error parsing simulation box parameters.");
                else  std::cout << "Error parsing simulation box parameters."  << std::endl;
                throw std::runtime_error("Error parsing box parameters.");
              }
              box = std::make_shared<Box>(Box(lx,ly,lz));
              if (defined["messages"])   msg->msg(Messenger::INFO,"Simulation box is "+box_data.type+" with size (lx,ly,lz) = ("+lexical_cast<string>(lx)+","+lexical_cast<string>(ly)+","+lexical_cast<string>(lz)+").");
              else  std::cout << "Simulation box is "+box_data.type+" with size (lx,ly,lz) = ("+lexical_cast<string>(lx)+","+lexical_cast<string>(ly)+","+lexical_cast<string>(lz)+")." << std::endl;
            }
            else
            {
              std::cerr << "Error parsing box command at line : " << current_line << std::endl;
              throw std::runtime_error("Error parsing command script.");
            }
          }
          else if (command_data.command == "messages")       // if command is messages, handle the global messenger 
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), input_parser, qi::space))  // Note that Messenger uses the same parser as the input files parser
            {
              msg = std::shared_ptr<Messenger>(new Messenger(input_data.name));
              msg->msg(Messenger::INFO,"Messages will be sent to "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse messenger command." << std::endl;
              throw std::runtime_error("Could not parse messenger command.");
            }
          }
          else if (command_data.command == "input")       // if command is input, parse it and read in the system coordinates, if possible
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), input_parser, qi::space))  
            {
              if (!defined["messages"])  // If messenger is not defined, send to default messenger defined at the top of this file
              {
                string msg_name = DEFAULT_MESSENGER;
                msg = std::shared_ptr<Messenger>(new Messenger(msg_name));
                msg->msg(Messenger::WARNING,"Messenger was not defined prior to the reading in data. If not redefined all messages will be sent to "+msg_name+".");
              }
              if (!defined["box"])   // Cannot run without a box
              {
                msg->msg(Messenger::ERROR,"Simulation box has not been defined. Please define it before reading in coordinates.");
                throw std::runtime_error("Simulation box not defined.");
              }
              sys = std::make_shared<System>(System(input_data.name,msg,box));
              sys->set_periodic(periodic);
              msg->msg(Messenger::INFO,"Finished reading system coordinates from "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse input command." << std::endl;
              throw std::runtime_error("Could not parse input command.");
            }
          }
          else if (command_data.command == "read_bonds")       // if command is read_bonds, parse it and read in the bonds information, if possible
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), input_parser, qi::space))  
            {
              if (!defined["input"])  // We need to have system defined before we can add any bond potentials 
              {
                if (defined["messages"])
                  msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before reading bond information.");
                else
                  std::cerr << "System has not been defined. Please define system using \"input\" command before reading bond information." << std::endl;
                throw std::runtime_error("System not defined.");
              }
              sys->read_bonds(input_data.name);
              msg->msg(Messenger::INFO,"Finished reading bond data from "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse read_bonds command." << std::endl;
              throw std::runtime_error("Could not parse read_bonds command.");
            }
          }
          else if (command_data.command == "read_angles")       // if command is read_angles, parse it and read in the angles information, if possible
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), input_parser, qi::space))  
            {
              if (!defined["input"])  // We need to have system defined before we can add any angle potentials
              {
                if (defined["messages"])
                  msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before reading angle information.");
                else
                  std::cerr << "System has not been defined. Please define system using \"input\" command before reading angle information." << std::endl;
                throw std::runtime_error("System not defined.");
              }
              sys->read_angles(input_data.name);
              msg->msg(Messenger::INFO,"Finished reading angle data from "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse read_angles command." << std::endl;
              throw std::runtime_error("Could not parse read_angles command.");
            }
          }
          else if (command_data.command == "read_cell_boundary")       // if command is read_cell_boundary, parse it and read in the cell boundary connectivity information, if possible
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), input_parser, qi::space))  
            {
              if (!defined["input"])  // We need to have system defined before we can add any angle potentials
              {
                if (defined["messages"])
                  msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before reading cell boundary connectivity information.");
                else
                  std::cerr << "System has not been defined. Please define system using \"input\" command before reading cell boundary connectivity infromation." << std::endl;
                throw std::runtime_error("System not defined.");
              }
              sys->read_boundary_neighbours(input_data.name);
              msg->msg(Messenger::INFO,"Finished reading boundary connectivity data from "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse read_cell_boundary command." << std::endl;
              throw std::runtime_error("Could not parse read_cell_boundary command.");
            }
          }
          else if (command_data.command == "pair_potential")       // if command is pair potential, parse it and add this pair potential to the list of pair potentials
          {
            if (!defined["input"])  // We need to have system defined before we can add any external potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any pair potentials.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any pair potentials." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), potential_parser, qi::space))
            {
              if (!has_potential) 
              {
                pot = std::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(potential_data.params.begin(), potential_data.params.end(), param_parser, qi::space, parameter_data))
              {
                if (!defined["nlist"]) // Create a default neighbour list if none has been defined
                {
                  msg->msg(Messenger::WARNING,"No neighbour list defined. Some pair potentials (e.g. Lennard-Jones) need it. We are making one with default cutoff and padding.");
                  nlist = std::make_shared<NeighbourList>(NeighbourList(sys,msg,DEFAULT_CUTOFF,DEFAULT_PADDING,parameter_data));
                  defined["nlist"] = true;
                }
                std::string phase_in = "constant";
                if (parameter_data.find("phase_in") != parameter_data.end())
                    phase_in = parameter_data["phase_in"];                
                pot->add_pair_potential(potential_data.type, pair_potentials[potential_data.type](
                                                                                                  sys,
                                                                                                  msg,
                                                                                                  nlist,
                                                                                                  std::shared_ptr<Value>(values[phase_in](msg,parameter_data)),
                                                                                                  parameter_data
                                                                                                 ));
                msg->msg(Messenger::INFO,"Added "+potential_data.type+" to the list of pair potentials.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse pair potential parameters for potential type "+potential_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing pair potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing pair_potential command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing pair_potential command.");
            }
            
          }
          else if (command_data.command == "external")       // if command is external potential, parse it and add this external potential to the list of external potentials
          {
            if (!defined["input"])  // We need to have system defined before we can add any external potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any external potentials.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any external potentials." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), external_parser, qi::space))
            {
              if (!has_potential) 
              {
                pot = std::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(external_data.params.begin(), external_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_external_potential(external_data.type, external_potentials[external_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Added "+external_data.type+" to the list of external potentials.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse external potential parameters for potential type "+external_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing external potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing external command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing external command.");
            }
          }
          else if (command_data.command == "bond")       // if command is bond, parse it and add this bond potential to the list of bond potentials
          {
            if (!defined["input"])  // We need to have system defined before we can add any external potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any bond potentials.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any angle potentials." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (!defined["bond"])  // We need to have bonds defined before we can add any bond potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"Bond file has not been specified. Please define system using \"read_bonds\" command before adding any bond potentials.");
              else
                std::cerr << "Bonds have not been defined. Please define system using \"read_bonds\" command before adding any bond potentials." << std::endl;
              throw std::runtime_error("Bonds not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), bond_parser, qi::space))
            {
              if (!has_potential) 
              {
                pot = std::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(bond_data.params.begin(), bond_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_bond_potential(bond_data.type, bond_potentials[bond_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Added "+bond_data.type+" to the list of bond potentials.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse bond potential parameters for bond type "+bond_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing bond potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing bond command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing bond command.");
            }
            
          }
          else if (command_data.command == "angle")       // if command is angle, parse it and add this angle potential to the list of angle potentials
          {
            if (!defined["input"])  // We need to have system defined before we can add any external potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any angle potentials.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any angle potentials." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (!defined["angle"])  // We need to have angles defined before we can add any bond potential 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"Angle file has not been specified. Please define system using \"read_angles\" command before adding any angle potentials.");
              else
                std::cerr << "Angles have not been defined. Please define system using \"read_angles\" command before adding any angle potentials." << std::endl;
              throw std::runtime_error("Angles not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), angle_parser, qi::space))
            {
              if (!has_potential) 
              {
                pot = std::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(angle_data.params.begin(), angle_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_angle_potential(angle_data.type, angle_potentials[angle_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Added "+angle_data.type+" to the list of angle potentials.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse angle potential parameters for angle type "+angle_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing angle potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing angle command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing angle command.");
            }
            
          }
          else if (command_data.command == "constraint")       // if command is constraint, parse it and add create appropriate constraint object
          {
            if (!defined["input"])  // we need to have system defined before we can add any constraints
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding constraint.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding constraint." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), constraint_parser, qi::space))
            {
              if (!has_constraints) 
              {
                constraint = std::make_shared<Constrainer>(Constrainer(sys,msg));
                has_constraints = true;
              }
              if (qi::phrase_parse(constraint_data.params.begin(), constraint_data.params.end(), param_parser, qi::space, parameter_data))
              {
                //constraint = std::shared_ptr<Constraint>(constraints[constraint_data.type](sys,msg,parameter_data));  // dirty workaround shared_ptr and inherited classes
                constraint->add_constraint(constraint_data.type,constraints[constraint_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Adding constraint of type "+constraint_data.type+".");
                // Enforce constraint so we make sure all particles lie on it
                for (int i = 0; i < sys->size(); i++)
                  constraint->enforce(sys->get_particle(i));
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for constraint "+constraint_data.type+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing constraint parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse constraint command at line "+lexical_cast<string>(current_line));
              throw std::runtime_error("Error parsing constraint line.");
            }
          }
          else if (command_data.command == "nlist")       // if command is nlist, parse it and create appropriate neighbour list object
          {
            if (!defined["input"])  // we need to have system defined before we can set up neighbour list
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding neighbour list.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding neighbour list." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), param_parser, qi::space, parameter_data))
            {
              double rcut = DEFAULT_CUTOFF;
              double pad = DEFAULT_PADDING;
              if (parameter_data.find("rcut") != parameter_data.end()) rcut = lexical_cast<double>(parameter_data["rcut"]);
              if (parameter_data.find("pad") != parameter_data.end())  pad = lexical_cast<double>(parameter_data["pad"]);
              nlist = std::make_shared<NeighbourList>(NeighbourList(sys,msg,rcut,pad,parameter_data));
              msg->msg(Messenger::INFO,"Created neighbour list with cutoff "+lexical_cast<string>(rcut)+" and padding distance "+lexical_cast<string>(pad)+".");
              defined["nlist"] = true;
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse neighbour list parameters at line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing neighbour list parameters.");
            }
          }
          else if (command_data.command == "dump")   // parse dump commands by adding list of dumps
          {
            if (!defined["input"])  // We need to have system defined before we can add any dumps
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any dump.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any dumps." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), log_dump_parser, qi::space))
            {
              if (qi::phrase_parse(log_dump_data.params.begin(), log_dump_data.params.end(), param_parser, qi::space, parameter_data))
              {
                dump.push_back(std::shared_ptr<Dump>(new Dump(sys,msg,nlist,log_dump_data.name,parameter_data)));
                msg->msg(Messenger::INFO,"Adding dump to file "+log_dump_data.name+".");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  if ((*it).second != "")
                  {
                    msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for dump "+log_dump_data.name+" is set to "+(*it).second+".");
                  }
                  else 
                  {
                    msg->msg(Messenger::INFO,"Setting flag " + (*it).first + " for dump "+log_dump_data.name+".");
                  }
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for dump "+log_dump_data.name+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing dump parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing dump command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing dump line.");
            }
          }
          else if (command_data.command == "log")   // parse dump commands by adding list of logs
          {
            if (!defined["input"])  // We need to have system defined before we can add any logs
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any logs.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any logs." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), log_dump_parser, qi::space))
            {
              if (qi::phrase_parse(log_dump_data.params.begin(), log_dump_data.params.end(), param_parser, qi::space, parameter_data))
              {
                log.push_back(std::shared_ptr<Logger>(new Logger(sys,msg,pot,aligner,log_dump_data.name,parameter_data)));
                msg->msg(Messenger::INFO,"Adding logger to file "+log_dump_data.name+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for log "+log_dump_data.name+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing log parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing log command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing log line.");
            }
          }
          else if (command_data.command == "integrator")   // parse integrator command and construct integrator object
          {
            if (!defined["input"])  // We need to have system defined before we can add the integrator
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any dump.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any dumps." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            /*
            if (!has_potential)
            {
              msg->msg(Messenger::ERROR,"No potentials has been defined. Please define them using \"pair_potential\" or \"external\" command before adding integrator.");
              throw std::runtime_error("No potentials defined.");
            }
            */
            if (!defined["constraint"])
            {
              msg->msg(Messenger::ERROR,"Constraint has not been defined. Please define them using \"constraint\" command before adding integrator.");
              throw std::runtime_error("Constraint not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), integrator_parser, qi::space))
            {
              if (qi::phrase_parse(integrator_data.params.begin(), integrator_data.params.end(), param_parser, qi::space, parameter_data))
              {
                std::string group_name;
                if (parameter_data.find("group") == parameter_data.end())
                  group_name = "all";
                else
                  group_name = parameter_data["group"];
                if (integrator.find(integrator_data.type+"_"+group_name) == integrator.end())
                {
                  std::string temperature_control;
                  if (parameter_data.find("temperature_control") == parameter_data.end())
                    temperature_control = "constant";
                  else
                    temperature_control = parameter_data["temperature_control"];
                  integrator[integrator_data.type+"_"+group_name] = std::shared_ptr<Integrator>(
                                                                                                  integrators[integrator_data.type]
                                                                                                  (
                                                                                                    sys,
                                                                                                    msg,
                                                                                                    pot,
                                                                                                    aligner,
                                                                                                    nlist,
                                                                                                    constraint,
                                                                                                    std::shared_ptr<Value>(values[temperature_control](msg,parameter_data)),
                                                                                                    parameter_data
                                                                                                  )
                                                                                                 );
                  msg->msg(Messenger::INFO,"Adding integrator of type "+integrator_data.type+" (for group "+group_name+").");
                }
                else
                {
                  msg->msg(Messenger::ERROR,"Integrator of type "+integrator_data.type+" (for group "+group_name+") already exists. You need to disable the old one first.");
                  throw std::runtime_error("Integrator already exists.");
                }
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for integrator "+integrator_data.type+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing integrator parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing integrator command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing integrator line.");
            }
          }
          else if (command_data.command == "pair_param")       // parse pair parameters for potential
          {
            if (!defined["pair_potential"])  // We need to have pair potentials defined before we can change their parameters
            {
              msg->msg(Messenger::ERROR,"No pair potentials have been defined. Please define them using \"pair_potential\" command before modifying any pair parameters.");
              throw std::runtime_error("No pair potentials defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), potential_parser, qi::space))
            {
              if (qi::phrase_parse(potential_data.params.begin(), potential_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_pair_potential_parameters(potential_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+potential_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse pair potential parameters for potential type "+potential_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing pair potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing pair_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing pair_param command.");
            }
            
          }
          else if (command_data.command == "pair_type_param")       // parse pair type parameters for potential
          {
            if (!defined["pair_potential"])  // We need to have pair potentials defined before we can change their parameters
            {
              msg->msg(Messenger::ERROR,"No pair potentials have been defined. Please define them using \"pair_potential\" command before modifying any pair type parameters.");
              throw std::runtime_error("No pair potentials defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), potential_parser, qi::space))
            {
              if (qi::phrase_parse(potential_data.params.begin(), potential_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_pair_type_parameters(potential_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+potential_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse pair type parameters for potential type "+potential_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing pair type parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing pair_type_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing pair_type_param command.");
            }
            
          }
          else if (command_data.command == "external_param")       // parse pair parameters for potential
          {
            if (!defined["external"])  // We need to have at least one external potential defined before we can change the parameters
            {
              msg->msg(Messenger::ERROR,"No external potentials have been defined. Please define them using \"external\" command before modifying any parameters.");
              throw std::runtime_error("No pair potentials defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), external_parser, qi::space))
            {
              if (qi::phrase_parse(external_data.params.begin(), external_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_external_potential_parameters(external_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+external_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse external potential parameters for potential type "+external_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing external potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing external_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing external_param command.");
            }
          }
          else if (command_data.command == "bond_param")       // parse bond parameters for a bond potential
          {
            if (!defined["read_bonds"])  // We need to have at least one bond potential defined before we can change the parameters
            {
              msg->msg(Messenger::ERROR,"No bond potentials have been defined. Please define them using \"bond\" command before modifying any parameters.");
              throw std::runtime_error("No bond potentials defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), bond_parser, qi::space))
            {
              if (qi::phrase_parse(bond_data.params.begin(), bond_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_bond_potential_parameters(bond_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+bond_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse bond potential parameters for bond potential type "+bond_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing bond potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing bond_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing bond_param command.");
            }
          }
          else if (command_data.command == "angle_param")       // parse angle parameters for an angle potential
          {
            if (!defined["read_angles"])  // We need to have at least one angle potential defined before we can change the parameters
            {
              msg->msg(Messenger::ERROR,"No angle potentials have been defined. Please define them using \"angle\" command before modifying any parameters.");
              throw std::runtime_error("No angle potentials defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), angle_parser, qi::space))
            {
              if (qi::phrase_parse(angle_data.params.begin(), angle_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_angle_potential_parameters(angle_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+angle_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse angle potential parameters for bond potential type "+angle_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing angle potential parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing angle_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing angle_param command.");
            }
          }
          else if (command_data.command == "run")       // run actual simulation
          {
            if (!defined["input"])  // We need to have system defined before we can run simulation
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any dump.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any dumps." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            /*
            if (!has_potential)
            {
              msg->msg(Messenger::ERROR,"No potentials has been defined. Please define them using \"pair_potential\" or \"external\" command before running simulation.");
              throw std::runtime_error("No potentials defined.");
            }
            */
            if (!defined["constraint"])
            {
              msg->msg(Messenger::ERROR,"Constraint has not been defined. Please define it using \"constraint\" command before running simulation.");
              throw std::runtime_error("Constraint not defined.");
            }
            if (!defined["integrator"])
            {
              msg->msg(Messenger::ERROR,"Integrator has not been defined. Please define it using \"integrator\" command before running simulation.");
              throw std::runtime_error("Integrator not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), run_parser, qi::space))
            {
              msg->msg(Messenger::INFO,"Starting simulation run for "+lexical_cast<string>(run_data.steps)+" steps.");
              for (std::map<std::string, IntegratorPtr>::iterator it_integ = integrator.begin(); it_integ != integrator.end(); it_integ++)
              {
                std::string integrator_name = (*it_integ).first;
                msg->add_config("run."+integrator_name+".steps",lexical_cast<string>(run_data.steps));
              }
              // Precompute forces and torques
              if (pot)
                pot->compute(1e-3);  // Some value to make sure phase in is working.
              if (aligner)
                aligner->compute();
              int nlist_builds = 0;     // Count how many neighbour list builds we had during this run
              for (int t = 0; t <= run_data.steps; t++)
              {
                sys->set_step(time_step);
                sys->set_run_step(t);
                if (constraint->rescale())
                {
                  nlist->build();
                  nlist_builds++;
                }
                for (vector<DumpPtr>::iterator it_d = dump.begin(); it_d != dump.end(); it_d++)
                  (*it_d)->dump(time_step);
		for (vector<LoggerPtr>::iterator it_l = log.begin(); it_l != log.end(); it_l++)
                  (*it_l)->log();
                for (std::map<std::string, IntegratorPtr>::iterator it_integ = integrator.begin(); it_integ != integrator.end(); it_integ++)
                  (*it_integ).second->integrate();
                if (has_population)
                {
                  for (map<string,PopulationPtr>::iterator it_pop = population.begin(); it_pop != population.end(); it_pop++)
                  {
                    (*it_pop).second->divide(time_step);
                    if (sys->get_force_nlist_rebuild() && ((pot && pot->need_nlist()) || (aligner && aligner->need_nlist())))
                    {
                      nlist->build();
                      nlist_builds++;
                      sys->set_force_nlist_rebuild(false);
                    }
                    (*it_pop).second->remove(time_step);
                    if (sys->get_force_nlist_rebuild() && ((pot && pot->need_nlist()) || (aligner && aligner->need_nlist())))
                    {
                      nlist->build();
                      nlist_builds++;
                      sys->set_force_nlist_rebuild(false);
                    }
                    (*it_pop).second->grow(time_step);
                    if (sys->get_force_nlist_rebuild() && ((pot && pot->need_nlist()) || (aligner && aligner->need_nlist())))
                    {
                      if (sys->get_nlist_rescale() != 1.0)
                        nlist->rescale_cutoff(sys->get_nlist_rescale());
                      else
                        nlist->build();
                      nlist_builds++;
                      sys->set_force_nlist_rebuild(false);
                      sys->set_nlist_rescale(1.0);
                    }
                    (*it_pop).second->elongate(time_step);
                    if (sys->get_force_nlist_rebuild() && ((pot && pot->need_nlist()) || (aligner && aligner->need_nlist())))
                    {
                      if (sys->get_nlist_rescale() != 1.0)
                        nlist->rescale_cutoff(sys->get_nlist_rescale());
                      else
                        nlist->build();
                      nlist_builds++;
                      sys->set_force_nlist_rebuild(false);
                      sys->set_nlist_rescale(1.0);
                    }
                  }
                }
                // Check the neighbour list rebuild only if necessary 
                if ((pot && pot->need_nlist()) || (aligner && aligner->need_nlist()))
                {
                  bool nlist_rebuild = false;
                  if (sys->get_force_nlist_rebuild())
                  {
                    nlist_rebuild = true;
                    sys->set_force_nlist_rebuild(false);
                  }
                  else
                  {
                    for (int i = 0; i < sys->size(); i++)
                      if (nlist->need_update(sys->get_particle(i)))
                      {
                        nlist_rebuild = true;
                        break;
                      }
                  }
                  if (nlist_rebuild)
                  {
                    nlist->build();
                    nlist_builds++;
                  }
                }
                if (t % PRINT_EVERY == 0)
                  std::cout << "Time step: " << t <<"/" << run_data.steps << "   cumulative time step : " << time_step<< std::endl;
                time_step++;
              }
              msg->msg(Messenger::INFO,"Built neighbour list "+lexical_cast<string>(nlist_builds)+" time. Average number of steps between two builds : "+lexical_cast<string>(static_cast<double>(run_data.steps)/nlist_builds)+".");
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse number of run steps in line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing number of run steps.");
            }
          }
          else if (command_data.command == "pair_align")       // if command is pair align, parse it and add this pair alignment to the list of pairwise aligners
          {
            if (!defined["input"])  // We need to have system defined before we can add any aligners 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any pairwise aligners.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any pairwise aligners." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), align_parser, qi::space))
            {
              if (!has_aligner) 
              {
                aligner = std::make_shared<Aligner>(Aligner(sys,msg));
                has_aligner = true;
              }
              if (qi::phrase_parse(pair_align_data.params.begin(), pair_align_data.params.end(), param_parser, qi::space, parameter_data))
              {
                if (!defined["nlist"]) // Create a default neighbour list if non has been defined
                {
                  msg->msg(Messenger::WARNING,"No neighbour list defined. Some pair aligners (e.g. mf) need it. We are making one with default cutoff and padding.");
                  nlist = std::make_shared<NeighbourList>(NeighbourList(sys,msg,DEFAULT_CUTOFF,DEFAULT_PADDING, parameter_data));
                  defined["nlist"] = true;
                }
                aligner->add_pair_align(pair_align_data.type, pair_aligners[pair_align_data.type](sys,msg,nlist,parameter_data));
                msg->msg(Messenger::INFO,"Added "+pair_align_data.type+" to the list of pairwise aligners.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse pairwise alignment parameters for aligner type "+pair_align_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing pair aligner parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing pair_potential command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing pair_potential command.");
            }
          }
          else if (command_data.command == "align_param")       // parse pair parameters for potential
          {
            if (!has_aligner)  // We need to have pair aligners defined before we can change their parameters
            {
              msg->msg(Messenger::ERROR,"No pair aligners have been defined. Please define them using \"pair_align\" command before modifying any pair parameters.");
              throw std::runtime_error("No pair aligners defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), align_parser, qi::space))
            {
              if (qi::phrase_parse(pair_align_data.params.begin(), pair_align_data.params.end(), param_parser, qi::space, parameter_data))
              {
                aligner->add_pair_align_parameters(pair_align_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+pair_align_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse pair alignment parameters for potential type "+pair_align_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing pair alignment parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing align_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing align_param command.");
            }
            
          }
          else if (command_data.command == "external_align")       // if command is pair align, parse it and add this pair alignment to the list of pairwise aligners
          {
            if (!defined["input"])  // We need to have system defined before we can add any aligners 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any external aligners.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any external aligners." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), external_align_parser, qi::space))
            {
              if (!has_aligner) 
              {
                aligner = std::make_shared<Aligner>(Aligner(sys,msg));
                has_aligner = true;
                defined["external_aligner"] = true;
              }
              if (qi::phrase_parse(external_align_data.params.begin(), external_align_data.params.end(), param_parser, qi::space, parameter_data))
              {
                aligner->add_external_align(external_align_data.type, external_aligners[external_align_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Added "+external_align_data.type+" to the list of external aligners.");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse external alignment parameters for aligner type "+external_align_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing external aligner parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing external_align command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing external_align command.");
            }
          }
          else if (command_data.command == "external_align_param")       // parse parameters for external aligner
          {
            if (!defined["external_aligner"])  // We need to have external aligners defined before we can change their parameters
            {
              msg->msg(Messenger::ERROR,"No aligners have been defined. Please define them using \"external_align\" command before modifying any parameters.");
              throw std::runtime_error("No external aligners defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), external_align_parser, qi::space))
            {
              if (qi::phrase_parse(external_align_data.params.begin(), external_align_data.params.end(), param_parser, qi::space, parameter_data))
              {
                aligner->add_external_align_parameters(external_align_data.type, parameter_data);
                msg->msg(Messenger::INFO,"Setting new parameters for "+external_align_data.type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse external alignment parameters for aligner type "+external_align_data.type+" in line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing external alignment parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing external_align_param command at line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing external_align_param command.");
            }
            
          }
          else if (command_data.command == "group")       // if command is group, parse it and add this group to the list of all groups
          {
            if (!defined["input"])  // We need to have system defined before we can add any groups
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding any particle groups.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding any particle groups." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), group_parser, qi::space))
            {
              if (qi::phrase_parse(group_data.params.begin(), group_data.params.end(), param_parser, qi::space, parameter_data))
              {
                // Protect against groups boundary and tissue being redefined 
                if (group_data.name == "boundary" || group_data.name == "tissue")
                {
                  msg->msg(Messenger::ERROR,"Groups \"tissue\" and \"boundary\" are reserved and cannot be used.");
                  throw std::runtime_error("Trying to redefine reserved group name."); 
                }
                sys->make_group(group_data.name,parameter_data);
                msg->msg(Messenger::INFO,"Adding group "+group_data.name+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for group "+group_data.name+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing group parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing group command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing group line.");
            }
          }
          else if (command_data.command == "disable")       // if command is disable, parse it and disable given integrator 
          {
            if (!defined["integrator"])
            {
              msg->msg(Messenger::ERROR,"Integrator has not been defined. Please define it using \"integrator\" command before trying to disable it.");
              throw std::runtime_error("Integrator not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), disable_parser, qi::space))
            {
              if (qi::phrase_parse(disable_data.params.begin(), disable_data.params.end(), param_parser, qi::space, parameter_data))
              {
                std::string group_name;
                if (parameter_data.find("group") == parameter_data.end())
                  group_name = "all";
                else
                  group_name = parameter_data["group"];
                std::string to_disable = disable_data.type+"_"+group_name;
                if (integrator.find(to_disable) == integrator.end())
                {
                  msg->msg(Messenger::ERROR,"Integrator "+to_disable+" has not been defined. Cannot disable it.");
                  throw std::runtime_error("Trying to disable non-existent integrator.");
                }
                else
                {
                  msg->msg(Messenger::INFO,"Disabling integrator "+to_disable+".");
                  integrator.erase(to_disable);
                }
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for integrator disable "+disable_data.type+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing disable parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing disable command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing disable line.");
            }
          }
          else if (command_data.command == "population")       // if command is population, parse it and add create appropriate population object
          {
            if (!defined["input"])  // we need to have system defined before we can add any population 
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before adding population.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before adding population." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), population_parser, qi::space))
            {
              if (qi::phrase_parse(population_data.params.begin(), population_data.params.end(), param_parser, qi::space, parameter_data))
              {
                std::string group_name;
                if (parameter_data.find("group") == parameter_data.end())
                  group_name = "all";
                else
                  group_name = parameter_data["group"];
                if (population.find(population_data.type+"_"+group_name) == population.end())
                {
                  population[population_data.type+"_"+group_name] = std::shared_ptr<Population>(populations[population_data.type](sys,msg,parameter_data));  // dirty workaround shared_ptr and inherited classes
                  if (defined["nlist"])
                    population[population_data.type+"_"+group_name]->set_nlist(nlist);
                  msg->msg(Messenger::INFO,"Adding population of type "+population_data.type+" (for group "+group_name+").");
                  has_population = true;
                }
                else
                {
                  msg->msg(Messenger::ERROR,"Population control "+population_data.type+" (for type "+group_name+") already exists. Please disable it first.");
                  throw std::runtime_error("Population already exists.");
                }
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for population "+population_data.type+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing population parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse population command at line "+lexical_cast<string>(current_line));
              throw std::runtime_error("Error parsing population line.");
            }
          }
          else if (command_data.command == "disable_population")       // if command is disable_population, parse it and disable given population 
          {
            if (!defined["population"])
            {
              msg->msg(Messenger::ERROR,"Population has not been defined. Please define it using \"population\" command before trying to disable it.");
              throw std::runtime_error("Population not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), disable_population_parser, qi::space))
            {
              if (qi::phrase_parse(population_disable_data.params.begin(), population_disable_data.params.end(), param_parser, qi::space, parameter_data))
              {
                std::string group_name;
                if (parameter_data.find("group") == parameter_data.end())
                  group_name = "all";
                else
                  group_name = parameter_data["group"];
                std::string to_disable = population_disable_data.type+"_"+group_name;
                if (population.find(to_disable) == population.end())
                {
                  msg->msg(Messenger::ERROR,"Population "+to_disable+" has not been defined. Cannot disable it.");
                  throw std::runtime_error("Trying to disable non-existent population control.");
                }
                else
                {
                  msg->msg(Messenger::INFO,"Disabling population "+to_disable+".");
                  population.erase(to_disable);
                }
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for population disable "+disable_data.type+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing disable population parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing disable_population command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing disable_population line.");
            }
          }
          else if (command_data.command == "zero_momentum")       // if command is used to zero total momentum of a group of particles
          {
            if (!defined["input"])  // we need to have system defined before we can zero its momentum
            {
              if (defined["messages"])
                msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before trying to zero out momentum.");
              else
                std::cerr << "System has not been defined. Please define system using \"input\" command before trying to zero out momentum." << std::endl;
              throw std::runtime_error("System not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), param_parser, qi::space, parameter_data))
            {
              string particle_group = "all";
              if (parameter_data.find("group") != parameter_data.end()) particle_group = parameter_data["group"];
              sys->zero_cm_momentum(particle_group);
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse zero_momentum parameters at line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing zero_momentum parameters.");
            }
          }
          else if (command_data.command == "ntypes")       // run actual simulation
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), run_parser, qi::space))
            {
              msg->msg(Messenger::INFO,"Reading dummy keyword ntypes.");
              msg->write_config("ntypes",lexical_cast<string>(run_data.steps));
            }
            else
            {
              msg->msg(Messenger::ERROR,"Could not parse number of types for the dummy command ntypes in line : "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing number of types.");
            }
          }
          else if (command_data.command == "config")   // parse dump commands by adding list of logs
          {
            if (!defined["messages"])
            {
              std::cerr << "Messages file needs to be defined before we can write out configurations. Please use \"messages\" command before adding configs." << std::endl;
              throw std::runtime_error("Messages not defined.");
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), log_dump_parser, qi::space))
            {
              if (qi::phrase_parse(log_dump_data.params.begin(), log_dump_data.params.end(), param_parser, qi::space, parameter_data))
              {
                std::string config_file_type = "xml";
                if (parameter_data.find("type") != parameter_data.end())
                 config_file_type = parameter_data["type"];
                msg->add_config_file(log_dump_data.name,config_file_type);
                msg->msg(Messenger::INFO,"Adding config file "+log_dump_data.name+" of type "+config_file_type+".");
              }
              else
              {
                msg->msg(Messenger::ERROR,"Could not parse parameters for config file "+log_dump_data.name+" at line "+lexical_cast<string>(current_line)+".");
                throw std::runtime_error("Error parsing confing parameters.");
              }
            }
            else
            {
              msg->msg(Messenger::ERROR,"Error parsing config command as line "+lexical_cast<string>(current_line)+".");
              throw std::runtime_error("Error parsing config line.");
            }
          }
          else if (command_data.command == "timestep")       // if command is timestep, handle the global integrator timestep 
          {
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), timestep_parser, qi::space))  
            {
              if (!defined["input"])  // we need to have system defined before we can add any constraints
              {
                if (defined["messages"])
                  msg->msg(Messenger::ERROR,"System has not been defined. Please define system using \"input\" command before setting the integrator time step.");
                else
                  std::cerr << "System has not been defined. Please define system using \"input\" command before setting integrator time step." << std::endl;
                throw std::runtime_error("System not defined.");
              }
              else
              {
                if (defined["messages"])
                  msg->msg(Messenger::INFO,"Setting global integrator timestep to "+lexical_cast<string>(timestep_data.dt)+".");
                sys->set_integrator_step(timestep_data.dt);
              }
            }
            else
            {
              std::cerr << "Could not parse timestep command." << std::endl;
              throw std::runtime_error("Could not parse timestep command.");
            }
          } 
        }
        else
        {
          std::cerr << "Error parsing line : " << current_line << std::endl;
          std::cerr << "Unknown command : " << command_line << std::endl;
          throw std::runtime_error("Error parsing command script.");
        }
        
      }
      current_line++;
    }
  }
  else
  {
    std::cerr << "Could not open file : " << argv[1] << " for reading." << std::endl;
  }
  
  // Record time at the end
  std::time_t end_time = std::time(NULL);
  int run_time = static_cast<int>(std::difftime(end_time,start_time));
  int hours = run_time/3600;
  int minutes = (run_time % 3600)/60;
  int seconds = (run_time % 3600)%60;
  
  msg->msg(Messenger::INFO,"Simulation took "+lexical_cast<string>(hours)+" hours, "+lexical_cast<string>(minutes)+" minutes and "+lexical_cast<string>(seconds)+" seconds ("+lexical_cast<string>(run_time)+" seconds).");
  msg->msg(Messenger::INFO,"Average "+lexical_cast<string>(static_cast<double>(time_step)/run_time)+" time steps per second.");
  
  command_file.close();
  
  
  return 0;
}

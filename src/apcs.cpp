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
 * \file apcs.cpp
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

#include "defaults.hpp"
#include "messenger.hpp"
#include "dump.hpp"
#include "logger.hpp"
#include "parse_command.hpp"
#include "parse_constraint.hpp"
#include "parse_external.hpp"
#include "parse_input.hpp"
#include "parse_run.hpp"
#include "parse_potential.hpp"
#include "parse_aux.hpp"
#include "parse_parameters.hpp"
#include "parse_rng_seed.hpp"
#include "parse_box.hpp"
#include "parse_integrator.hpp"
#include "parse_log_dump.hpp"
#include "parse_align.hpp"
#include "parse_external_align.hpp"
#include "parse_group.hpp"
#include "parse_disable.hpp"
#include "parse_population.hpp"
#include "constraint.hpp"
#include "constraint_sphere.hpp"
#include "constraint_plane.hpp"
#include "constraint_plane_walls.hpp"
#include "constraint_cylinder.hpp"
#include "rng.hpp"
#include "particle.hpp"
#include "box.hpp"
#include "system.hpp"
#include "neighbour_list.hpp"
#include "external_potential.hpp"
#include "external_gravity_potential.hpp"
#include "pair_potential.hpp"
#include "pair_coulomb_potential.hpp"
#include "pair_soft_potential.hpp"
#include "pair_lj_potential.hpp"
#include "pair_gaussian_potential.hpp"
#include "pair_morse_potential.hpp"
#include "potential.hpp"
#include "integrator_brownian.hpp"
#include "integrator_vicsek.hpp"
#include "integrator_nve.hpp"
#include "integrator_nematic.hpp"
#include "integrator.hpp"
#include "aligner.hpp"
#include "pair_aligner.hpp"
#include "pair_polar_aligner.hpp"
#include "pair_nematic_aligner.hpp"
#include "pair_vicsek_aligner.hpp"
#include "external_aligner.hpp"
#include "population.hpp"
#include "population_random.hpp"


typedef boost::function< ConstraintPtr(SystemPtr, MessengerPtr, pairs_type&) > constraint_factory;                                                             //!< Type that handles all constraints
typedef boost::function< PairPotentialPtr(SystemPtr, MessengerPtr, NeighbourListPtr, pairs_type&) > pair_potential_factory;                                    //!< Type that handles all pair potentials
typedef boost::function< ExternalPotentialPtr(SystemPtr, MessengerPtr, pairs_type&) > external_potential_factory;                                              //!< Type that handles all external potentials
typedef boost::function< IntegratorPtr(SystemPtr, MessengerPtr, PotentialPtr, AlignerPtr, NeighbourListPtr, ConstraintPtr, pairs_type&) > integrator_factory;  //!< Type that handles all integrators
typedef boost::function< PairAlignPtr(SystemPtr, MessengerPtr, NeighbourListPtr, pairs_type&) > pair_aligner_factory;                                          //!< Type that handles all pairwise alignment
typedef boost::function< ExternalAlignPtr(SystemPtr, MessengerPtr, pairs_type&) > external_aligner_factory;                                                    //!< Type that handles all external alignment
typedef boost::function< PopulationPtr(SystemPtr, MessengerPtr, pairs_type&) > population_factory;                                                             //!< Type that handles all populations

int main(int argc, char* argv[])
{
  // Parser data
  CommandData         command_data;
  BoxData             box_data;
  InputData           input_data;
  ExternalData        external_data;
  LogDumpData         log_dump_data;
  PotentialData       potential_data;
  IntegratorData      integrator_data;
  ConstraintlData     constraint_data;
  RunData             run_data;
  AlignData           pair_align_data;
  ExternalAlignData   external_align_data;
  GroupData           group_data;
  DisableData         disable_data;
  PopulationData      population_data;
  pairs_type          parameter_data;   // All parameters for different commands
  
  // Parser grammars
  command_grammar         command_parser(command_data);
  box_grammar             box_parser(box_data);
  input_grammar           input_parser(input_data);
  external_grammar        external_parser(external_data); 
  log_dump_grammar        log_dump_parser(log_dump_data);
  potential_grammar       potential_parser(potential_data);
  integrator_grammar      integrator_parser(integrator_data);  
  constraint_grammar      constraint_parser(constraint_data);
  run_grammar             run_parser(run_data);
  align_grammar           align_parser(pair_align_data);
  external_align_grammar  external_align_parser(external_align_data);
  group_grammar           group_parser(group_data);
  disable_grammar         disable_parser(disable_data);
  population_grammar      population_parser(population_data);
  key_value_sequence      param_parser;
  
  // Class factories 
  std::map<std::string, constraint_factory> constraints;
  std::map<std::string, pair_potential_factory> pair_potentials;
  std::map<std::string, external_potential_factory> external_potentials;
  std::map<std::string, integrator_factory> integrators;
  std::map<std::string, pair_aligner_factory> pair_aligners;
  std::map<std::string, external_aligner_factory> external_aligners;
  std::map<std::string, population_factory> populations;
  
  std::ifstream command_file;    // File with simulation parameters and controls
  std::string command_line;      // Line from the command file
  
  MessengerPtr msg;                                // Messenger object
  BoxPtr box;                                      // Handles simulation box
  SystemPtr sys;                                   // System object
  PotentialPtr pot;                                // Handles all potentials
  ConstraintPtr constraint;                        // Handles the constraint to the manifold
  NeighbourListPtr nlist;                          // Handles global neighbour list
  std::map<std::string,IntegratorPtr> integrator;  // Handles the integrator
  AlignerPtr    aligner;                           // Handles all aligners
  ExternalAlignPtr external_aligner;               // Handles all external aligners
  vector<DumpPtr> dump;                            // Handles all different dumps
  vector<LoggerPtr> log;                           // Handles all different logs
  PopulationPtr  population;                       // Handles population
  
  bool periodic = false;       // If true, use periodic boundary conditions
  bool has_potential = false;  // If false, potential handling object (Potential class) has not be initialized yet
  bool has_aligner = false;    // If false, pairwise aligners are not defined
  bool has_population = false; // If false, population is not defined
  int time_step = 0;           // Counts current time step
  
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
  
    
  // Register spherical constraint with the constraint class factory
  constraints["sphere"] = boost::factory<ConstraintSpherePtr>();
  // Register xy plane constraint with the constraint class factory
  constraints["plane"] = boost::factory<ConstraintPlanePtr>();
  // Register plane walls constraint with the constraint class factory
  constraints["walls"] = boost::factory<ConstraintPlaneWallsPtr>();
  // Register cylindrical constraint with the constraint class factory
  constraints["cylinder"] = boost::factory<ConstraintCylinderPtr>();
  
  // Register Lennard-Jones pair potential with the pair potentials class factory
  pair_potentials["lj"] = boost::factory<PairLJPotentialPtr>();
  // Register Coulomb pair potential with the pair potentials class factory
  pair_potentials["coulomb"] = boost::factory<PairCoulombPotentialPtr>();
  // Register soft pair potential with the pair potentials class factory
  pair_potentials["soft"] = boost::factory<PairSoftPotentialPtr>();
  // Register Gaussian pair potential with the pair potentials class factory
  pair_potentials["gaussian"] = boost::factory<PairGaussianPotentialPtr>();  
  // Register Morse pair potential with the pair potentials class factory
  pair_potentials["morse"] = boost::factory<PairMorsePotentialPtr>();  
  
  // Register gravity to the external potential class factory
  external_potentials["gravity"] = boost::factory<ExternalGravityPotentialPtr>();
  
  // Register Brownian dynamics integrator with the integrators class factory
  integrators["brownian"] = boost::factory<IntegratorBrownianPtr>();
  // Register Vicsek dynamics integrator with the integrators class factory
  integrators["vicsek"] = boost::factory<IntegratorVicsekPtr>();
  // Register NVE integrator with the integrators class factory
  integrators["nve"] = boost::factory<IntegratorNVEPtr>();
  // Register nematic integrator with the integrators class factory
  integrators["nematic"] = boost::factory<IntegratorNematicPtr>();
  
  // Register polar aligner with the pairwise aligner class factory
  pair_aligners["polar"] = boost::factory<PairPolarAlignPtr>();
  // Register nematic aligner with the pairwise aligner class factory
  pair_aligners["nematic"] = boost::factory<PairNematicAlignPtr>();
  // Register Vicsek aligner with the pairwise aligner class factory
  pair_aligners["vicsek"] = boost::factory<PairVicsekAlignPtr>();
  
  // Register random population control with the class factory
  populations["random"] = boost::factory<PopulationRandomPtr>();
  
  
  if (argc < 2)
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << "    apcs <file_name>" << std::endl;
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
              box = boost::make_shared<Box>(Box(lx,ly,lz));
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
              msg = boost::shared_ptr<Messenger>(new Messenger(input_data.name));
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
                msg = boost::shared_ptr<Messenger>(new Messenger(msg_name));
                msg->msg(Messenger::WARNING,"Messenger was not defined prior to the reading in data. If not redefined all messages will be sent to "+msg_name+".");
              }
              if (!defined["box"])   // Cannot run without a box
              {
                msg->msg(Messenger::ERROR,"Simulation box has not been defined. Please define it before reading in coordinates.");
                throw std::runtime_error("Simulation box not defined.");
              }
              sys = boost::make_shared<System>(System(input_data.name,msg,box));
              sys->set_periodic(periodic);
              msg->msg(Messenger::INFO,"Finished reading system coordinates from "+input_data.name+".");
            }
            else
            {
              std::cerr << "Could not parse messenger command." << std::endl;
              throw std::runtime_error("Could not parse messenger command.");
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
            if (!defined["nlist"]) // Create a default neighbour list if none has been defined
            {
              msg->msg(Messenger::WARNING,"No neighbour list defined. Some pair potentials (e.g. Lennard-Jones) need it. We are making one with default cutoff and padding.");
              nlist = boost::make_shared<NeighbourList>(NeighbourList(sys,msg,DEFAULT_CUTOFF,DEFAULT_PADDING));
              defined["nlist"] = true;
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), potential_parser, qi::space))
            {
              if (!has_potential) 
              {
                pot = boost::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(potential_data.params.begin(), potential_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_pair_potential(potential_data.type, pair_potentials[potential_data.type](sys,msg,nlist,parameter_data));
                msg->msg(Messenger::INFO,"Added "+potential_data.type+" to the list of pair potentials.");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for pair potential "+potential_data.type+" is set to "+(*it).second+".");
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
                pot = boost::make_shared<Potential>(Potential(sys,msg));
                has_potential = true;
              }
              if (qi::phrase_parse(external_data.params.begin(), external_data.params.end(), param_parser, qi::space, parameter_data))
              {
                pot->add_external_potential(external_data.type, external_potentials[external_data.type](sys,msg,parameter_data));
                msg->msg(Messenger::INFO,"Added "+external_data.type+" to the list of external potentials.");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for pair potential "+external_data.type+" is set to "+(*it).second+".");
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
              if (qi::phrase_parse(constraint_data.params.begin(), constraint_data.params.end(), param_parser, qi::space, parameter_data))
              {
                constraint = boost::shared_ptr<Constraint>(constraints[constraint_data.type](sys,msg,parameter_data));  // dirty workaround shared_ptr and inherited classes
                msg->msg(Messenger::INFO,"Adding constraint of type "+constraint_data.type+".");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for constraint "+constraint_data.type+" is set to "+(*it).second+".");
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
              nlist = boost::make_shared<NeighbourList>(NeighbourList(sys,msg,rcut,pad));
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
                dump.push_back(boost::shared_ptr<Dump>(new Dump(sys,msg,log_dump_data.name,parameter_data)));
                msg->msg(Messenger::INFO,"Adding dump to file "+log_dump_data.name+".");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for dump "+log_dump_data.name+" is set to "+(*it).second+".");
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
                log.push_back(boost::shared_ptr<Logger>(new Logger(sys,msg,pot,aligner,log_dump_data.name,parameter_data)));
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
                if (integrator.find(integrator_data.type) == integrator.end())
                {
                  std::string group_name;
                  if (parameter_data.find("group") == parameter_data.end())
                    group_name = "all";
                  else
                    group_name = parameter_data["group"];
                  integrator[integrator_data.type+"_"+group_name] = boost::shared_ptr<Integrator>(integrators[integrator_data.type](sys,msg,pot,aligner,nlist,constraint, parameter_data));
                  msg->msg(Messenger::INFO,"Adding integrator of type "+integrator_data.type+".");
                  for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                    msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for integrator "+integrator_data.type+" is set to "+(*it).second+".");
                }
                else
                {
                  msg->msg(Messenger::ERROR,"Integrator of type "+integrator_data.type+" already exists. You need to disable the old one first.");
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
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for pair potential "+potential_data.type+" is set to "+(*it).second+".");
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
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for external potential "+external_data.type+" is set to "+(*it).second+".");
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
              // Precompute forces and torques
              if (pot)
                pot->compute();
              if (aligner)
                aligner->compute();
              int nlist_builds = 0;     // Count how many neighbour list builds we had during this run
              for (int t = 0; t <= run_data.steps; t++)
              {
                sys->set_step(time_step);
                for (vector<DumpPtr>::iterator it_d = dump.begin(); it_d != dump.end(); it_d++)
                  (*it_d)->dump(time_step);
                for (vector<LoggerPtr>::iterator it_l = log.begin(); it_l != log.end(); it_l++)
                  (*it_l)->log();
                for (std::map<std::string, IntegratorPtr>::iterator it_integ = integrator.begin(); it_integ != integrator.end(); it_integ++)
                  (*it_integ).second->integrate();
                if (has_population)
                {
                  population->divide(time_step);
                  population->remove(time_step);
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
            if (!defined["nlist"]) // Create a default neighbour list if non has been defined
            {
              msg->msg(Messenger::WARNING,"No neighbour list defined. Some pair aligners (e.g. mf) need it. We are making one with default cutoff and padding.");
              nlist = boost::make_shared<NeighbourList>(NeighbourList(sys,msg,DEFAULT_CUTOFF,DEFAULT_PADDING));
              defined["nlist"] = true;
            }
            if (qi::phrase_parse(command_data.attrib_param_complex.begin(), command_data.attrib_param_complex.end(), align_parser, qi::space))
            {
              if (!has_aligner) 
              {
                aligner = boost::make_shared<Aligner>(Aligner(sys,msg));
                has_aligner = true;
              }
              if (qi::phrase_parse(pair_align_data.params.begin(), pair_align_data.params.end(), param_parser, qi::space, parameter_data))
              {
                aligner->add_pair_align(pair_align_data.type, pair_aligners[pair_align_data.type](sys,msg,nlist,parameter_data));
                msg->msg(Messenger::INFO,"Added "+pair_align_data.type+" to the list of pairwise aligners.");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for pairwise alignment "+pair_align_data.type+" is set to "+(*it).second+".");
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
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for pair aligner "+pair_align_data.type+" is set to "+(*it).second+".");
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
                population = boost::shared_ptr<Population>(populations[population_data.type](sys,msg,parameter_data));  // dirty workaround shared_ptr and inherited classes
                msg->msg(Messenger::INFO,"Adding population of type "+population_data.type+".");
                for(pairs_type::iterator it = parameter_data.begin(); it != parameter_data.end(); it++)
                  msg->msg(Messenger::INFO,"Parameter " + (*it).first + " for population "+population_data.type+" is set to "+(*it).second+".");
                has_population = true;
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

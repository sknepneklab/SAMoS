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

#include <boost/functional/factory.hpp>
#include <boost/function.hpp>


#include "messenger.hpp"
#include "dump.hpp"
#include "parse_command.hpp"
#include "parse_constraint.hpp"
#include "parse_extrenal.hpp"
#include "parse_input.hpp"
#include "parse_run.hpp"
#include "parse_potential.hpp"
#include "parse_aux.hpp"
#include "parse_parameters.hpp"
#include "parse_rng_seed.hpp"
#include "parse_box.hpp"
#include "parse_integrator.hpp"
#include "parse_log_dump.hpp"
#include "constraint.hpp"
#include "constraint_sphere.hpp"
#include "constraint_plane.hpp"
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
#include "potential.hpp"
#include "integrator_brownian.hpp"
#include "integrator.hpp"

typedef boost::function< ConstraintPtr(SystemPtr, MessengerPtr, pairs_type&) > constraint_factory;                                                 //!< Type that handles all constraints
typedef boost::function< PairPotentialPtr(SystemPtr, MessengerPtr, NeighbourListPtr, pairs_type&) > pair_potential_factory;                        //!< Type that handles all pair potentials
typedef boost::function< ExternalPotentialPtr(SystemPtr, MessengerPtr, pairs_type&) > external_potential_factory;                                  //!< Type that handles all external potentials
typedef boost::function< IntegratorPtr(SystemPtr, MessengerPtr, PotentialPtr, NeighbourListPtr, ConstraintPtr, pairs_type&) > integrator_factory;  //!< Type that handles all integrators


int main(int argc, char* argv[])
{
  
  std::ifstream command_file;    // File with simulation parameters and controls
  std::string command_line;      // Line from the command file
  
  MessengerPtr msg;         // Messenger object
  BoxPtr box;               // Handles simulation box
  SystemPtr sys;            // System object
  PotentialPtr pot;         // Handles all potentials
  ConstraintPtr constraint; // Handles the constraint to the manifold
  NeighbourListPtr nlist;   // Handles global neighbour list
  vector<DumpPtr> dump;     // Handles all different dumps
  
  // flags that ensure that system the simulation is not in a bad state
  bool messenger_defined = false;    // If false, system has no messenger (set it to default "messages.msg")
  bool box_defined = false;          // If false, no simulation box defined (cannot run simulation)
  bool system_defined = false;       // If false, no system object defines (cannot run the simulation)
  bool potential_defined = false;    // If false, no potentials have been specified (cannot run simulation)
  bool constraint_defined = false;   // If false, no constraints have been defined (a full 3d simulation with a rather slow neighbour list building algorithm)
  bool nlist_defined = false;        // If false, system has no neighbour list
  
  // Class factories 
  std::map<std::string, constraint_factory> constraints;
  std::map<std::string, pair_potential_factory> pair_potentials;
  std::map<std::string, external_potential_factory> external_potentials;
  std::map<std::string, integrator_factory> integrators;
  
  // Register spherical constraint with the constraint class factory
  constraints["sphere"] = boost::factory<ConstraintSpherePtr>();
  // Register xy plane constraint with the constraint class factory
  constraints["plane"] = boost::factory<ConstraintPlanePtr>();
  
  // Register Lennard-Jones pair potential with the pair potentials class factory
  pair_potentials["lj"] = boost::factory<PairLJPotentialPtr>();
  // Register Coulomb pair potential with the pair potentials class factory
  pair_potentials["coulomb"] = boost::factory<PairCoulombPotentialPtr>();
  // Register soft pair potential with the pair potentials class factory
  pair_potentials["soft"] = boost::factory<PairSoftPotentialPtr>();
  
  // Register gravity to the external potential class factory
  external_potentials["gravity"] = boost::factory<ExternalGravityPotentialPtr>();
  
  // Register Brownian dynamics integrator with the integrators class factory
  integrators["brownian"] = boost::factory<IntegratorBrownianPtr>();
  
  if (argc < 1)
  {
    std::cerr << "Usage: " << std::endl;
    std::cerr << "    apcs <file_name>" << std::endl;
    return -1;
  }
  
  command_file.open(argv[1],std::ifstream::in);
  if (command_file)
  {
    while (std::getline(command_file, command_line)
    {
      // ... todo
    }
  }
  else
  {
    std::cerr << "Could not open file : " << anrg[1] << " for reading." << std::endl;
  }
  
  command_file.close();
  
  
  return 0;
}

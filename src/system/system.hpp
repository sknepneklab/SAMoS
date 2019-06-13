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
 * \file system.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2013
 * \brief Declaration of System class.
 */ 

#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__

#include <vector>
#include <string>
#include <stdexcept>
#include <exception>
#include <fstream>
#include <cmath>
#include <map>
#include <list>
#include <algorithm>
#include <memory>

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/lexical_cast.hpp>

#include "particle.hpp"
#include "bond.hpp"
#include "angle.hpp"
#include "box.hpp"
#include "messenger.hpp"
#include "group.hpp"
#include "defaults.hpp"
#include "mesh.hpp"

#include "parse_parameters.hpp"

using std::vector;
using std::string;
using std::runtime_error;
using std::ifstream;
using std::exception;
using std::sqrt;
using std::map;
using std::list;
using std::find;
using std::min_element;
using std::shared_ptr;
using std::make_shared;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::split_regex;
using boost::regex;
using namespace boost::algorithm;

/*! This class handles collection of all particles, i.e. the entire system.
 */

class System
{
public:
  
  //! Construct the system 
  System(const string&, MessengerPtr, BoxPtr);
  
  ~System() { m_particles.clear(); m_bonds.clear(); m_angles.clear(); }
  
  //! Get system size
  int size() { return m_particles.size(); } //!< \return Number of particles in the system.
  
  //! Get number of bonds
  int num_bonds() { return m_bonds.size(); } //!< \return Number of bonds in the system.
  
  //! Get number of angles
  int num_angles() { return m_angles.size(); } //!< \return Number of angles in the system.
  
  //! Get particle 
  //! \param i index of the particle to return 
  Particle& get_particle(int i) { return m_particles[i]; }  
  
  //! Get a bond
  //! \param i index of the bond to return
  Bond& get_bond(int i) { return m_bonds[i]; }
  
  //! Get an angle
  //! \param i index of the angle to return 
  Angle& get_angle(int i) { return m_angles[i]; }
  
  
  //! Get simulation box
  BoxPtr get_box() { return m_box; }
  
  //! Get mesh 
  Mesh& get_mesh() { return m_mesh; }
  
  //! Get the value of periodic boundary conditions flag
  bool get_periodic() { return m_periodic; }
  
  //! Set current time step
  //! \param step current time step
  void set_step(int step) { m_time_step = step; }
  
  //! Set current run step
  //! \param step current run step
  void set_run_step(int step) { m_run_step = step; }
  
  //! Get current time step
  int get_step() { return m_time_step; }  //!< \return current time step
  
  //! Get current run step
  int get_run_step() { return m_run_step; }  //!< \return current run step
  
  //! Reset all forces to zero
  void reset_forces()
  {
    // Reset all forces (accelerations) to zero
    for (int i = 0; i < this->size(); i++)
    {
      Particle& p = m_particles[i];
      p.fx = 0.0; p.fy = 0.0; p.fz = 0.0;
      p.s_xx = 0.0; p.s_xy = 0.0; p.s_xz = 0.0;
      p.s_yx = 0.0; p.s_yy = 0.0; p.s_yz = 0.0;
      p.s_zx = 0.0; p.s_zy = 0.0; p.s_zz = 0.0;
    }
  }
  
  //! Reset all torques to zero
  void reset_torques()
  {
    for (int i = 0; i < this->size(); i++)
    {
      Particle& p = m_particles[i];
      p.tau_x = 0.0; p.tau_y = 0.0; p.tau_z = 0.0;
    }
  }
  
  
  //! Set the periodic boundary conditions flag
  //! \param periodic value of the periodic boundary conditions flag
  void set_periodic(bool periodic) 
  { 
    if (periodic)
      m_msg->msg(Messenger::INFO,"Periodic flag for the system is set to true.");
    else
      m_msg->msg(Messenger::INFO,"Periodic flag for the system is set to false.");
    m_periodic = periodic; 
  }
   
  //! Get particle group
  //! \param name Group name
  GroupPtr get_group(const string name)  
  { 
    if (!(this->has_group(name)))
    {
      m_msg->msg(Messenger::ERROR,"Group "+name+" does not exist.");
      throw runtime_error("Non-existent particle group.");
    }
    return m_group[name];
  }

  //! Check is neighbour list needs forced rebuild (like after adding or removing particles
  bool get_force_nlist_rebuild() { return m_force_nlist_rebuild; }
  
  //! Set the force_nlist_rebuild flag
  //! \param val new value of the flag
  void set_force_nlist_rebuild(bool val) { m_force_nlist_rebuild = val; }
  
  //! Generate a group of particles
  void make_group(const string, pairs_type&);
  
  //! Check if a group exists
  //! \param name Group name 
  bool has_group(const string name)  
  {  
    if (m_group.find(name) == m_group.end())
      return false;
    return true;
  }
  
  //! Get the number of particle types
  int get_ntypes() const { return m_n_types; }
  
  //! Get the number of bond types
  int get_n_bond_types() const { return m_n_bond_types; }
  
  //! Get the number of angle types
  int get_n_angle_types() const { return m_n_angle_types; }
  
  //! Add particle to the system
  void add_particle(Particle&);
  
  //! Remove particle from the system
  void remove_particle(int);
  
  //! Move particle from one group to the other
  void change_group(int, const string&, const string&);
    
  //! Enable per particle energy tracking
  void enable_per_particle_eng() { m_compute_per_particle_eng = true; }
  
  //! disable per particle energy tracking
  void disable_per_particle_eng() { m_compute_per_particle_eng = false; }
  
  //! Compute per particle energy 
  bool compute_per_particle_energy() { return m_compute_per_particle_eng; }
  
  //! Zero centre of mass momentum
  void zero_cm_momentum(const string&);
  
  //! Read in all bonds from the bonds file
  void read_bonds(const string&);
  
  //! Read in all angles from the angles file
  void read_angles(const string&);
  
  //! Read in all boundary neighbours in tissue simulations 
  void read_boundary_neighbours(const string&);
  
  //! Get list of exclusions for a given particle
  //! \param i particle index
  vector<int>& get_exclusions(int i) { return m_exclusions[i]; }
 
  //! Return true is there are exclusions in the system
  bool has_exclusions() { return m_has_exclusions; }
 
  //! Check if a particle j is in the exclusion list of particle i
  //! \param i particle whose exclusions to check
  //! \param j check if this particle is in the exclusion list
  bool in_exclusion(int i, int j)
  {
    if (find(m_exclusions[i].begin(), m_exclusions[i].end(), j) == m_exclusions[i].end())
      return false;
    return true;
  }
  
  //! Compute tangent in the direction of the neighbour with the smallest index
  void compute_tangent(int, double&, double&, double&);
  
  //! Computes total system area for cell systems
  double compute_area();
  
  //! Computes average perimeter of all cells
  double compute_average_perimeter();
  
  //! Apply period boundary conditions on a quantity 
  void apply_periodic(double&, double&, double&);
  
  //! Make sure that all group information on particles matches group information in the lists
  bool group_ok(const string&);
  
  //! enforce period boundary conditions
  void enforce_periodic(Particle&);
  
  //! Set the neighbour list cutoff rescale parameter
  //! \param scale scale parameter
  void set_nlist_rescale(double scale) { m_nlist_rescale = scale; }
  
  //! Get value of the n_list rescale parameter
  double get_nlist_rescale() { return m_nlist_rescale; }
  
  //! Update mesh information for tissue simulations
  void update_mesh();
  
  //! Set the value of the integrator time step
  //! \param dt step size
  void set_integrator_step(double dt)  { m_dt = dt; }
  
  //! Get the value of the integrator time step
  double get_integrator_step() { return m_dt; }
  
  //! Set the number of mesh iteration 
  //! \param iter number of iteration 
  void set_max_mesh_iterations(int iter) { m_max_mesh_iter = iter; }

  //! Set value of the boundary type of particles in cell simulations 
  //! \param type boundary particle type
  void set_boundary_type(int type) { m_boundary_type = type; }
  
  //! Get list of all boundary particles
  vector<int>& get_boundary() { return m_boundary; }
  
  //! Add particle to the list of boundary particles
  //! \param id particle id
  void add_boundary(int id) { m_boundary.push_back(id); }
  
  //! Return true if system has boundary neighbours
  bool has_boundary_neighbours() { return m_has_boundary_neighbours; }
 
  //! Return list (vector) of all particles in a given molecule 
  //! \param mol_id id of the molecule
  vector<int>& get_mol_particles(int mol_id) { return m_molecules[mol_id]; }

  //! Return total number of molecules in the system
  int number_of_molecules() { return m_molecules.size(); }

  // Compute centre of mass of a molecule
  void molecule_cm(int, double&, double&, double&);

  // Set the record_force_type flag
  //! \param flag new value of the record_force_type flag
  void set_record_force_type(bool flag) { m_record_force_type = flag; }

  //! Get value of the record_force_type flag
  bool record_force_type() { return m_record_force_type; }
    
private:
  
  vector<Particle> m_particles;         //!< Contains all particle in the system 
  vector<Bond> m_bonds;                 //!< Contains all bonds in the system
  vector<Angle> m_angles;                //!< Contains all angles in the system
  MessengerPtr m_msg;                   //!< Handles messages sent to output
  BoxPtr m_box;                         //!< Simulation box object
  map<string, GroupPtr> m_group;        //!< All groups in the system 
  Mesh m_mesh;                          //!< Mesh is the data structure that holds tessallation infomration for tissue simulations
  bool m_periodic;                      //!< If true, we use periodic boundary conditions 
  int m_time_step;                      //!< Current time step
  int m_run_step;                       //!< Time step for the current run
  bool m_compute_per_particle_eng;      //!< If true, compute per particle potential and alignment energy (we need to be able to turn it on and off since it is slow - STL map in the inner loop!)
  int m_num_groups;                     //!< Total number of groups in the system
  bool m_force_nlist_rebuild;           //!< Forced rebuilding of neighbour list
  double m_nlist_rescale;               //!< Rescale neighbour list cutoff by this much
  int m_n_types;                        //!< Number of different particle types (used to set pair parameters) 
  int m_n_bond_types;                   //!< Number of different bond types
  int m_n_angle_types;                  //!< Number of different angle types
  int m_current_particle_flag;          //!< Keeps track of the last particle flag (distinct immutable id) of all particles. For bookkeeping. Clumsy as hell!  
  double m_dt;                          //!< This is a global integrator step used by all integrators (\note: it can be overwritten by a specific integrator)
  bool m_has_exclusions;                //!< If true, there are bonded interactions in the system and therefore those are accompanied with exclusions
  int m_max_mesh_iter;                  //!< Maximum number of iterations when cleaning up boundaries in the tissue simulations
  vector<vector<int> > m_exclusions;    //!< Which particles to be excluded from computing non bonded interactions (basically everything in bonds and angles)
  vector<int> m_boundary;               //!< Contains all particles that belong to the boundary (for tissue simulations)
  int m_boundary_type;                  //!< Type of the boundary particles that are added for cell simulations
  bool m_has_boundary_neighbours;       //!< If true, systems contains boundary neighbours (used in cells simulations)
  vector<vector<int> > m_molecules;     //!< List of all particles ids in a given molecule
  bool m_record_force_type;             //!< If true, for each particle record each force type that acts on it
   
};

typedef shared_ptr<System> SystemPtr;

#endif

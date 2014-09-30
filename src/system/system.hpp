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

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include "particle.hpp"
#include "box.hpp"
#include "messenger.hpp"
#include "group.hpp"

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
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::split_regex;
using boost::regex;
using boost::shared_ptr;
using boost::make_shared;
using namespace boost::algorithm;

/*! This class handles collection of all particles, i.e. the entire system.
 */

class System
{
public:
  
  //! Construct the system 
  System(const string&, MessengerPtr, BoxPtr);
  
  ~System() { m_particles.clear(); }
  
  //! Get system size
  int size() { return m_particles.size(); } //!< \return Number of particles in the system.
  
  //! Get particle 
  //! \param i index of the particle to return 
  Particle& get_particle(int i) { return m_particles[i]; }  
  
  //! Get simulation box
  BoxPtr get_box() { return m_box; }
  
  //! Get the value of periodic boundary conditions flag
  bool get_periodic() { return m_periodic; }
  
  //! Set current time step
  //! \param step current time step
  void set_step(int step) { m_time_step = step; }
  
  //! Get current time step
  int get_step() { return m_time_step; }  //!< \return current time step
  
  //! Reset all forces to zero
  void reset_forces()
  {
    // Reset all forces (accelerations) to zero
    for (int i = 0; i < this->size(); i++)
    {
      Particle& p = m_particles[i];
      p.fx = 0.0; p.fy = 0.0; p.fz = 0.0;
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
  
  //! Add particle to the system
  void add_particle(Particle&);
  
  //! Remove particle from the system
  void remove_particle(Particle&);
  
  //! Enable per particle energy tracking
  void enable_per_particle_eng() { m_compute_per_particle_eng = true; }
  
  //! disable per particle energy tracking
  void disable_per_particle_eng() { m_compute_per_particle_eng = false; }
  
  //! Compute per particle energy 
  bool compute_per_particle_energy() { return m_compute_per_particle_eng; }
  
private:
  
  vector<Particle> m_particles;         //!< Contains all particle in the system 
  MessengerPtr m_msg;                   //!< Handles messages sent to output
  BoxPtr m_box;                         //!< Simulation box object
  map<string, GroupPtr> m_group;        //!< All groups in the system 
  bool m_periodic;                      //!< If true, we use periodic boundary conditions 
  int m_time_step;                      //!< Current time step
  bool m_compute_per_particle_eng;      //!< If true, compute per particle potential and alignment energy (we need to be able to turn it on and off since it is slow - STL map in the inner loop!)
  int m_num_groups;                     //!< Total number of groups in the system
  bool m_force_nlist_rebuild;           //!< Forced rebuilding of neighbour list
  int m_n_types;                        //!< Number of different particle types (used to set pair parameters) 
    
};

typedef shared_ptr<System> SystemPtr;

#endif
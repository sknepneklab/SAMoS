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

#include <boost/regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>

#include "particle.hpp"
#include "box.hpp"
#include "messenger.hpp"

using std::vector;
using std::string;
using std::runtime_error;
using std::ifstream;
using std::exception;
using boost::lexical_cast;
using boost::bad_lexical_cast;
using boost::split_regex;
using boost::regex;
using boost::shared_ptr;
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
  
  //! Reset all forces to zero
  void reset_forces()
  {
    // Reset all forces (accelerations) to zero
    for (int i = 0; i < this->size(); i++)
    {
      Particle& p = m_particles[i];
      p.ax = 0.0; p.ay = 0.0; p.az = 0.0;
    }
  }
  
  //! Set the periodic boundary conditions flag
  //! \param periodic value of the periodic boundary conditions flag
  void set_periodic(bool periodic) { m_periodic = periodic; }
  
private:
  
  vector<Particle> m_particles;         //!< Contains all particle in the system 
  MessengerPtr m_msg;                   //!< Handles messages sent to output
  BoxPtr m_box;                         //!< Simulation box object
  bool m_periodic;                      //!< If true, we use periodic boundary conditions 
  
};

typedef shared_ptr<System> SystemPtr;

#endif
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
 * \file particle.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2013
 * \brief Declaration of Particle class.
 */ 

#ifndef __PARTICLE_HPP__
#define __PARTICLE_HPP__

#include <iostream>
#include <list>
#include <vector>
#include <map>
#include <string>

#include <boost/format.hpp>

using std::ostream;
using boost::format;
using std::map;
using std::string;
using std::list;
using std::vector;

const int NUM_PART_ATTRIB = 10;  //!< Number of particle attributes

/*! Particle class
 *  This is a basic unit in the simulation.
 *  It contains relevant parameters that describe each active particle,
 *  such as its position, velocity, type, etc.
 *  \note particle coordinates are given in the external, 3d Euclidean space.
 *  Its confinement to a given manifold (curved space) is enforced 
 *  via constraints.
 */
class Particle
{
public:
  //! Construct a Particle object
  //! \param id particle id
  //! \param type particle type
  //! \param r particle radius
  Particle(int id, int type, double r) : m_id(id), m_type(type), m_r(r) 
  { 
    fx = 0.0; fy = 0.0; fz = 0.0; 
    tau_x = 0.0; tau_y = 0.0; tau_z = 0.0;
    age = 0.0;
    m_flag = 0;
  }
  
  //! Get particle id
  int get_id() const { return m_id; } //!< \return particle id (m_id)
  
  //! Get particle type
  int get_type() const { return m_type; } //!< \return particle type (m_type)
  
  //! Get particle radius (size)
  double get_radius() const { return m_r; } //!< \return particle radius (m_r)
  
  //! Set id (for particle relabelling purposes)
  //! \param id new id
  void set_id(int id) { m_id = id; }
  
  //! Set type (for group changes and birth/death)
  //! \param type new type
  void set_type(int type) { m_type = type; }
  
  //! Set radius (for group changes and birth/death)
  //! \param r new type
  void set_radius(double r) { m_r = r; }
  
  
  //! Return potential energy of a given type
  //! \param type potential energy type to return 
  double get_pot_energy(const string& type)
  {
    if (m_pot_eng.find(type) != m_pot_eng.end())
      return m_pot_eng[type];
    else
      return 0.0;      
  }
  
  //! Set value of the potential energy of a given type
  //! \param type potential energy type (such as "soft")
  //! \param val value
  void set_pot_energy(const string& type, double val)
  {
    m_pot_eng[type] = val;
  }
  
  //! Add value of the potential energy of a given type
  //! \param type potential energy type (such as "soft")
  //! \param val value
  void add_pot_energy(const string& type, double val)
  {
    m_pot_eng[type] += val;
  }
  
  //! Return alignment potential energy of a given type
  //! \param type alignment potential energy type to return 
  double get_align_energy(const string& type)
  {
    if (m_align_eng.find(type) != m_align_eng.end())
      return m_align_eng[type];
    else
      return 0.0;      
  }
  
  //! Set value of the alignment potential energy of a given type
  //! \param type alignment potential energy type (such as "polar")
  //! \param val value
  void set_align_energy(const string& type, double val)
  {
    m_align_eng[type] = val;
  }
  
  //! Add value of the alignment potential energy of a given type
  //! \param type alignment potential energy type (such as "polar")
  //! \param val value
  void add_align_energy(const string& type, double val)
  {
    m_align_eng[type] += val;
  }
  
  //! Set flag parameter
  //! \param flag value of the flag 
  void set_flag(int flag) { m_flag = flag; }
  
  //! Get flag parameter
  int get_flag() const { return m_flag; }
  
  ///@{
  double x, y, z;              //!< Position in the embedding 3d flat space
  //@}
  ///@{
  double vx, vy, vz;           //!< Components of the particle velocity
  //@}
  ///@{
  double fx, fy, fz;           //!< Components of the force acting on the particle 
  //@}
  ///@{
  double tau_x, tau_y, tau_z;  //!< Components of the vector that handles director updates (torque in a sense?)
  //@}
  ///@{
  double nx, ny, nz;           //!< Particle direction vector (not necessarily equal to velocity direction)
  //@}
  double omega;                //!< Magnitude of the angular velocity (in the direction of the normal to the surface)
  double age;                  //!< Particle age (used when deciding to remove and split the particle)
  list<string> groups;         //!< List of all groups particle belongs to.
  vector<int> bonds;           //!< List of all bonds this particle belongs to
  vector<int> angles;          //!< List of all angles this particle belongs to
  
private:  // Make these attributes immutable 
  
  int m_id;                //!< Unique id
  int m_flag;              //!< This is a unique id which does not change regardless of removal of old and addition of new particles (note: clumsy, needs redesign)
  int m_type;              //!< Particle type
  double m_r;              //!< Particle radius 
  map<string,double> m_pot_eng;   //!< Holds current value of the potential energy of a given type 
  map<string,double> m_align_eng; //!< Holds alignment potential energy of a given type
    
};

ostream& operator<<(ostream&, const Particle&);


#endif
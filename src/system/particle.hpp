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

#include <boost/format.hpp>

using std::ostream;
using boost::format;

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
  }
  
  //! Get particle id
  int get_id() const { return m_id; } //!< \return particle id (m_id)
  
  //! Get particle type
  int get_type() const { return m_type; } //!< \return particle type (m_type)
  
  //! Get particle radius (size)
  double get_radius() const { return m_r; } //!< \return particle radius (m_r)
  
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
  
private:  // Make these attributes immutable 
  
  int m_id;                //!< Unique id
  int m_type;              //!< Particle type
  double m_r;              //!< Particle radius 
    
};

ostream& operator<<(ostream&, const Particle&);


#endif
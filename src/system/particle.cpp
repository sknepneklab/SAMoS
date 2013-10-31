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
 * \file particle.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2013
 * \brief Auxiliary functions for Particle class
 */ 

#include "particle.hpp"

ostream& operator<<(ostream& out, const Particle& p)
{
  out << format("id : %d\n") % p.get_id()            
  << format("type : %d\n") % p.get_type()
  << format("radius : %15.9e\n") % p.get_radius()            
  << format("(x,y,z) : (%15.9e,%15.9e,%15.9e)\n") % p.x % p.y % p.z
  << format("(vx,vy,vz) : (%15.9e,%15.9e,%15.9e)\n") % p.vx % p.vy % p.vz
  << format("(fx,fy,fz) : (%15.9e,%15.9e,%15.9e)\n") % p.fx % p.fy % p.fz
  << format("omega : %15.9e\n") % p.omega;            
  return out;
}
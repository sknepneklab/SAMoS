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
 * \file particle.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2013
 * \brief Auxiliary functions for Particle class
 */ 

#include "particle.hpp"

ostream& operator<<(ostream& out, const Particle& p)
{
  out << "-------- PARTICLE ------------" << endl;
  out << format("id : %d\n") % p.get_id()     
  << format("flag : %d\n") % p.get_flag()
  << format("type : %d\n") % p.get_type()
  << format("radius : %15.9e\n") % p.get_radius()       
  << format("age : %15.9e\n") % p.age    
  << format("parent : %d\n") % p.get_parent()
  << format("(x,y,z) : (%15.9e,%15.9e,%15.9e)\n") % p.x % p.y % p.z
  << format("(vx,vy,vz) : (%15.9e,%15.9e,%15.9e)\n") % p.vx % p.vy % p.vz
  << format("(nx,ny,nz) : (%15.9e,%15.9e,%15.9e)\n") % p.nx % p.ny % p.nz
  << format("(nvx,nvy,nvz) : (%15.9e,%15.9e,%15.9e)\n") % p.Nx % p.Ny % p.Nz
  << format("(fx,fy,fz) : (%15.9e,%15.9e,%15.9e)\n") % p.fx % p.fy % p.fz
  << format("omega : %15.9e\n") % p.omega
  << format("coordination : %d\n") % p.coordination
  << format("native cell area : %d\n") % p.A0
  << format("groups : ");
  for (list<string>::const_iterator it = p.groups.begin(); it != p.groups.end(); it++)
    out << format(" %s ") % (*it);
  if (p.bonds.size() > 0)
  {
    out << endl << format("bonds : ");
    for (vector<int>::const_iterator it_b = p.bonds.begin(); it_b != p.bonds.end(); it_b++)
      out << format(" %d ") % (*it_b);
  }
  if (p.angles.size() > 0)
  {
    out << endl << format("angles : ");
    for (vector<int>::const_iterator it_a = p.angles.begin(); it_a != p.angles.end(); it_a++)
      out << format(" %d ") % (*it_a);
  }
  out << endl;
  if (p.boundary)
    out << "particle belongs to the tissue boundary" << endl;
  out << "----------------------------";
  out << endl;
  return out;
}


    

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
 * \file constrainer.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Sep-2015
 * \brief Declaration of Constrainer class.
 */ 

#ifndef __CONSTRAINER_HPP__
#define __CONSTRAINER_HPP__

#include <map>
#include <string>

using std::map;
using std::string;

#include "messenger.hpp"
#include "system.hpp"

#include "constraint.hpp"

/*! Constainer class handles all constraints presnt in the system. All constraints are stored in a STL vector.
 * The integrator part of the code calls member functions of this class
 *  to apply actual constraints.
*/
class Constrainer
{
public:
  
  //! Construct Constrainer object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  Constrainer(SystemPtr sys, const MessengerPtr msg) : m_system(sys), m_msg(msg) { }
  
  //! Destructor
  ~Constrainer()
  {
    m_constraints.clear();
  }
  
  //! Add constraint to the list of all constraints
  //! \param name Unique name of the constraint
  //! \param constraint Pointer to the constraint
  void add_constraint(const string& name, ConstraintPtr constraint)
  {
    m_constraints.push_back(constraint);
    m_msg->msg(Messenger::INFO,"Added  : " + name + " to the list of constraints.");
  }
  
  //! Apply all constaints 
  void enforce(Particle& p)
  {
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      (*it_c)->enforce(p);
  }
  
  //! Rotate director around normal vector to the surface
  void rotate_director(Particle& p, double phi)
  {
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      (*it_c)->rotate_director(p,phi);
  }
  
  //! Rotate velocity around normal vector to the surface
  void rotate_velocity(Particle& p, double phi)
  {
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      (*it_c)->rotate_velocity(p,phi);
  }
  
  //! Project torque onto normal vector and return rotation angle change
  double project_torque(Particle& p)
  {
    double val = 0.0;
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      val += (*it_c)->project_torque(p);
    return val;
  }
  
  //! Computes normal to the surface
  //! \note That way thing are set up know only the normal to last constraint will be applied
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz)
  {
    Nx = 0.0; Ny = 0.0; Nz = 0.0;
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      if (find(p.groups.begin(),p.groups.end(),(*it_c)->get_group()) != p.groups.end())
      {
        double nx = 0.0, ny = 0.0, nz = 0.0;
        (*it_c)->compute_normal(p,nx,ny,nz);
        Nx += nx; Ny += ny;  Nz += nz;
      }
    double len_N = sqrt(Nx*Nx + Ny*Ny + Nz*Nz);
    Nx /= len_N;  Ny /= len_N;  Nz /= len_N;
  }
  
  bool rescale() 
  {
    bool res = false;
    for (vector<ConstraintPtr>::iterator it_c = m_constraints.begin(); it_c != m_constraints.end(); it_c++)
      res = res || (*it_c)->rescale();
    return res;
  }
  
private:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  
  vector<ConstraintPtr>  m_constraints;  //!< List of all constraints.
   
};

typedef shared_ptr<Constrainer> ConstrainerPtr;

#endif

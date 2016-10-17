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
 * \file constraint_tetrahedron.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Aug-2015
 * \brief Declaration of ConstraintTetrahedron class.
 */ 

#ifndef __CONSTRAINT_TETRAHERDON_HPP__
#define __CONSTRAINT_TETRAHERDON_HPP__

#include <cmath>

#include "system.hpp"
#include "vector3d.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sqrt;

const Vector3d tetra_a1 = Vector3d( 1.0, 0.0,-0.7071067811865475);
const Vector3d tetra_a2 = Vector3d(-1.0, 0.0,-0.7071067811865475);
const Vector3d tetra_a3 = Vector3d( 0.0, 1.0, 0.7071067811865475);
const Vector3d tetra_a4 = Vector3d( 0.0,-1.0, 0.7071067811865475);

/*! Enforces all particles to move on a surface with tetragonal symmetry.
 *  This surface is generated as the equipotential surface of 4 unit charges
 *  placed at the corners of a tetrahedron. The user can control their separation 
 *  and the level-set they wish to choose.
*/
class ConstraintTetrahedron : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., tetrahedron radius)
  ConstraintTetrahedron(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("scale") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Tetrahedral constraint. No scale set. Assuming 10.");
      m_s = 10.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Tetrahedral constraint. Scale set to "+param["scale"]+".");
      m_s = lexical_cast<double>(param["scale"]);
    }
    m_msg->write_config("constraint.tetrahedron.scale",lexical_cast<string>(m_scale));
    m_a1 = m_s*tetra_a1;
    m_a2 = m_s*tetra_a2;
    m_a3 = m_s*tetra_a3;
    m_a4 = m_s*tetra_a4;
    if (param.find("value") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Tetrahedral constraint. No level set value provided. Assuming 1.0");
      m_value = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Tetrahedral constraint. Level set value set at "+param["value"]+".");
      m_value = lexical_cast<double>(param["value"]);
    }
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle&, double&, double&, double&);
  
  // Computer gradient at a point
  void compute_gradient(Particle&, double&, double&, double&); 
  
  // Value of the constraint
  double constraint_value(Particle&);
  
  // Rescale constraint
  bool rescale();
      
private:
  
  double m_s;         //!< Tetrahedron size (essentially acts as the radius)
  double m_value;     //!< Level set (value of the potential)
  Vector3d m_a1;      //!< First vector of the four particles generating the surface
  Vector3d m_a2;      //!< Second vector of the four particles generating the surface
  Vector3d m_a3;      //!< Third vector of the four particles generating the surface
  Vector3d m_a4;      //!< Fourth vector of the four particles generating the surface
  
};

typedef shared_ptr<ConstraintTetrahedron> ConstraintTetrahedronPtr;  //!< Shared pointer to the Constraint object

#endif

/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file constraint_plane_walls.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 13-Aug-2014
 * \brief Declaration of ConstraintPlaneWalls class.
 */ 

#ifndef __CONSTRAINT_PLANE_WALLS_HPP__
#define __CONSTRAINT_PLANE_WALLS_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sin;
using std::cos;

/*! Enforces all particles to lay on the xy plane of size Lx x Ly
 *  Periodic boundary in the y-direction and repulsive walls at
 *  -L and L in the x-direction. z-coordinate set to zero.
*/
class ConstraintPlaneWalls : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintPlaneWalls(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("l") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Plan walls constraint. No wall position set. Assuming 1.0.");
      m_l = 1.0;
    }
    else
    {
      double l = lexical_cast<double>(param["l"]);
      if (l > m_system->get_box()->xhi)
      {
        m_msg->msg(Messenger::WARNING,"Plan walls constraint. Wall position beyond box limits. Setting to box limits.");
        m_l = m_system->get_box()->xhi;
      }
      else
        m_l = l;
      m_msg->msg(Messenger::INFO,"Plane walls constraint. Wall position set to "+param["l"]+".");
    }
    m_msg->write_config("constraint.plane_wall.l",lexical_cast<string>(m_l));
  }
  
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the plane (z axis)
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the plane (z axis)
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector to the plane (z axis) and return rotation angle change
  double project_torque(Particle&);
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz) { Nx = 0.0; Ny = 0.0; Nz = 1.0;  }
  
  // Computer gradient at a point
  void compute_gradient(Particle& p, double& gx, double& gy, double& gz) { }
  
  // Value of the constraint
  double constraint_value(Particle& p) { return 0.0; }
    
private:
  
  double m_l;     //!< position of the wall in x-direction (other wall will be at -m_l)
  
};

typedef shared_ptr<ConstraintPlaneWalls> ConstraintPlaneWallsPtr;  //!< Shared pointer to the ConstraintPlaneWalls object

#endif
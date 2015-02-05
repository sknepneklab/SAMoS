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
 * \file constraint_actomyo.hpp
 * \author Amit Das, dosamit@gmail.com
 * \date 2--Nov-2014
 * \brief Declaration of ConstraintActomyo class.
 */ 

#ifndef __CONSTRAINT_ACTOMYO_HPP__
#define __CONSTRAINT_ACTOMYO_HPP__

#include <cmath>

#include "system.hpp"
#include "parse_parameters.hpp"
#include "constraint.hpp"

using std::sin;
using std::cos;


/*! Enforces all particles of a actin to lay on a plane parallel to the xy direction of size Lx x Ly and
 *   particles of myosin to move at another plane of same size above the former containing the actins
 *  Also makes sure periodic boundary conditions are enforced
*/
class ConstraintActomyo : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintActomyo(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    m_lx = m_system->get_box()->Lx;
    m_ly = m_system->get_box()->Ly;
    m_lz = m_system->get_box()->Lz;
    m_msg->msg(Messenger::INFO,"Actomyosin constraint. Setting box dimensions to lx = "+lexical_cast<string>(m_lx)+", ly = "+lexical_cast<string>(m_ly)+", and lz ="+lexical_cast<string>(m_lz)+".");
  }
  
  //! Computes normal to the surface
  void compute_normal(Particle& p, double& Nx, double& Ny, double& Nz) 
  {
    Nx = 0.0; Ny = 0.0; Nz = 1.0; 
  }
  
  //! Enforce constraint
  void enforce(Particle&);
  
  //! Rotate director around normal vector to the plane (z axis)
  void rotate_director(Particle&, double);
  
  //! Rotate velocity around normal vector to the plane (z axis)
  void rotate_velocity(Particle&, double);
  
  //! Project torque onto normal vector to the plane (z axis) and return rotation angle change
  double project_torque(Particle&);
    
private:
  
  double m_lx;     //!< box size in x direction 
  double m_ly;     //!< box size in y direction 
  double m_lz;     //!< box size in z direction 
  
};

typedef shared_ptr<ConstraintActomyo> ConstraintActomyoPtr;  //!< Shared pointer to the Constraint object

#endif

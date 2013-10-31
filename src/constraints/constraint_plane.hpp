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
 * \file constraint_plane.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Declaration of ConstraintPlane class.
 */ 

#ifndef __CONSTRAINT_PLANE_HPP__
#define __CONSTRAINT_PLANE_HPP__

#include <cmath>

#include "system.hpp"

#include "parse_parameters.hpp"


/*! Enforces all particles to lay on the xy plane of size Lx x Ly
 *  Also makes sure periodic boundary conditions are enforced
*/
class ConstraintPlane : public Constraint
{
public:
  
  //! Constructor
  //! \param id unique constraint id
  //! \param sys pointer to the system object
  //! \param msg Pointer to the internal state messenger
  //! \param param parameters that define the manifolds (e.g., sphere radius)
  ConstraintPlane(SystemPtr sys, MessengerPtr msg, pairs_type& param) : Constraint(sys,msg,param)
  { 
    if (param.find("r") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Spherical constraint. No radius set. Assuming 1.");
      m_r = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Spherical constraint. Radius set to "+param["r"]+".");
      m_r = lexical_cast<double>(param["r"]);
    }
    m_lx = m_system->get_box()->Lx;
    m_ly = m_system->get_box()->Ly;
  }
  
  //! Enforce constraint
  void enforce();
    
private
  
  double m_lx;     //!< box size in x direction 
  double m_ly;     //!< box size in y direction 
  
};

typedef shared_ptr<ConstraintPlane> ConstraintPlanePtr;  //!< Shared pointer to the Constraint object

#endif
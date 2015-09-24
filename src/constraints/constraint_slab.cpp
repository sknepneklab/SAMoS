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
 * \file constraint_slab.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Sep-2015
 * \brief Implementation of the slab constraint
 */ 

#include "constraint_slab.hpp"

/*! Force all particles to be confined move between planes at \f$ z = z_{lo} \f$ and 
 *  \f$ z = z_{hi} \f$.
 */
void ConstraintSlab::enforce(Particle& p)
{
  bool apply = false;
  if (m_group == "all")
    apply = true;
  else
    apply = (find(p.groups.begin(),p.groups.end(),m_group) != p.groups.end());
  if (apply)
  {
    if (p.z < m_z_lo)
    {
      p.z = m_z_lo;
      p.vz = -p.vz;
      p.nz = -p.nz;
    }
    if (p.z > m_z_hi)
    {
      p.z = m_z_hi;
      p.vz = -p.vz;
      p.nz = -p.nz;
    }
  }
  if (m_system->get_periodic())
    m_system->enforce_periodic(p);
}




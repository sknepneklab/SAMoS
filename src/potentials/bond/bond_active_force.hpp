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
 * \file bond_active_force.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 12-Aug-2015
 * \brief Declaration of BondActiveForce class
 */ 

#ifndef __BOND_HARMONIC_FORCE_HPP__
#define __BOND_HARMONIC_FORCE_HPP__

#include <cmath>

#include "bond_potential.hpp"

using std::sqrt;

//! Structure that handles parameters for the harmonic bond
struct BondActiveParameters
{
  double f;   // Intensity of active force
};

/*! BondActiveForce implements actie force along the chain. Force is along the 
 *  tangent of the chain, which is in this case along the bond. It is
 *  \f$ \vec F^{(i)}_{active} = f\vec\tau_i \f$ where \f$ f \f$ is the strength of the force
 *  and \f$ \tau_i = (\vec r_{i+1} - \vec r_i)/|\vec r_{i+1} - \vec r_i| \f$ is the direction 
 *  force is pointing in (i.e. along the bond connecting beads \f$ \vec r_i \f$ and \f$ r_{i+j} \f$.
 *  \note Both beads will get the same force (this is extrnal active force!) and the force on each 
 *  bead will be sum of the forces acting on all bonds that are connected to it.
 *  \note This class assumes that beads in consecutive bonds are ordered, i.e., that the chain has
 *  clearly defined head and the tail.
 *  \note Potential would here make very little sense, so we do not define it.
 */
class BondActiveForce : public BondPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (k and l0)
  BondActiveForce(SystemPtr sys, MessengerPtr msg,  pairs_type& param) : BondPotential(sys, msg, param)
  {
    int ntypes = m_system->get_n_bond_types();
    if (param.find("f") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Bond Active force: No active force intensity given. Setting it to 1.");
      m_f = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Bond Active Force: Active force strength (f) set to "+param["f"]+".");
      m_f = lexical_cast<double>(param["f"]);
    }
    m_msg->write_config("potential.bond.active.f",lexical_cast<string>(m_f));
    m_bond_params = new BondActiveParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_bond_params[i].f = m_f;
    }
  }
  
  virtual ~BondActiveForce()
  {
    delete [] m_bond_params;
  }
                                                                                                                
  //! Set parameters data for specific bond type    
  void set_bond_parameters(pairs_type& bond_param)
  {
    map<string,double> param;
    
    int type;
    
    if (bond_param.find("type") == bond_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Bond type has not been defined for active bond force.");
      throw runtime_error("Missing key for active bond parameters.");
    }
    type = lexical_cast<int>(bond_param["type"]);
        
    if (bond_param.find("f") != bond_param.end())
    {
      m_msg->msg(Messenger::INFO,"Active bond force. Setting active force strength "+bond_param["f"]+" for bond of types "+lexical_cast<string>(type)+".");
      param["f"] = lexical_cast<double>(bond_param["f"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Active bond force. Using default active force strength ("+lexical_cast<string>(m_f)+") for bonds of types "+lexical_cast<string>(type)+".");
      param["f"] = m_f;
    }
    m_msg->write_config("potential.bond.active.type_"+bond_param["type"]+".f",lexical_cast<string>(param["f"]));
        
    m_bond_params[type-1].f = param["f"];
         
    m_has_bond_params = true;
  }
  
   
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_f;                       //!< active force strength
  BondActiveParameters* m_bond_params;   //!< type specific bond parameters 
    
};

typedef shared_ptr<BondActiveForce> BondActiveForcePtr;

#endif

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
 * \file bond_fene_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Sept-2015
 * \brief Declaration of BondFenePotential class
 */ 

#ifndef __BOND_FENE_POTENTIAL_HPP__
#define __BOND_FENE_POTENTIAL_HPP__

#include <cmath>

#include "bond_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the fene bond
struct BondFeneParameters
{
  double k;    // Spring constant
  double r0;   // maximum bond lenght
};

/*! BondFenePotential implements standard FEME bond potential 
 *  Potential is given as \f$ U_{fene}\left(r_{ij}\right) = -\frac{1}{2}kr_0^2\log\left[1-\left(\frac{\left|\vec r_{ij}\right|}{r_0}\right)^2\right] + WCA(r_{ij})\f$,
 *  where \f$ k \f$ is the spring constant, \f$ r_0 \f$ is the the maximum bond length and \f$ \vec r_{ij} = \vec r_i - \vec r_j \f$. WCA is the repulsive part
 *  of Lennard-Jones potential (so called, Weeks-Chandler-Andersen potential), which is *NOT* included here. Therfore, the user has to make sure to disable
 *  exclusions in bonds in order to properly treat the short-range part of the potential.
 */
class BondFenePotential : public BondPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (k and l0)
  BondFenePotential(SystemPtr sys, MessengerPtr msg,  pairs_type& param) : BondPotential(sys, msg, param)
  {
    int ntypes = m_system->get_n_bond_types();
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No spring constant (k) specified for fene bond potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global spring constant (k) for fene bond potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.bond.fene.k",lexical_cast<string>(m_k));
    if (param.find("r0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No maximum bond length (r0) specified for fene bond potential. Setting it to 1.5.");
      m_r0 = 1.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global maximum bond length (r0) for fene bond potential is set to "+param["r0"]+".");
      m_r0 = lexical_cast<double>(param["r0"]);
    }
    m_msg->write_config("potential.bond.fene.r0",lexical_cast<string>(m_r0));
    m_bond_params = new BondFeneParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_bond_params[i].k = m_k;
      m_bond_params[i].r0 = m_r0;
    }
  }
  
  virtual ~BondFenePotential()
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
      m_msg->msg(Messenger::ERROR,"Bond type has not been defined for fene bond potential.");
      throw runtime_error("Missing key for fene bond parameters.");
    }
    type = lexical_cast<int>(bond_param["type"]);
        
    if (bond_param.find("k") != bond_param.end())
    {
      m_msg->msg(Messenger::INFO,"Fene bond potential. Setting spring constant to "+bond_param["k"]+" for bond of types "+lexical_cast<string>(type)+".");
      param["k"] = lexical_cast<double>(bond_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Fene bond potential. Using default spring constant ("+lexical_cast<string>(m_k)+") for bonds of types "+lexical_cast<string>(type)+".");
      param["k"] = m_k;
    }
    if (bond_param.find("r0") != bond_param.end())
    {
      m_msg->msg(Messenger::INFO,"Fene bond potential. Setting maximum bond length to "+bond_param["r0"]+" for bond of types "+lexical_cast<string>(type)+".");
      param["r0"] = lexical_cast<double>(bond_param["r0"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Fene bond potential. Using default maximum bond length ("+lexical_cast<string>(m_r0)+") for bonds of types "+lexical_cast<string>(type)+".");
      param["r0"] = m_r0;
    }
    m_msg->write_config("potential.bond.fene.type_"+bond_param["type"]+".k",lexical_cast<string>(param["k"]));
    m_msg->write_config("potential.bond.fene.type_"+bond_param["type"]+".r0",lexical_cast<string>(param["r0"]));
    
    m_bond_params[type-1].k = param["k"];
    m_bond_params[type-1].r0 = param["r0"];
         
    m_has_bond_params = true;
  }
  
   
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_k;       //!< spring constant
  double m_r0;      //!< maximum bond length
  BondFeneParameters* m_bond_params;   //!< type specific bond parameters 
    
};

typedef shared_ptr<BondFenePotential> BondFenePotentialPtr;

#endif

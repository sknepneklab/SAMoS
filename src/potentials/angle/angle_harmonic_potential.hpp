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
 * \file angle_harmonic_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Declaration of AngleHarmonicPotential class
 */ 

#ifndef __ANGLE_HARMONIC_POTENTIAL_HPP__
#define __ANGLE_HARMONIC_POTENTIAL_HPP__

#include <cmath>

#include "angle_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::acos;

//! Structure that handles parameters for the harmonic angle
struct AngleHarmoicParameters
{
  double k;    // stiffness
  double t0;   // equilibrium angle
};

/*! AngleHarmonicPotential implements standard harmonic angle potential 
 *  Potential is given as \f$ U_{harm}\left(\theta\right) = k_\theta\left(\theta-\theta_0\right)^2 \f$,
 *  where \f$ k_\theta \f$ is the stiffness, \f$ \theta \f$ is the angle between two bonds originating from the particle
 *  and \f$ \theta_0 \f$ is the rest (equilibrium) value of the angle. Note that we absorbed usual prefactor of \f$ 1/2 \f$
 *  into the definition of \f$ k_\theta \f$.
 */
class AngleHarmonicPotential : public AnglePotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (k and t0)
  AngleHarmonicPotential(SystemPtr sys, MessengerPtr msg,  pairs_type& param) : AnglePotential(sys, msg, param)
  {
    int ntypes = m_system->get_n_angle_types();
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No stiffness (k) specified for harmonic angle potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global stiffness (k) for harmonic angle potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.angle.harmonic.k",lexical_cast<string>(m_k));
    if (param.find("t_eq") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No equilibrium angle (t0) specified for harmonic angle potential. Setting it to pi.");
      m_t0 = M_PI;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global equilibrium angle (t0) for harmonic angle potential is set to "+param["t_eq"]+".");
      m_t0 = lexical_cast<double>(param["t_eq"]);
    }
    m_msg->write_config("potential.angle.harmonic.t_eq",lexical_cast<string>(m_t0));
    m_angle_params = new AngleHarmoicParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_angle_params[i].k = m_k;
      m_angle_params[i].t0 = m_t0;
    }
  }
  
  virtual ~AngleHarmonicPotential()
  {
    delete [] m_angle_params;
  }
                                                                                                                
  //! Set parameters data for specific angle type    
  void set_angle_parameters(pairs_type& angle_param)
  {
    map<string,double> param;
    
    int type;
    
    if (angle_param.find("type") == angle_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Angle type has not been defined for harmonic angle potential.");
      throw runtime_error("Missing key for harmonic angle parameters.");
    }
    type = lexical_cast<int>(angle_param["type"]);
        
    if (angle_param.find("k") != angle_param.end())
    {
      m_msg->msg(Messenger::INFO,"Harmonic angle potential. Setting stiffness to "+angle_param["k"]+" for angle of types "+lexical_cast<string>(type)+".");
      param["k"] = lexical_cast<double>(angle_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Harmonic angle potential. Using default stiffness ("+lexical_cast<string>(m_k)+") for angles of types "+lexical_cast<string>(type)+".");
      param["k"] = m_k;
    }
    if (angle_param.find("t_eq") != angle_param.end())
    {
      m_msg->msg(Messenger::INFO,"Harmonic angle potential. Setting equilibrium angle to "+angle_param["t_eq"]+" for angle of types "+lexical_cast<string>(type)+".");
      param["t_eq"] = lexical_cast<double>(angle_param["t_eq"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Harmonic angle potential. Using default equilibrium angle ("+lexical_cast<string>(m_t0)+") for angles of types "+lexical_cast<string>(type)+".");
      param["t_eq"] = m_t0;
    }
    m_msg->write_config("potential.angle.harmonic.type_"+angle_param["type"]+".k",lexical_cast<string>(param["k"]));
    m_msg->write_config("potential.angle.harmonic.type_"+angle_param["type"]+".t_eq",lexical_cast<string>(param["t_eq"]));
    
    m_angle_params[type-1].k = param["k"];
    m_angle_params[type-1].t0 = param["t_eq"];
         
    m_has_angle_params = true;
  }
  
   
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_k;       //!< stiffness
  double m_t0;      //!< equilibrium angle
  AngleHarmoicParameters* m_angle_params;   //!< type specific angle parameters 
    
};

typedef shared_ptr<AngleHarmonicPotential> AngleHarmonicPotentialPtr;

#endif

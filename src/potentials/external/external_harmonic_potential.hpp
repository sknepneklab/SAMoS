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
 * \file external_harmonic_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 19-Aug-2015
 * \brief Declaration of ExternalHarmonicPotential class
 */ 

#ifndef __EXTERNAL_HARMONIC_POTENTIAL_HPP__
#define __EXTERNAL_HARMONIC_POTENTIAL_HPP__

#include "external_potential.hpp"


/*! ExternalHarmonicPotential handles external harmonic potential with repsect
 *  to a plane parallel to the xy-plane.
 *  In other words,
 *  \f$ V_{harm}\left(r\right) = \frac{1}{2} k (z-z0)^2 \f$ and the force is \f$ \vec F = - k (z - z0) \vec e_z \f$,
 *  where \f$ k \f$ is the stiffness constant and \f$ z_0 \f$ is the position of the plane along the z-axis.
*/
class ExternalHarmonicPotential : public ExternalPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (harmonic strength)
  ExternalHarmonicPotential(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalPotential(sys, msg, param)  
  {
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No stiffness specified for the external harmonic potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Stiffness of the external harmonic potential set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.external.harmonic.k",lexical_cast<string>(m_k));
    if (param.find("z0") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No position of the zero for the external harmonic potential given. Setting it to 0.");
      m_z0 = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Position of the zero for the external harmonic potential set to "+param["z0"]+".");
      m_z0 = lexical_cast<double>(param["z0"]);
    }
    m_msg->write_config("potential.external.harmonic.z0",lexical_cast<string>(m_z0));
  }
                                                                                                                
  //! Get the total potential energy
  double get_potential_energy() { return m_potential_energy; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  void set_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for external gravitational potential.");
      throw runtime_error("Missing key for external potential parameters.");
    }
    
    int type = lexical_cast<int>(pair_param["type"]);
        
    if (pair_param.find("k") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"External harmonic potential. Setting k to "+pair_param["k"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["k"] = lexical_cast<double>(pair_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External harmonic potential. Using default k ("+lexical_cast<string>(m_k)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["k"] = m_k;
    }
    m_msg->write_config("potential.external.harmonic.type_"+pair_param["type"]+".k",lexical_cast<string>(param["k"]));
        
    if (pair_param.find("z0") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"External harmonic potential. Setting z0 to "+pair_param["z0"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["z0"] = lexical_cast<double>(pair_param["z0"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External harmonic potential. Using default z0 ("+lexical_cast<string>(m_z0)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["z0"] = m_z0;
    }
    m_msg->write_config("potential.external.harmonic.type_"+pair_param["type"]+".z0",lexical_cast<string>(param["z0"]));
    
    
    m_type_params[type-1]["k"] = param["k"];
    m_type_params[type-1]["z0"] = param["z0"];
    
    m_has_params = true;
  }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_k;     //!< stiffness of external harmonic potential
  double m_z0;    //!< position of the zero of the potential 
  
};

typedef shared_ptr<ExternalHarmonicPotential> ExternalHarmonicPotentialPtr;

#endif

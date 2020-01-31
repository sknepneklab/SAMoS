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
 * \file external_boundary_pull.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Sep-2016
 * \brief Declaration of ExternalBoundaryPull class
 */ 

#ifndef __EXTERNAL_BOUNDARY_PULL_HPP__
#define __EXTERNAL_BOUNDARY_PULL_HPP__

#include "external_potential.hpp"


/*! ExternalBoundaryPull handles pulling force onto the boundary of a tissue layer. 
 *  For a boundary particle \f$ i \f$ we apply a force of magnitude \f$ \alpha \f$ in 
 *  the direction \f$ -\left(\vec r{ji} + \vec r_{ki}\right) \f$, where \f$ \vec r_{ji} \f$ 
 *  in the vecror along the edge connecting particle \f$ i \f$ and it boundary neighbour
 *  \f$ j \f$. Similarly, \f$ \vec r_{ki} $\f is the vector connecting particle \f$ i \f$
 *  with its other boundary neighbour \f$ k \f$. Both vectors point away from particle \f$ i \f$. 
 *  \note This is not a potential, but to avoid code bloating we place it
 *  in the same class structure as regular external potential.
 *  \note The force always acts in the direction \f$ -\left(\vec r{ji} + \vec r_{ki}\right) \f$. For 
 *  a non-convex region of the boundary can lead to the inward pointing force. 
*/
class ExternalBoundaryPull : public ExternalPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (gravity strength)
  ExternalBoundaryPull(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalPotential(sys, msg, param)  
  {
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No force magnitude specified for boundary pull. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"force magnitude in boundary pull set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.external.boundary_pull.alpha",lexical_cast<string>(m_alpha));
  }
                                                                                                                
  //! Get the total potential energy (this is not a potential so remove 0!)
  double get_potential_energy() { return 0.0; } //!< \return value of the total potential energy
  
  //! Set pair parameters data for pairwise interactions    
  void set_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Particle type has not been defined for boundary pull.");
      throw runtime_error("Missing key for boundary pull.");
    }
    
    int type = lexical_cast<int>(pair_param["type"]);
        
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"boundary pull. Setting alpha to "+pair_param["alpha"]+" for particle pair of type "+lexical_cast<string>(type)+".");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"boundary pull. Using default alpha ("+lexical_cast<string>(m_alpha)+") for particle pair of type "+lexical_cast<string>(type)+").");
      param["alpha"] = m_alpha;
    }
    m_msg->write_config("potential.external.boundary_pull.type_"+pair_param["type"]+".alpha",lexical_cast<string>(param["alpha"]));

    m_type_params[type-1]["alpha"] = param["alpha"];
    
    m_has_params = true;
  }
  
  //! Computes forces for all particles
  void compute();
  
  
private:
       
  double m_alpha;    //!< force magnitude
  
};

typedef shared_ptr<ExternalBoundaryPull> ExternalBoundaryPullPtr;

#endif

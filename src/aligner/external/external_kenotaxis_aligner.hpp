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

// cloning and modifying external_shape_align.hpp
// Aligner for boundary particles to orient perpendicular to the tissue boundary

#ifndef __EXTERNAL_KENOTAXIS_ALIGN_HPP__
#define __EXTERNAL_KENOTAXIS_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

struct KenoParameters
{
  double J;
};


class ExternalKenotaxisAlign: public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalKenotaxisAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("J") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (J) specified for cell kenotaxis alignment. Setting it to 1.");
      m_J = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (J) for cell kenotaxis alignment is set to "+param["J"]+".");
      m_J = lexical_cast<double>(param["J"]);
    }
    m_msg->write_config("aligner.external.keno.J",lexical_cast<string>(m_J));
    
    m_type_params = new KenoParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].J = m_J;
    }
  }
  
  virtual ~ExternalKenotaxisAlign()
  {
    delete [] m_type_params;
  }
      
   //! Set pair parameters data for self alignment   
  void set_parameters(pairs_type& self_param)
  {
    map<string,double> param;
    int type;
    
    if (self_param.find("J") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External cell kenotaxis aligner. Setting coupling constant to "+self_param["J"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["J"] = lexical_cast<double>(self_param["J"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External cell kenotaxis aligner. Using default strength ("+lexical_cast<string>(m_J)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["J"] = m_J;
    }
    
    m_type_params[type-1].J = param["J"];
        
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_J;                             //!< Coupling constant
  KenoParameters* m_type_params;        //!< type specific pair parameters 
     
};

typedef shared_ptr<ExternalKenotaxisAlign> ExternalKenotaxisAlignPtr;

#endif


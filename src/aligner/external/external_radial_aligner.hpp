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
 * \file external_radial_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Feb-2018
 * \brief Declaration of ExternalRadialAlign class
 */ 

#ifndef __EXTERNAL_RADIAL_ALIGN_HPP__
#define __EXTERNAL_RADIAL_ALIGN_HPP__

#include <cmath>
#include <vector>

#include "external_aligner.hpp"

using std::make_pair;
using std::sqrt;
using std::vector;

//! Structure that handles parameters for the radial alignment 
struct RadialAlignParameters
{
  double J;      // alignment strength 
  double pos_x;  // x coordinate of the attraction centre
  double pos_y;  // y coordinate of the attraction centre
};


/*! ExternalRadialAlign implements the alignment to the external field that points
 *  radially towards a centre at attraction. This alignment is intended for wound healing
 *  simulations and we implicitly assume radial symmetry (i.e., that we have a planar system).
 *  Alignment filed is given as \f$ \vec H = J\vec e_r \f$, where \f$ \vec e_r \f$ is the
 *  radius vector pointing from the user-specified centre \f$ \vec r_c \f$ towards the particle 
 *  \f$ \vec r_i \f$. In this case, torque is \f$ \vec \tau_i = \vec H \times \vec n_i \f$.
 *  User can specify that the centre of attraction follows the centre of mass. 
 */
class ExternalRadialAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalRadialAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param), 
                                                                            m_use_cm(false)
  {
    int ntypes = m_system->get_ntypes();
    // Set alignment coupling strength J
    if (param.find("J") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No coupling constant (J) specified for cell radial alignment. Setting it to 1.");
      m_J = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global coupling constant (J) for cell radial alignment is set to "+param["J"]+".");
      m_J = lexical_cast<double>(param["J"]);
    }
    m_msg->write_config("aligner.external.radial.J",lexical_cast<string>(m_J));
    
    // Set the center of attraction 
    if (param.find("pos_x") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Radial field aligner. No x component of the attraction centre offset given. Assuming 0.");
      m_pos_x = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Radial field aligner. x component the attraction centre offset set to "+param["pos_x"]+".");
      m_pos_x = lexical_cast<double>(param["pos_x"]);
    }
    m_msg->write_config("aligner.external.radial.pos_x",lexical_cast<string>(m_pos_x));
    if (param.find("pos_y") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Radial field aligner. No y component of the attraction centre offset given. Assuming 0.");
      m_pos_y = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Radial field aligner. y component the attraction centre offset set to "+param["pos_y"]+".");
      m_pos_y = lexical_cast<double>(param["pos_y"]);
    }
    m_msg->write_config("aligner.external.radial.pos_y",lexical_cast<string>(m_pos_y));

    if (param.find("use_centre") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Centre of attraction will be computed relative to the geometric centre of the system.");
      m_use_cm = true;
      m_msg->write_config("aligner.external.radial.use_cm","true");
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Centre of attraction is relative to the fixed origin.");
      m_use_cm = false;
      m_msg->write_config("aligner.external.radial.use_cm","false");
    }
    
    m_type_params = new RadialAlignParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
        m_type_params[i].J = m_J;
        m_type_params[i].pos_x = m_pos_x;
        m_type_params[i].pos_y = m_pos_y; 
    }
  }
  
  virtual ~ExternalRadialAlign()
  {
    delete [] m_type_params;
  }
      
   //! Set pair parameters data for self alignment   
  void set_parameters(pairs_type& self_param)
  {
    map<string,double> param;
    int type;
    
    if (self_param.find("type") == self_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external radial alignment.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    
    if (self_param.find("J") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External radial aligner. Setting coupling strength of the external radial field to "+self_param["J"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["J"] = lexical_cast<double>(self_param["J"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External field aligner. Using default coupling strength J for the external radial field alignment ("+lexical_cast<string>(m_J)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["J"] = m_J;
    }
    m_msg->write_config("aligner.external.radial.type."+self_param["type"]+".J",lexical_cast<string>(param["J"]));

    if (self_param.find("pos_x") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External radial aligner. Setting x component of the centre of attraction offset for the external radial field to "+self_param["pos_x"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["pos_x"] = lexical_cast<double>(self_param["pos_x"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External radial aligner. Using default x component of the centre of attraction offset for the external radial field alignment ("+lexical_cast<string>(m_pos_x)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["pos_x"] = m_pos_x;
    }
    m_msg->write_config("aligner.external.radial.type."+self_param["type"]+".pos_x",lexical_cast<string>(param["pos_x"]));
    
    if (self_param.find("pos_y") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External radial aligner. Setting y component of the centre of attraction offset for the external radial field to "+self_param["pos_y"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["pos_y"] = lexical_cast<double>(self_param["pos_y"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External radial aligner. Using default y component of the centre of attraction offset for the external radial field alignment ("+lexical_cast<string>(m_pos_y)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["pos_x"] = m_pos_x;
    }
    m_msg->write_config("aligner.external.radial.type."+self_param["type"]+".pos_y",lexical_cast<string>(param["pos_y"]));
    
    m_type_params[type-1].J = param["J"];
    m_type_params[type-1].pos_x = param["pos_x"];
    m_type_params[type-1].pos_y = param["pos_y"];
    
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:
       
  double m_J;                              // alignment coupling strength 
  double m_pos_x;                          //!< x component of the offset of the centre of attraction 
  double m_pos_y;                          //!< y component of the offset of the centre of attraction
  bool m_use_cm;                           //!< offset relative to the geometirc centre of the system
  RadialAlignParameters* m_type_params;   //!< type specific alignment  
     
};

typedef shared_ptr<ExternalRadialAlign> ExternalRadialAlignPtr;

#endif


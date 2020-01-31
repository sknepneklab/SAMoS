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
 * \file external_piv_aligner.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Jan-2020
 * \brief Declaration of ExternalPIVAlign class
 */ 

#ifndef __EXTERNAL_PIV_ALIGN_HPP__
#define __EXTERNAL_PIV_ALIGN_HPP__

#include <cmath>
#include <vector>
#include <string>
#include <sstream>

#include <boost/filesystem.hpp>

#include "external_aligner.hpp"

using std::istringstream;
using std::make_pair;
using std::sqrt;
using std::string;
using std::vector;
using namespace boost::filesystem;

//! Structure that handles parameters for the polar alignment 
struct PIVData
{
  vector<vector<double>> piv_data;
};

struct PIVAlignParameters
{
  double gamma;  // Alignment strength
}

/*! ExternalPIVAlign implements the alignment to the velocity field obtained from 
 *  experiments in embryos. Given a local vector vector \f$ \vec v_i \f$. Torque on particle \f$ i \f$ is then 
 *  \f$ \vec \tau_i = \gamma \vec v_i \times \vec n_i \f$, where \f$ \gamma \f$ is a constant. We note that
 *  valued of \f$ \vec v_i \f$ are interpolated between two consecutive time steps (current and the one ahead).
 *  All cells are binned into boxes (given by the PIV resolution). Each cell is aligned to the vector \f$ \vec v_i \f$ 
 *  of the box it falls into.
 *  \note: One could do a more sophisticaed spatial and temporal interpolation. 
 */
class ExternalPIVAlign : public ExternalAlign
{
public:
  
   //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters 
  ExternalPIVAlign(SystemPtr sys, MessengerPtr msg, pairs_type& param) : ExternalAlign(sys,msg,param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("xfile") == param.end())
    {
      m_msg->msg(Messenger::ERROR,"External PIV aligner. No x component of file external PIV given.");
      throw runtime_error("No x component file provided for the external PIV alignment.");
    }
    else
    {
      path p(param["xfile"]);
      if (exists(p))
      {
        this->__read_csv(param["xfile"], m_x);
      }
      else
        throw runtime_error("File " + param["xfile"] + " does not exist.");
    }
    if (param.find("hx") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No x component of external PIV given. Assuming 0.");
      m_hx = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. x component of the external PIV set to "+param["hx"]+".");
      m_hx = lexical_cast<double>(param["hx"]);
    }
    m_msg->write_config("aligner.external.PIV.hx",lexical_cast<string>(m_hx));
    if (param.find("hy") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No y component of external PIV given. Assuming 0.");
      m_hy = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. y component of the external PIV set to "+param["hy"]+".");
      m_hy = lexical_cast<double>(param["hy"]);
    }
    m_msg->write_config("aligner.external.PIV.hy",lexical_cast<string>(m_hy));
    if (param.find("hz") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No z component of external PIV given. Assuming 0.");
      m_hz = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. z component of the external PIV set to "+param["hz"]+".");
      m_hz = lexical_cast<double>(param["hz"]);
    }
    m_msg->write_config("aligner.external.PIV.hz",lexical_cast<string>(m_hz));
    
    m_type_params = new PIVAlignParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_type_params[i].gamma = m_gamma;
    }
  }
  
  virtual ~ExternalPIVAlign()
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
      m_msg->msg(Messenger::ERROR,"Type has not been defined for external PIV alignment.");
      throw runtime_error("Missing key for external alignment parameters.");
    }
    
    type = lexical_cast<int>(self_param["type"]);
    
    
    if (self_param.find("hx") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Setting x component of external filed to "+self_param["hx"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hx"] = lexical_cast<double>(self_param["hx"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Using default x component of external filed ("+lexical_cast<string>(m_hx)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hx"] = m_hx;
    }
    m_msg->write_config("aligner.external.PIV.type."+self_param["type"]+".hx",lexical_cast<string>(param["hx"]));
    if (self_param.find("hy") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Setting y component of external filed to "+self_param["hy"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hy"] = lexical_cast<double>(self_param["hy"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Using default y component of external filed ("+lexical_cast<string>(m_hy)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hy"] = m_hy;
    }
    m_msg->write_config("aligner.external.PIV.type."+self_param["type"]+".hy",lexical_cast<string>(param["hy"]));
    if (self_param.find("hz") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Setting z component of external filed to "+self_param["hz"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["hz"] = lexical_cast<double>(self_param["hz"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Using default z component of external filed ("+lexical_cast<string>(m_hz)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["hz"] = m_hz;
    }
    m_msg->write_config("aligner.external.PIV.type."+self_param["type"]+".hz",lexical_cast<string>(param["hz"]));
    
    m_type_params[type-1].hx = param["hx"];
    m_type_params[type-1].hy = param["hy"];
    m_type_params[type-1].hz = param["hz"];
    
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:

  PIVData  m_x;          // x components of the grid
  PIVData  m_y;          // y components of the grid
  vector<PIVData> m_u;   // Holds x components of the velocity over time    
  vector<PIVData> m_v;   // Holds y components of the velocity over time
  int m_max_time;        // Total number of time steps
  int m_step_gap;        // Number of steps between two consecutive PIV snapshots
  double m_gamma;        // Alignment strength

  void __read_csv(const string &, PIVData &);
};

typedef shared_ptr<ExternalPIVAlign> ExternalPIVAlignPtr;

#endif


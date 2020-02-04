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
#include <algorithm>

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>


#include "external_aligner.hpp"

using std::istringstream;
using std::make_pair;
using std::sqrt;
using std::string;
using std::vector;
using std::sort;
using namespace boost::filesystem;

//! Structure that handles parameters for the polar alignment 
struct PIVData
{
  vector<vector<double>> piv_data;
};

struct PIVAlignParameters
{
  double gamma;  // Alignment strength
};

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

    if (param.find("gamma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No alignment strength given. Assuming 1.0.");
      m_gamma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Setting alignent strength to " + param["gamma"] + ".");
      m_gamma = lexical_cast<double>(param["gamma"]);
    }
    m_msg->write_config("aligner.external.piv.gamma",lexical_cast<string>(m_gamma));
    if (param.find("xmin") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No minimum value of x of the grid set. Assuming -0.5.");
      m_xlow = -0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Setting minimum value of x of the grid to " + param["xmin"] + ".");
      m_xlow = lexical_cast<double>(param["xmin"]);
    }
    m_msg->write_config("aligner.external.piv.xlow",lexical_cast<string>(m_xlow));
    if (param.find("xmax") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No maximum value of x of the grid set. Assuming 0.5.");
      m_xhi = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Setting maximum value of x of the grid to " + param["xmax"] + ".");
      m_xhi = lexical_cast<double>(param["xmax"]);
    }
    m_msg->write_config("aligner.external.piv.xhi",lexical_cast<string>(m_xhi));
    if (param.find("ymin") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No minimum value of y of the grid set. Assuming -0.5.");
      m_ylow = -0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Setting minimum value of y of the grid to " + param["ymin"] + ".");
      m_ylow = lexical_cast<double>(param["ymin"]);
    }
    m_msg->write_config("aligner.external.piv.ylow",lexical_cast<string>(m_ylow));
    if (param.find("ymax") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No maximum value of y of the grid set. Assuming 0.5.");
      m_yhi = 0.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Setting maximum value of y of the grid to " + param["ymax"] + ".");
      m_yhi = lexical_cast<double>(param["ymax"]);
    }
    m_msg->write_config("aligner.external.piv.yhi",lexical_cast<string>(m_yhi));
    if (param.find("piv_path") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No path to PIV files given. Assuming current directory.");
      m_path = "./";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Path of PIV files set to "+param["piv_path"]+".");
      m_path = param["piv_path"];
    }
    m_msg->write_config("aligner.external.piv.piv_path",m_path);
    if (param.find("piv_ext") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. No extension for PIV files given. Assuming csv.");
      m_ext = "csv";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Extension for PIV files set to "+param["piv_ext"]+".");
      m_ext = param["piv_ext"];
    }
    m_msg->write_config("aligner.external.piv.piv_ext",m_ext);
    if (param.find("piv_steps") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"External PIV aligner. Number of simulation steps between two consecutive PIV snapshots not given. Assuming 1000.");
      m_n_piv_steps = 1000;
    }
    else
    {
      m_msg->msg(Messenger::INFO, "External PIV aligner. Number of simulation steps between two consecutive PIV snapshots set to " + param["piv_steps"] + ".");
      m_n_piv_steps = lexical_cast<int>(param["piv_steps"]);
    }

    if (param.find("ufiles") == param.end())
    {
      m_msg->msg(Messenger::ERROR,"External PIV aligner. No base name of x velocity components given.");
      throw runtime_error("No u component file provided for the external PIV alignment.");
    }
    else
    {
      vector<string> filelist;
      this->__get_file_list(param["ufiles"], filelist);
      for (auto f : filelist)
      {
        PIVData piv;
        this->__read_csv(f, piv);
        m_u.push_back(piv);
      }
    }
    
    if (param.find("vfiles") == param.end())
    {
      m_msg->msg(Messenger::ERROR,"External PIV aligner. No base name of y velocity components given.");
      throw runtime_error("No v component file provided for the external PIV alignment.");
    }
    else
    {
      vector<string> filelist;
      this->__get_file_list(param["vfiles"], filelist);
      for (auto f : filelist)
      {
        PIVData piv;
        this->__read_csv(f, piv);
        m_v.push_back(piv);
      }
    }

    // Do some bookkeeping
    m_nx = m_u[0].piv_data[0].size();
    m_ny = m_u[0].piv_data.size();

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
    
    
    if (self_param.find("gamma") != self_param.end())
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Setting alignment strength to "+self_param["gamma"]+" for particle pair of type ("+lexical_cast<string>(type)+").");
      param["gamma"] = lexical_cast<double>(self_param["gamma"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"External PIV aligner. Using default alignment strength ("+lexical_cast<string>(m_gamma)+") for particle pair of types ("+lexical_cast<string>(type)+")");
      param["gamma"] = m_gamma;
    }
    m_msg->write_config("aligner.external.PIV.type."+self_param["type"]+".gamma",lexical_cast<string>(param["gamma"]));
    
    m_type_params[type-1].gamma = param["gamma"];
    
    m_has_params = true;
  }
 
  
  //! Computes "torques"
  void compute();
  
  
private:

  double m_xlow;         // Minimum value of x of the grid
  double m_xhi;          // Maximum value of x of the grid
  double m_ylow;         // Minimum value of y of the grid
  double m_yhi;          // Maximum value of y of the grid
  int m_nx;              // Number of grid points along x
  int m_ny;              // Number of grid points along y

  string m_path;         // Path to PIV data
  string m_ext;          // PIV files extension

  int m_n_piv_steps;     // Number of simulation steps between two consecutive PIV frames.

  vector<PIVData> m_u;   // Holds x components of the velocity over time
  vector<PIVData> m_v;   // Holds y components of the velocity over time
  
  int m_max_time;        // Total number of time steps
  int m_step_gap;        // Number of steps between two consecutive PIV snapshots
  double m_gamma;        // Alignment strength

  int m_piv_frame_counter;
  int m_step_counter;

  PIVAlignParameters* m_type_params;   //!< type specific alignment  

  void __read_csv(const string &, PIVData &);
  void __get_file_list(const string &, vector<string> &);
};

typedef shared_ptr<ExternalPIVAlign> ExternalPIVAlignPtr;

#endif


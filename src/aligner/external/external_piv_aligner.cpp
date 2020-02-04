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
 * \file external_piv_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 11-Apr-2015
 * \brief Declaration of ExternalPIVAlign class
 */ 


#include "external_piv_aligner.hpp"

void ExternalPIVAlign::compute()
{
  int N = m_system->size();
  double gamma = m_gamma;
  double dx = (m_xhi - m_xlow)/(m_nx-1), dy = (m_yhi - m_ylow)/(m_ny - 1);

  double t = static_cast<double>(m_step_counter) / static_cast<double>(m_n_piv_steps);
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_has_params)
      gamma = m_type_params[pi.get_type()-1].gamma;
    if (pi.x > m_xhi || pi.x < m_xlow || pi.y > m_yhi || pi.y < m_ylow)
      throw runtime_error("One of the cells is outside the PIV mesh range. Please confien the system or make a larger PIV mesh");

    int I = int((pi.x - m_xlow) / dx);
    int J = int((pi.y - m_ylow) / dy);

    double u_c = m_u[m_piv_frame_counter].piv_data[I][J], u_n = m_u[m_piv_frame_counter + 1].piv_data[I][J];
    double v_c = m_v[m_piv_frame_counter].piv_data[I][J], v_n = m_v[m_piv_frame_counter + 1].piv_data[I][J];

    double u = t * u_c + (1.0 - t) * u_n;
    double v = t * v_c + (1.0 - t) * v_n;

    double vel = sqrt(u * u + v * v);
    if (vel > 0.0)
    {
      u *= gamma/vel;
      v *= gamma/vel;
    }
    pi.tau_z += u*pi.ny - v*pi.nx;
  
  }
  m_step_counter++;
  if (m_step_counter % m_n_piv_steps == 0)
  {
    m_piv_frame_counter++;
    m_step_counter = 0;
  }
  if (m_piv_frame_counter >= m_u.size())
    throw runtime_error("Simulation time exceeded avaliable PIV snapshots.");
}

// Private member function
void ExternalPIVAlign::__read_csv(const string & fname, PIVData & data)
{
  std::ifstream f(fname.c_str());
  if (f)
  {
    while (f.good())
    {
      string line;
      std::getline(f, line);
      if (line.length() > 0)
      {
        std::istringstream buffer(line);
        string sval;
        vector<double> fline;
        while (getline(buffer, sval, ','))
          fline.push_back(stod(sval));
        data.piv_data.push_back(fline);
      }
    }
  }
  else
  {
    throw runtime_error("Could not open file " + fname + ".");
  }
}

// Get a list of all files matching base name
void ExternalPIVAlign::__get_file_list(const string & base_name, vector<string> & all_matching_files)
{
  const boost::regex filter(base_name + ".*\\."+m_ext);
  directory_iterator end_itr;
  for(directory_iterator i(m_path); i != end_itr; ++i)
  { 
     if( !is_regular_file(i->status()) ) continue;
     boost::smatch what;
     if( !boost::regex_match(i->path().filename().string(), what, filter) ) continue;
     all_matching_files.push_back(m_path+i->path().filename().string());
  } 
  sort(all_matching_files.begin(), all_matching_files.end());
}

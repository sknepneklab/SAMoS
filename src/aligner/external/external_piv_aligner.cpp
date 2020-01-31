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
  double hx = m_hx, hy = m_hy, hz = m_hz;  
  
  for  (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_has_params)
    {
      hx = m_type_params[pi.get_type()-1].hx;
      hy = m_type_params[pi.get_type()-1].hy;
      hz = m_type_params[pi.get_type()-1].hz;
    }
    pi.tau_x += hy*pi.nz - hz*pi.ny;
    pi.tau_y += hz*pi.nx - hx*pi.nz;
    pi.tau_z += hx*pi.ny - hy*pi.nx;
  }
}

// Private member function
void ExternalPIVAlign::__read_csv(const string & fame, PIVData & data)
{
  ifstream f(fname.c_str());
  while (f.good())
  {
    string line;
    getline(f, line);
    if (line.length() > 0)
    {
      istringstream buffer(line);
      string sval;
      vector<double> fline;
      while (getline(buffer, sval, ','))
        fline.push_back(stod(sval));
      data.piv_data.push_back(fline);
    }
  }
}

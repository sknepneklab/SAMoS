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
 * \file population_actomyosin.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 08-Jun-2016
 * \brief Implementation of PopulationActomyosin class.
 */ 

#include "population_actomyosin_molecule.hpp"

/*! This function controls detachment of the entire myosin molecule. This is a
 *  achieved via changing types of myosin head groups from attached (type "A") to detached (type "D") a
 *  and vice versa. 
 *
 *  \note In order to avoid head groups, e.g. changing  their type from "A" to "D" and back to "A" in the 
 *  single time step, we implement both processes in the same function using lists of indices.
 *  \param t current time step
 *  
*/
void PopulationActomyosinMolecule::divide(int t)
{
  if ((m_freq > 0) && (t % m_freq == 0) && (m_detach_rate > 0.0))  // Attempt D to A transition only at certain time steps
  { 
    if (!m_system->group_ok(m_group_name))
    {
      cout << "Before actomyosin molecule detachment: Group info mismatch for group : " << m_group_name << endl;
      throw runtime_error("Group mismatch.");
    }
    double detach_prob = m_detach_rate*m_freq*m_system->get_integrator_step();
    BoxPtr box = m_system->get_box();
    double Lx = box->xhi - box->xlo, Ly = box->yhi - box->ylo;
    int nbin = static_cast<int>(Lx/m_gridsize);    

    myMatrix cmyo(nbin,vector<double>(nbin));
    double dbinx = (1.0/static_cast<double>(nbin))*Lx, dbiny = (1.0/static_cast<double>(nbin))*Ly;

    vector<double> rbinx(nbin);
    vector<double> rbiny(nbin);

    for (int j = 0; j < nbin; j++)
    {
      rbinx[j] = -0.5*Lx + (j+0.5)*dbinx;
      rbiny[j] = -0.5*Ly + (j+0.5)*dbiny;

      for (int k = 0; k < nbin; k++)
        cmyo[j][k]=0.0;
    }
    int Nmol = m_system->number_of_molecules();
    // loop over all molecules in the system
    for (int mol = 0; mol < Nmol; mol++)
    {
      vector<int>& particles = m_system->get_mol_particles(mol);  // get all particle ids that belong to this molecule
      for (unsigned int i = 0; i < particles.size(); i++)   // loop over all particles in the molecule
      {
        Particle& pi = m_system->get_particle(particles[i]);
        // determining the locations of a myosin bead in the grid defined on the system
        if (pi.get_type() != m_type_actin)
        {
          int binx = nbin/2 + floor(pi.x/dbinx);
          int biny = nbin/2 + floor(pi.y/dbiny);
          cmyo[binx][biny] += 1.0;
        }
      }
    }
    // Probability of detachment for a given particle is rate per particle multiplied with time,
    // where time is equal to m_freq*integrator_time_step.   
    for (int mol = 0; mol < Nmol; mol++)
    {
      vector<int>& particles = m_system->get_mol_particles(mol);  // get all particle ids that belong to this molecule
      bool detach = false;
      double myo_comx = 0.0, myo_comy = 0.0, myo_comz = 0.0;
      int imyo = 0;
      int ii = 0;
      int iout = 0;
      for (unsigned int i = 0; i < particles.size(); i++)   // loop over all particles in the molecule
      {
        Particle& pi = m_system->get_particle(particles[i]);
        // determining the locations of a myosin in the grid defined on the system
        if (pi.get_type() != m_type_actin)
        {
          double px = pi.x + pi.ix*Lx, py = pi.y + pi.iy*Ly;    
          if (pi.ix !=0 || pi.iy !=0)  iout++;
          myo_comx += px; 
          myo_comy += py;
          myo_comz += pi.z;
          imyo++;
        }
        if (pi.z > 1.5)    // if at least one of the beads is of type "attached", mark for detachment. NOTE: Why is thi hard coded?
          ii++;
      }
      if (imyo == 18)  // Why is this hard coded?
      {
        myo_comx /= imyo;
        myo_comy /= imyo;
        myo_comz /= imyo;          
        int binx = nbin/2 + floor(myo_comx/dbinx);
        int biny = nbin/2 + floor(myo_comy/dbiny);
        imyo=0;
      }
      if (ii == 18)   // Why is this hard coded?
      {                        
        detach = true;
        ii = 0;
      }
      // if marked for detachment and random number below detachment probability, make transition attached --> detached on all beads
      if (detach)
      {
        double dprob = m_rng->drnd();
        if (dprob < detach_prob)
        {        
          vector<int> to_movex;
          vector<int> to_movey;
          int ind = 0;
          for (int j = 0; j < nbin; j++)
            for (int k = 0; k < nbin; k++)
              if (cmyo[j][k] == 0.0)
              {
                to_movex.push_back(j); 
                to_movey.push_back(k);
                ind++;
              }
          if (ind>0)
          {
            int nmove = static_cast<int>(ind*m_rng->drnd()); 
            int jj = to_movex[nmove], kk = to_movey[nmove];
            double dx = rbinx[jj] - myo_comx, dy = rbiny[kk] - myo_comy, dz = 2.0 - myo_comz; 
            if (iout==1)
            {
              cout << "moving  distance " << dx << "," << dy << "," << dz << endl;          
              cout << "moving  bin " << rbinx[jj] << "," << rbiny[kk] << endl;          
              cout << "myo loc " << myo_comx << "," << myo_comy << "," << myo_comz << endl;          
              cout << "conc at new box " << cmyo[jj][kk] << endl;          
            }
            for (unsigned int ij = 0; ij < particles.size(); ij++)   // loop over all particles in the molecule
            {
              Particle& pj = m_system->get_particle(particles[ij]);
              //cout<<pj.z<<endl;
              double px=pj.x+pj.ix*Lx, py=pj.y+pj.iy*Ly;    
              pj.x = px+dx;                   
              pj.y = py+dy;                   
              pj.z += dz;
              pj.ix = 0;
              pj.iy = 0;
              cmyo[jj][kk]++;
              int binx = nbin/2 + floor(pj.x/dbinx);
              int biny = nbin/2 + floor(pj.y/dbiny);
            }
            //cout<<"myosin moved to box "<<jj<<","<<kk<<"\n"<<endl;
          }
        }
        else
        {
//          cout<<"cannot move myosin "<<dprob<<endl;        
        }
      }
      myo_comx = 0.0;
      myo_comy = 0.0;    
      myo_comz = 0.0;    
      //cout<<"molecule "<<mol<<" done"<<endl;                  
    }
  }
}

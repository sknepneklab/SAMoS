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
 * \file pair_ljrod_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Mar-2015
 * \brief Declaration of PairLJRodPotential class
 */ 

#ifndef __PAIR_LJROD_POTENTIAL_HPP__
#define __PAIR_LJROD_POTENTIAL_HPP__

#include <cmath>
#include <string>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::fabs;
using std::string;

const double LJ_HARD_CORE_DISTANCE = 0.5;   // Distance (in units of sigma) at which we assume potential is effectively infinite

//! Structure that handles parameters for the rod pair potential
struct LJRodParameters
{
  double eps;
  double sigma;
  double rcut;
};

/*! PairLJRodPotential implements the anisotropic Lennard-Jones potential between rods of length \f$ l \f$ 
 *  and radius \f$ \sigma \f$. Details of the interaction potential, force and torques are straightforward
 *  to derive, but lengthy. They are given in a separate set of notes in the docs directory. 
 * 
 *  \note LJRod length is encoded within the potential and not as a part of the Particle class as
 *  for the time being we would like to keep the core of the code work with spherical objects only. If the 
 *  anisotropy turns out to be important, this may change and the Particle class might start carrying
 *  quaternion describing its orientation. This would allow of a uniform interface for implementing 
 *  arbitrary anisotropic potentials. 
 * 
 *  \todo Add brief description of the model here.
 */
class PairLJRodPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairLJRodPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    if (!m_nlist)
    {
      m_msg->msg(Messenger::ERROR,"Lennard Jones rod pair potential requires neighbour list. None given.");
      throw runtime_error("Neighbour list required by Lennard Jones rod potential, but non specified.");
    }
    if (param.find("epsilon") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential depth (epsilon) specified for the rod Lennard Jones pair potential. Setting it to 1.");
      m_eps = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential depth (epsilon) for rod Lennard Jones pair potential is set to "+param["epsilon"]+".");
      m_eps = lexical_cast<double>(param["epsilon"]);
    }
    m_msg->write_config("potential.pair.ljrod.epsilon",lexical_cast<string>(m_eps));
    
    if (param.find("sigma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No particle diameter (sigma) specified for the rod Lennard Jones pair potential. Setting it to 1.");
      m_sigma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global particle diameter (sigma) for rod Lennard Jones pair potential is set to "+param["sigma"]+".");
      m_sigma = lexical_cast<double>(param["sigma"]);
    }
    m_msg->write_config("potential.pair.ljrod.sigma",lexical_cast<string>(m_sigma));
    
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the rod Lennard Jones pair potential. Setting it to 3.0.");
      m_rcut = 3.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for rod Lennard Jones pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    m_msg->write_config("potential.pair.ljrod.rcut",lexical_cast<string>(m_rcut));
    
    if (param.find("WCA") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Rod Lennard-Jones pair potential. Using WCA repulsive potential.");
      m_rcut = 1.122462048309373*m_sigma;
      m_shifted = true;
      m_msg->write_config("potential.pair.ljrod.rcut",lexical_cast<string>(m_rcut));
      m_msg->write_config("potential.pair.ljrod.WCA","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.ljrod.phase_in","true");
    }    
    
    if (m_rcut > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      " is smaller than the rod Lennard-Jones cuttof distance ("+lexical_cast<string>(m_rcut)+
      "). Results will not be reliable.");
    
    if (param.find("shifted") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard-Jones potential shifted to zero at cutoff.");
      m_shifted = true;
      m_msg->write_config("potential.pair.ljrod.shifted","true");
    }
    m_pair_params = new LJRodParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new LJRodParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].eps = m_eps;
        m_pair_params[i][j].sigma = m_sigma;
        m_pair_params[i][j].rcut = m_rcut;
      }
    }
  }
  
  virtual ~PairLJRodPotential()
  {
    for (int i = 0; i < m_ntypes; i++)
      delete [] m_pair_params[i];
    delete [] m_pair_params;
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in rod Lennard-Jones potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in rod Lennard-Jones potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    
    if (pair_param.find("epsilon") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Setting epsilon to "+pair_param["epsilon"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = lexical_cast<double>(pair_param["epsilon"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Using default epsilon ("+lexical_cast<string>(m_eps)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = m_eps;
    }
    m_msg->write_config("potential.pair.ljrod.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".epsilon",lexical_cast<string>(param["epsilon"]));
    if (pair_param.find("sigma") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Setting sigma to "+pair_param["sigma"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = lexical_cast<double>(pair_param["sigma"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Using default sigma ("+lexical_cast<string>(m_sigma)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = m_sigma;
    }
    m_msg->write_config("potential.pair.ljrod.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".sigma",lexical_cast<string>(param["sigma"]));
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Rod Lennard Jones pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }
    if (pair_param.find("WCA") != pair_param.end())
    {
      m_msg->msg(Messenger::WARNING,"Rod Lennard-Jones pair potential. Using WCA repulsive potential.");
      param["rcut"] = 1.122462048309373*param["sigma"];
      m_shifted = true;
      m_msg->write_config("potential.pair.ljrod.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".WCA","true");
      m_msg->write_config("potential.pair.ljrod.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".shifted","true");
    }
    m_msg->write_config("potential.pair.ljrod.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rcut",lexical_cast<string>(param["rcut"]));
       
    
    m_pair_params[type_1-1][type_2-1].eps = param["epsilon"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].eps = param["epsilon"];
    m_pair_params[type_1-1][type_2-1].sigma = param["sigma"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].sigma = param["sigma"];
    m_pair_params[type_1-1][type_2-1].rcut = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rcut = param["rcut"];
    
    if (param["rcut"] > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      "is smaller than the rod Lennard-Jones cuttof distance ("+lexical_cast<string>(m_rcut)+
      ") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+"). Results will not be reliable.");
    
    m_has_pair_params = true;
  }
  
  //! Returns true since rod potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_eps;    //!< potential depth 
  double m_sigma;  //!< particle diameter
  double m_rcut;   //!< cutoff distance
  LJRodParameters** m_pair_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<PairLJRodPotential> PairLJRodPotentialPtr;

#endif

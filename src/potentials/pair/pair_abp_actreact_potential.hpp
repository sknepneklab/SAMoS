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
 * \file pair_abp_actreact_potential.hpp
 * \author Silke Henkes, silkehenkes@gmail.com
 * \date 01-February-2024
 * \brief Declaration of PairABPActReactPotential class 
 */ 

#ifndef __PAIR_ABP_ACTREACT_POTENTIAL_HPP__
#define __PAIR_ABP_ACTREACT_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the abp action reaction pair "potential"
struct ABPActReactParameters
{
  double p;
  double r_int;
};

/*! Pair ABP ActReact is an action-reaction preserving implementation of polar active crawling motion.
* That makes it by symmetry a nematic active stress (!), and a reciprocal interaction (!). Individual polar directions of particles (n) evolve with standard ABP dynamics.
* However, they lead to pair forces F_ij = b 1/2 (p^i n_i - p^j n_j) between particles, where p is the polarisation amplitude (~v0), i.e. pushing against each other.
* The prefactor b goes from 0 at just touching to 1 at r_i = r_j, b = (r_int-r_ij)/r_ij for r_ij < r_int, b = 0 otherwise.
* Parameters: v_0 per type, set as individual parameter (not pair parameter), and r_int, the interaction radius. 
* For results that are as intended, set r_int to the interaction radius of the mechanical potential(s) in the system 
* in particular: r_int = 1+2 eps for soft_attractive, map to re_fact as r_int = 1 + 2*(re_fact-1)  [=1.3 for re_fact=1.15]
* and make use of options use_particle_radii and phase_in. 
 */
class PairABPActReactPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairABPActReactPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("p");
    m_known_params.push_back("r_int");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for ABP action reaction pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in ABP action reaction pair potential.");
    }
    if (param.find("p") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No polarisation force specified for ABP action reaction pair potential. Setting it to 0.1.");
      m_p = 0.1;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global polarisation force specified for ABP action reaction pair potential "+param["p"]+".");
      m_p = lexical_cast<double>(param["p"]);
    }
    m_msg->write_config("potential.pair.abp_actreact.p",lexical_cast<string>(m_p));
    if (param.find("r_int") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential range (r_int) specified for ABP action reaction pair potential. Setting it to 1.3.");
      m_r_int = 1.3;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential range (r_int) for ABP action reaction pair potential is set to "+param["r_int"]+".");
      m_r_int = lexical_cast<double>(param["r_int"]);
    }
    m_msg->write_config("potential.pair.abp_actreact.r_int",lexical_cast<string>(m_r_int));
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"ABP action reaction pair potential is set to use particle radii to control its range. Parameter r_int will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.abp_actreact.use_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"ABP action reaction pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.abp_actreact.phase_in","true");
    }    
    
    m_particle_params = new ABPActReactParameters[m_ntypes];
    std::cout << "Initial constructor parameter values:" << endl;
    for (int i = 0; i < m_ntypes; i++)
    {
      m_particle_params[i].p = m_p;
      m_particle_params[i].r_int = m_r_int;
      std::cout << "type " << i << " p " << m_p << " r_int " << m_r_int << endl;
    }  
    
  }

  virtual ~PairABPActReactPotential()
  {
    //for (int i = 0; i < m_ntypes; i++)
    //  delete [] m_pair_params[i];
    //delete [] m_pair_params;
    delete [] m_particle_params;
  }
  
  //! Set type parameters data for individual particles
  void set_type_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type;
    
    if (pair_param.find("type") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type has not been defined for type specific parameters in vABP action reaction pair potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    
    type = lexical_cast<int>(pair_param["type"]);
    std::cout << "Setting pair parametes of type:" << type << endl;
        
    if (pair_param.find("p") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"ABP action reaction pair potential. Setting polarisation force "+pair_param["p"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["p"] = lexical_cast<double>(pair_param["p"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"ABP action reaction pair potential. Using default polarisation force ("+lexical_cast<string>(m_p)+") for particles of type "+lexical_cast<string>(type)+".");
      param["p"] = m_p;
    }
    m_msg->write_config("potential.pair.abp_actreact.type_"+pair_param["type"]+".push",lexical_cast<string>(param["p"]));
    if (pair_param.find("r_int") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"ABP action reaction pair potential. Setting interaction radious "+pair_param["r_int"]+" for particles of type "+lexical_cast<string>(type)+".");
      param["r_int"] = lexical_cast<double>(pair_param["r_int"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"ABP action reaction pair potential. Using default interaction radius ("+lexical_cast<string>(m_r_int)+") for particles of type "+lexical_cast<string>(type)+".");
      param["r_int"] = m_r_int;
    }
    m_msg->write_config("potential.pair.act_react.type_"+pair_param["type"]+".push",lexical_cast<string>(param["r_int"]));

        
    m_particle_params[type-1].p = param["p"];
    m_particle_params[type-1].r_int = param["r_int"];
    std::cout << "type " << type << " p " << param["p"] << " r_int " << param["r_int"] << endl;
        
    m_has_part_params = true;
  }

  //! Set pair parameters data for pairwise interactions   .. and do nothing. Otherwise the virtual void deities in the base class are unhappy.
  void set_pair_parameters(pairs_type& pair_param)
  {  
    std::cout << "Going through empty pair parameter setting" << endl;
  }
                                                                                                                
  
  
  //! Returns true since soft potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_p;                       //!< polarisation force
  double m_r_int;                       //!< potential range
  bool m_has_part_params;           //!< true if type specific particle parameters are given
  ABPActReactParameters*  m_particle_params;   //!< type specific particle parameters 
     
};

typedef shared_ptr<PairABPActReactPotential> PairABPActReactPotentialPtr;

#endif

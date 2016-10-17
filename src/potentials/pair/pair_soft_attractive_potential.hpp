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
 * \file pair_soft_attractive_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-May-2015
 * \brief Declaration of PairSoftAttractivePotential class 
 */ 

#ifndef __PAIR_SOFT_ATTRACTIVE_POTENTIAL_HPP__
#define __PAIR_SOFT_ATTRACTIVE_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the soft pair potential
struct SoftAttractiveParameters
{
  double k;
  double a;
  double fact;
};

/*! PairSoftAttractivePotential implements the pair interaction which is soft in nature with 
 *  an attractive part and soft repulsive core.
 *  Potential is given as \f$ U_{soft}\left(r_{ij}\right) = \frac{k}{2} \left(r_{ij} - r_e \right)^2 \f$, if
 *  \f$ r_{ij} \le r_{tunr} \f$,  \f$ U_{soft} =  -\frac{k}{2} \left(r_{ij} - r_f \right)^2 \f$ (with 
 *  \f$ r_f = 2 r_{turn} - r_e \f$)if \f$ r_e < r_{ij} < r_f \f$ and zero otherwise.
 *  \f$ k \f$ is the potential strength. Often, it is convenient to set \f$ r_e = a_i + a_j \f$, with
 *  \f$ a_i \f$ and \f$ a_j \f$ being the radii of the two particles. Finally, \f$ r_{ij} \f$ is the 
 *  Euclidean interparticle distance.
 *  \f$ r_{turn} = fr_e \f$, where \f$ f \f$ is a parameter given by the user.
 */
class PairSoftAttractivePotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairSoftAttractivePotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("k");
    m_known_params.push_back("a");
    m_known_params.push_back("re_fact");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for soft attractive pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in soft attractive potential.");
    }
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (k) specified for soft attractive pair potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (k) for soft attractive pair potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.pair.soft_attractive.k",lexical_cast<string>(m_k));
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential range (a) specified for soft attractive pair potential. Setting it to 2.");
      m_a = 2.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential range (a) for soft pair attractive potential is set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
    }
    m_msg->write_config("potential.pair.soft_attractive.a",lexical_cast<string>(m_a));
    if (param.find("re_fact") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No turn factor (re_fact) specified for soft attractive pair potential. Setting it to 1.25.");
      m_fact = 1.25;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global turn factor (re_fact) for soft pair attractive potential is set to "+param["re_fact"]+".");
      m_fact = lexical_cast<double>(param["re_fact"]);
    }
    m_msg->write_config("potential.pair.soft_attractive.re_fact",lexical_cast<string>(m_fact));
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Soft attractive pair potential is set to use particle radii to control its range. Parameter a will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.soft_attractive.use_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Soft pair attractive potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.soft_attractive.phase_in","true");
    }    
    
    m_pair_params = new SoftAttractiveParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new SoftAttractiveParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].k = m_k;
        m_pair_params[i][j].a = m_a;
        m_pair_params[i][j].fact = m_fact;
      }
    }
    
  }
  
  virtual ~PairSoftAttractivePotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in soft attractive potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in soft attractive potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("k") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Soft attractive pair potential. Setting strength to "+pair_param["k"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = lexical_cast<double>(pair_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Soft pair potential. Using default strength ("+lexical_cast<string>(m_k)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = m_k;
    }
    m_msg->write_config("potential.pair.soft_attractive.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".push",lexical_cast<string>(param["k"]));
    if (pair_param.find("a") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Soft attractive pair potential. Setting range to "+pair_param["a"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = lexical_cast<double>(pair_param["a"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Soft attractive pair potential. Using default range ("+lexical_cast<string>(m_a)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = m_a;
    }
    m_msg->write_config("potential.pair.soft_attractive.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".a",lexical_cast<string>(param["a"]));
    if (pair_param.find("re_fact") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Soft attractive pair potential. Setting turn factor "+pair_param["re_fact"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["re_fact"] = lexical_cast<double>(pair_param["re_fact"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Soft attractive pair potential. Using default turn factor ("+lexical_cast<string>(m_fact)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["re_fact"] = m_fact;
    }
    m_msg->write_config("potential.pair.soft_attractive.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".re_fact",lexical_cast<string>(param["re_fact"]));  

    
    m_pair_params[type_1-1][type_2-1].k = param["k"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].k = param["k"];
    m_pair_params[type_1-1][type_2-1].a = param["a"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].a = param["a"];
    m_pair_params[type_1-1][type_2-1].fact = param["re_fact"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].fact = param["re_fact"];
    
    m_has_pair_params = true;
  }
  
  //! Returns true since soft potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_k;                       //!< potential strength
  double m_a;                       //!< potential range
  double m_fact;                    //!< factor determining position of the potential turn r_turn
  SoftAttractiveParameters** m_pair_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<PairSoftAttractivePotential> PairSoftAttractivePotentialPtr;

#endif

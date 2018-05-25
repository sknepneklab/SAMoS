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
 * \file pair_yukawa_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-May-2018
 * \brief Declaration of PairYukawaPotential class
 */ 

#ifndef __PAIR_YUKAWA_POTENTIAL_HPP__
#define __PAIR_YUKAWA_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the Yukawa pair potential
struct YukawaParameters
{
  double g;
  double kappa;
  double rcut;
};

/*! PairYukawaPotential implements standard Yukawa potential. Potential is given as
 *  \f$ U_{Yukawa}\left(r_{ij}\right) = g\frac{e^{-\kappa r_{ij}}}{r_{ij}} \f$,
 *  where \f$ g \f$ is the potential strength, \f$ \kappa \f$ is the inverse potentail range and \f$ r_{ij} \f$ is the 
 *  interparticle distance.
 */
class PairYukawaPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (alpha and sigma)
  PairYukawaPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("g");
    m_known_params.push_back("kappa");
    m_known_params.push_back("rcut");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for Yukawa pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in Yukawa potential.");
    }
    if (param.find("g") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (g) specified for Yukawa pair potential. Setting it to 1.");
      m_g = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (g) for Yukawa pair potential is set to "+param["g"]+".");
      m_g = lexical_cast<double>(param["g"]);
    }
    m_msg->write_config("potential.pair.Yukawa.g",lexical_cast<string>(m_g));
    
    if (param.find("kappa") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No inverse potentail range (kappa) specified for Yukawa pair potential. Setting it to 1.");
      m_kappa = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global inverse potentail range (kappa) for Yukawa pair potential is set to "+param["kappa"]+".");
      m_kappa = lexical_cast<double>(param["kappa"]);
    }
    m_msg->write_config("potential.pair.Yukawa.kappa",lexical_cast<string>(m_kappa));
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the Yukawa pair potential. Setting it to 3.0.");
      m_rcut = 3.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for Yukawa pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    m_msg->write_config("potential.pair.yukawa.rcut",lexical_cast<string>(m_rcut));
    
    m_pair_params = new YukawaParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new YukawaParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].g = m_g;
        m_pair_params[i][j].kappa = m_kappa;
        m_pair_params[i][j].rcut = m_rcut;
      }
    }
  }
  
  virtual ~PairYukawaPotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in Yukawa potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in Yukawa potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("g") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Setting strength to "+pair_param["g"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["g"] = lexical_cast<double>(pair_param["g"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Using default strength ("+lexical_cast<string>(m_g)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["g"] = m_g;
    }
    if (pair_param.find("kappa") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Setting kappa to "+pair_param["kappa"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["kappa"] = lexical_cast<double>(pair_param["kappa"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Using default kappa ("+lexical_cast<string>(m_kappa)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["kappa"] = m_kappa;
    }
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Yukawa pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }
    m_msg->write_config("potential.pair.yukawa.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rcut",lexical_cast<string>(param["rcut"]));
    m_msg->write_config("potential.pair.Yukawa.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".g",lexical_cast<string>(param["g"]));
    m_msg->write_config("potential.pair.Yukawa.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".kappa",lexical_cast<string>(param["kappa"]));
        
    m_pair_params[type_1-1][type_2-1].g = param["g"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].g = param["g"];
    m_pair_params[type_1-1][type_2-1].kappa = param["kappa"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].kappa = param["kappa"];
    m_pair_params[type_1-1][type_2-1].rcut = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rcut = param["rcut"];
    
    m_has_pair_params = true;
  }
  
  //! Returns false since Yukawa potential does not need neighbour list
  bool need_nlist() { return false; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
        
  double m_g;           //!< potential strength
  double m_kappa;       //!< inverse potential range
  double m_rcut;        //!< cutoff distance
  YukawaParameters** m_pair_params;   //!< type specific pair parameters 
    
};

typedef shared_ptr<PairYukawaPotential> PairYukawaPotentialPtr;

#endif

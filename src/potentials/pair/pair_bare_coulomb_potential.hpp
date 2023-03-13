/* ***************************************************************************
 *
 *  Copyright (C) 2013-2022 University of Dundee
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
 * \file pair_coulomb_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Jun-2022
 * \brief Declaration of PairBareCoulombPotential class
 */ 

#ifndef __PAIR_BARE_COULOMB_POTENTIAL_HPP__
#define __PAIR_BARE_COULOMB_POTENTIAL_HPP__

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#else
#ifndef __OVERRIDE_OMP__
#define __OVERRIDE_OMP__
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num()  { return 0; }
inline omp_int_t omp_get_num_threads() { return 1; }
inline omp_int_t omp_get_max_threads() { return 1; }
#endif
#endif


#ifdef __AVX2__
#include <immintrin.h>
#endif

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;


//! Structure that handles parameters for the bare Coulomb pair potential
struct BareCoulombParameters
{
  double alpha;
};

/*! PairBareCoulombPotential implements standard Coulomb potential. Note that there is no hard-core term
 *  and the simulation will fail for system with particles with opposite charges. The potential is given as
 *  \f$ U_{Coul}\left(r_{ij}\right) = \frac{\alpha}{r_{ij}} \f$,
 *  where \f$ \alpha \f$ is the potential strength,  and \f$ r_{ij} \f$ is the 
 *  inter-particle distance.
 */
class PairBareCoulombPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (alpha and sigma)
  PairBareCoulombPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    m_known_params.push_back("alpha");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for bare Coulomb pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in Coulomb potential.");
    }
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (alpha) specified for bare Coulomb pair potential. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (alpha) for bare Coulomb pair potential is set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.pair.barecoulomb.alpha",lexical_cast<string>(m_alpha));
    
    m_pair_params = new BareCoulombParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new BareCoulombParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].alpha = m_alpha;
      }
    }
  }
  
  virtual ~PairBareCoulombPotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in bare Coulomb potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in bare Coulomb potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Coulomb pair potential. Setting strength to "+pair_param["alpha"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Coulomb pair potential. Using default strength ("+lexical_cast<string>(m_alpha)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = m_alpha;
    }
    m_msg->write_config("potential.pair.barecoulomb.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".alpha",lexical_cast<string>(param["alpha"]));
        
    m_pair_params[type_1-1][type_2-1].alpha = param["alpha"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].alpha = param["alpha"];
    
    m_has_pair_params = true;
  }
  
  //! Returns false since Coulomb potential does not need neighbour list
  bool need_nlist() { return false; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
        
  double m_alpha;       //!< potential strength
  BareCoulombParameters** m_pair_params;   //!< type specific pair parameters 
    
};

typedef shared_ptr<PairBareCoulombPotential> PairBareCoulombPotentialPtr;

#endif

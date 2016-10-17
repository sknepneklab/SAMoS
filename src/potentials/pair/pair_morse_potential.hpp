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
 * \file pair_morse_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 04-Oct-2014
 * \brief Declaration of PairMorsePotential class
 */ 

#ifndef __PAIR_MORSE_POTENTIAL_HPP__
#define __PAIR_MORSE_POTENTIAL_HPP__

#include "pair_potential.hpp"

using std::sqrt;
using std::exp;

//! Structure that handles parameters for the Morse pair potential
struct MorseParameters
{
  double D;
  double a;
  double re;
  double rcut;
};

/*! PairMorsePotential implements standard MOrse potential 
 * \f$ U_{Morse}\left(r_{ij}\right) = D \left(1 - e^{-a(r-r_e)} \right)^2 \f$,
 *  where \f$ D \f$ is the depth of the potential, \f$ a \f$ controls its width and \f$ r_e \f$ is the 
 *  the position of the minimum
 */
class PairMorsePotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (epsilon, sigma, and rcut)
  PairMorsePotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    if (!m_nlist)
    {
      m_msg->msg(Messenger::ERROR,"Morse pair potential requires neighbour list. None given.");
      throw runtime_error("Neighbour list required by Morse potential, but non specified.");
    }
    m_known_params.push_back("D");
    m_known_params.push_back("a");
    m_known_params.push_back("re");
    m_known_params.push_back("rcut");
    m_known_params.push_back("use_particle_radii");
    m_known_params.push_back("phase_in");
    m_known_params.push_back("shifted");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for Morse pair potential.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in Morse potential.");
    }
    if (param.find("D") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential depth (D) specified for the Morse pair potential. Setting it to 1.");
      m_D = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential depth (D) for Morse pair potential is set to "+param["D"]+".");
      m_D = lexical_cast<double>(param["D"]);
    }
    m_msg->write_config("potential.pair.morse.D",lexical_cast<string>(m_D));
    
    if (param.find("a") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No width parameter (a) specified for Morse pair potential. Setting it to 1.");
      m_a = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global width parameter (a) for morse pair potential is set to "+param["a"]+".");
      m_a = lexical_cast<double>(param["a"]);
    }
    m_msg->write_config("potential.pair.morse.a",lexical_cast<string>(m_a));
    
    if (param.find("re") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No position of the minimum (re) specified for Morse pair potential. Setting it to 1.");
      m_re = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global position of the minimum (re) for morse pair potential is set to "+param["re"]+".");
      m_re = lexical_cast<double>(param["re"]);
    }
    m_msg->write_config("potential.pair.morse.re",lexical_cast<string>(m_re));
    
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the Morse pair potential. Setting it to 2.5*re.");
      m_rcut = 2.5*m_re;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for Lennard Jones pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    m_msg->write_config("potential.pair.morse.rcut",lexical_cast<string>(m_rcut));
    
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Morse pair potential is set to use particle radii to control its range. Parameter re will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.morse.use_particle_radii","true");
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.morse.phase_in","true");
    }    
    
    if (m_rcut > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      " is smaller than the Morse cuttof distance ("+lexical_cast<string>(m_rcut)+
      "). Results will not be reliable.");
    
    if (param.find("shifted") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse potential shifted to zero at cutoff.");
      m_shifted = true;
    }
    m_pair_params = new MorseParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new MorseParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].D = m_D;
        m_pair_params[i][j].a = m_a;
        m_pair_params[i][j].re = m_re;
        m_pair_params[i][j].rcut = m_rcut;
      }
    }
  }
  
  virtual ~PairMorsePotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in Morse potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters inMorse potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    
    if (pair_param.find("D") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Setting D to "+pair_param["D"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["D"] = lexical_cast<double>(pair_param["D"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Using default D ("+lexical_cast<string>(m_D)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["D"] = m_D;
    }
    m_msg->write_config("potential.pair.morse.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".D",lexical_cast<string>(param["D"]));
    if (pair_param.find("a") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Setting a to "+pair_param["a"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = lexical_cast<double>(pair_param["a"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Using default a ("+lexical_cast<string>(m_a)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["a"] = m_a;
    }
    m_msg->write_config("potential.pair.morse.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".a",lexical_cast<string>(param["a"]));
    if (pair_param.find("re") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Setting re to "+pair_param["re"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["re"] = lexical_cast<double>(pair_param["re"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Using default re ("+lexical_cast<string>(m_re)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["re"] = m_re;
    }
    m_msg->write_config("potential.pair.morse.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".re",lexical_cast<string>(param["re"]));
    
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Morse pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }
    m_msg->write_config("potential.pair.morse.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rcut",lexical_cast<string>(param["rcut"]));
        
    m_pair_params[type_1-1][type_2-1].D = param["D"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].D = param["D"];
    m_pair_params[type_1-1][type_2-1].a = param["a"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].a = param["a"];
    m_pair_params[type_1-1][type_2-1].re = param["re"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].re = param["re"];
    m_pair_params[type_1-1][type_2-1].rcut = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rcut = param["rcut"];
    
    if (param["rcut"] > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      "is smaller than the Morse cuttof distance ("+lexical_cast<string>(m_rcut)+
      ") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+"). Results will not be reliable.");
    
    m_has_pair_params = true;
  }
  
  //! Returns true since Morse potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_D;    //!< potential depth 
  double m_a;    //!< potential well depth parameter
  double m_re;   //!< position of the potential minimum
  double m_rcut;   //!< cutoff distance
  MorseParameters** m_pair_params;   //!< type specific pair parameters 
    
};

typedef shared_ptr<PairMorsePotential> PairMorsePotentialPtr;

#endif

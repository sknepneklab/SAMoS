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
 * \file pair_gaussian_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 26-Sept-2017
 * \brief Declaration of PairGaussianPotential class
 */ 

#ifndef __PAIR_GAUSSIAN_POTENTIAL_HPP__
#define __PAIR_GAUSSIAN_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the Gaussian pair potential
struct GaussianParameters
{
  double A;
  double B;
  double alpha;
  double beta;
  double rA;
  double rB;
  double rcut;
};

/*! PairGaussianPotential implements a potential with two Gaussians shifted with respect to each other 
 * \f$ U_{Gauss}\left(r_{ij}\right) = A e^{-\alpha(r-r_A)^2} + B e^{-\beta(r-r_B)} \f$,
 *  where \f$ A \f$ and \f$ B \f$ control strengths of two parts, \f$ \alpha \f$ and \f$ \beta \f$
 *  their widths, and \f$ r_A \f$ and \f$ r_b \f$ position of the Gaussian peaks.
 */
class PairGaussianPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (epsilon, sigma, and rcut)
  PairGaussianPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    if (!m_nlist)
    {
      m_msg->msg(Messenger::ERROR,"Gaussian pair potential requires neighbour list. None given.");
      throw runtime_error("Neighbour list required by Gaussian potential, but non specified.");
    }
    if (param.find("A") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No first potential depth (A) specified for the Gaussian pair potential. Setting it to 2.");
      m_A = 2.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global first potential depth (A) for Gaussian pair potential is set to "+param["A"]+".");
      m_A = lexical_cast<double>(param["A"]);
    }
    m_msg->write_config("potential.pair.gaussian.A",lexical_cast<string>(m_A));
    if (param.find("B") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No second potential depth (B) specified for the Gaussian pair potential. Setting it to -1.");
      m_B = -1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global second potential depth (B) for Gaussian pair potential is set to "+param["B"]+".");
      m_B = lexical_cast<double>(param["B"]);
    }
    m_msg->write_config("potential.pair.gaussian.B",lexical_cast<string>(m_B));
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No first potential width (alpha) specified for the Gaussian pair potential. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global first potential width (alpha) for Gaussian pair potential is set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("potential.pair.gaussian.alpha",lexical_cast<string>(m_alpha));
    if (param.find("beta") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No second potential width (beta) specified for the Gaussian pair potential. Setting it to 3.");
      m_beta = 3.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global second potential width (beta) for Gaussian pair potential is set to "+param["beta"]+".");
      m_beta = lexical_cast<double>(param["beta"]);
    }
    m_msg->write_config("potential.pair.gaussian.beta",lexical_cast<string>(m_beta));
    if (param.find("rA") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No first potential peak position (rA) specified for the Gaussian pair potential. Setting it to 0.");
      m_rA = 0.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global first potential peak position (rA) for Gaussian pair potential is set to "+param["rA"]+".");
      m_rA = lexical_cast<double>(param["rA"]);
    }
    m_msg->write_config("potential.pair.gaussian.rA",lexical_cast<string>(m_rA));
    if (param.find("rB") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No second potential peak position (rB) specified for the Gaussian pair potential. Setting it to 1.5.");
      m_rB = 1.5;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global second potential peak position (rB) for Gaussian pair potential is set to "+param["rB"]+".");
      m_rB = lexical_cast<double>(param["rB"]);
    }
    m_msg->write_config("potential.pair.gaussian.rB",lexical_cast<string>(m_rB));
    
    if (param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the Gaussian pair potential. Setting it to 3.0.");
      m_rcut = 3.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for Gaussian pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    m_msg->write_config("potential.pair.gaussian.rcut",lexical_cast<string>(m_rcut));
    
    if (param.find("use_particle_radii") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Gaussian pair potential is set to use particle radii to control its range. Parameter rB will be ignored.");
      m_use_particle_radii = true;
      m_msg->write_config("potential.pair.gaussian.use_particle_radii","true");
    }
    
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
      m_msg->write_config("potential.pair.gaussian.phase_in","true");
    }   
    
    m_pair_params = new GaussianParameters*[m_ntypes];
    for (int i = 0; i < m_ntypes; i++)
    {
      m_pair_params[i] = new GaussianParameters[m_ntypes];
      for (int j = 0; j < m_ntypes; j++)
      {
        m_pair_params[i][j].A = m_A;
        m_pair_params[i][j].B = m_A;
        m_pair_params[i][j].alpha = m_alpha;
        m_pair_params[i][j].beta = m_beta;
        m_pair_params[i][j].rA = m_rA;
        m_pair_params[i][j].rB = m_rB;
        m_pair_params[i][j].rcut = m_rcut;
      }
    }
    
    if (m_rcut > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      " is smaller than the Gaussian cuttof distance ("+lexical_cast<string>(m_rcut)+
      "). Results will not be reliable.");
    
    if (param.find("shifted") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian potential shifted to zero at cutoff.");
      m_shifted = true;
    }
  }
  
  virtual ~PairGaussianPotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in Gaussian potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in Gaussian potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    
    if (pair_param.find("A") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting A to "+pair_param["A"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["A"] = lexical_cast<double>(pair_param["A"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default A ("+lexical_cast<string>(m_A)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["A"] = m_A;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".A",lexical_cast<string>(param["A"]));
    if (pair_param.find("B") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting B to "+pair_param["B"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["B"] = lexical_cast<double>(pair_param["B"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default B ("+lexical_cast<string>(m_B)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["B"] = m_B;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".B",lexical_cast<string>(param["B"]));
    if (pair_param.find("alpha") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting alpha to "+pair_param["alpha"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = lexical_cast<double>(pair_param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default alpha ("+lexical_cast<string>(m_alpha)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["alpha"] = m_alpha;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".alpha",lexical_cast<string>(param["alpha"]));
    if (pair_param.find("beta") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting beta to "+pair_param["beta"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["beta"] = lexical_cast<double>(pair_param["beta"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default beta ("+lexical_cast<string>(m_beta)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["beta"] = m_beta;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".beta",lexical_cast<string>(param["beta"]));
    if (pair_param.find("rA") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting rA to "+pair_param["rA"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rA"] = lexical_cast<double>(pair_param["rA"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default rA ("+lexical_cast<string>(m_rA)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rA"] = m_rA;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rA",lexical_cast<string>(param["rA"]));
    if (pair_param.find("rB") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting rB to "+pair_param["rB"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rB"] = lexical_cast<double>(pair_param["rB"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default rB ("+lexical_cast<string>(m_rB)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rB"] = m_rB;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rB",lexical_cast<string>(param["rB"]));
    
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Gaussian pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_rcut;
    }
    m_msg->write_config("potential.pair.gaussian.type_"+pair_param["type_1"]+"_and_type_"+pair_param["type_2"]+".rcut",lexical_cast<string>(param["rcut"]));
    
    m_pair_params[type_1-1][type_2-1].A = param["A"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].A = param["A"];
    m_pair_params[type_1-1][type_2-1].B = param["B"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].B = param["B"];
    m_pair_params[type_1-1][type_2-1].alpha = param["alpha"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].alpha = param["alpha"];
    m_pair_params[type_1-1][type_2-1].beta = param["beta"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].beta = param["beta"];
    m_pair_params[type_1-1][type_2-1].rA = param["rA"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rA = param["rA"];
    m_pair_params[type_1-1][type_2-1].rB = param["rB"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rB = param["rB"];
    m_pair_params[type_1-1][type_2-1].rcut = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].rcut = param["rcut"];
    
    if (param["rcut"] > m_nlist->get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<string>(m_nlist->get_cutoff())+
      "is smaller than the Gaussian cuttof distance ("+lexical_cast<string>(m_rcut)+
      ") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+"). Results will not be reliable.");
    
    m_has_pair_params = true;
  }
  
  //! Returns true since Gaussian potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_A;      //!< height of first potential
  double m_B;      //!< height of second potential
  double m_alpha;  //!< width of first potential
  double m_beta;   //!< width of second potential
  double m_rA;     //!< position of fist peak
  double m_rB;     //!< position of second peak
  double m_rcut;   //!< cutoff distance
  GaussianParameters** m_pair_params;   //!< type specific pair parameters   
};

typedef shared_ptr<PairGaussianPotential> PairGaussianPotentialPtr;

#endif

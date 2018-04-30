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
 * \file integrator_Sepulveda.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 30-Apr-2018
 * \brief Declaration of IntegratorSepulveda class
 */ 

#ifndef __INTEGRATOR_SEPULVEDA_H__
#define __INTEGRATOR_SEPULVEDA_H__

#include <cmath>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;
using std::exp;
using std::string;
using std::vector;

using boost::algorithm::to_lower_copy;

/*! IntegratorSepulveda class implements Sepulveda dynamics for the particle position.
    Particle director will not be integrated. 
    refeference: N. Sepulveda, L. Petitjean, O. Cochet, E. Grasland-Mongrain, P. Silberzan
    V. Hakim, Collective Cell Motion in an Epithelial Sheet Can Be Quantitatively Described by a 
    Stochastic Interacting Particle Model, PLoS Comput Biol 9(3): e1002944.
*/
class IntegratorSepulveda : public Integrator
{
public:
  
  //! Constructor
  //! \param sys Pointer to a System object containing all particles
  //! \param msg Internal message handler
  //! \param pot Pairwise and external interaction handler
  //! \param align Pairwise and external alignment handler
  //! \param nlist Neighbour list object
  //! \param cons Enforces constraints to the manifold surface
  //! \param temp Temperature control object
  //! \param param Contains information about all parameters 
  IntegratorSepulveda(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("alpha");
    m_known_params.push_back("seed");
    m_known_params.push_back("sigma");
    m_known_params.push_back("tau");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for Sepulveda integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in Sepulveda integrator.");
    }
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Sepulveda dynamics integrator for particle position. Friction constant not set. Using default value 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Sepulveda dynamics integrator for particle position. Setting friction constant to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    m_msg->write_config("integrator.Sepulveda.alpha",lexical_cast<string>(m_alpha));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Sepulveda dynamics integrator for particle position. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.Sepulveda.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Sepulveda dynamics integrator for particle position. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.Sepulveda.seed",param["seed"]);
    }
    if (param.find("sigma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Sepulveda dynamics integrator for particle position. Noise amplitude. Using default value 1.");
      m_sigma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Sepulveda dynamics integrator for particle position. Noise amplitude "+param["sigma"]+".");
      m_sigma = lexical_cast<double>(param["sigma"]);
    }
    m_msg->write_config("integrator.Sepulveda.sigma",lexical_cast<string>(m_sigma));
    if (param.find("tau") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Sepulveda dynamics integrator for particle position. Correlation time not set. Using default value 1.");
      m_tau = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Sepulveda dynamics integrator for particle position. Setting correlation time to "+param["tau"]+".");
      m_tau = lexical_cast<double>(param["tau"]);
    }
    m_msg->write_config("integrator.Sepulveda.tau",lexical_cast<string>(m_tau));
    // Set initial etas
    for (int i = 0; i < m_system->get_group(m_group_name)->get_size(); i++)
    {
      m_eta_x.push_back(m_rng->drnd());
      m_eta_y.push_back(m_rng->drnd());
      m_eta_z.push_back(m_rng->drnd());
    }
  }

  virtual ~IntegratorSepulveda()
  {
    m_eta_x.clear();  m_eta_y.clear(); m_eta_z.clear();
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;             //!< Random number generator 
  double  m_alpha;           //!< Friction constant 
  double  m_sigma;           //!< Noise amplitude 
  double  m_tau;             //!< Correlation time 
  vector<double>  m_eta_x;      //!< Normally distributed random numbers for BBK integrator
  vector<double>  m_eta_y;      //!< Normally distributed random numbers for BBK integrator
  vector<double>  m_eta_z;      //!< Normally distributed random numbers for BBK integrator


};

typedef shared_ptr<IntegratorSepulveda> IntegratorSepulvedaPtr;

#endif

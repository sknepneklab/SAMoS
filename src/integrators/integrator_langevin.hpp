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
 * \file integrator_langevin.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 29-Dec-2015
 * \brief Declaration of IntegratorLangevin class
 */ 

#ifndef __INTEGRATOR_LANGEVIN_H__
#define __INTEGRATOR_LANGEVIN_H__

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

/*! IntegratorLangevin class implements Langevin dynamics for the particle position.
    Particle director will not be integrated. 
*/
class IntegratorLangevin : public Integrator
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
  IntegratorLangevin(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param)
  { 
    m_known_params.push_back("gamma");
    m_known_params.push_back("seed");
    m_known_params.push_back("method");
    string param_test = this->params_ok(param);
    if (param_test != "")
    {
      m_msg->msg(Messenger::ERROR,"Parameter \""+param_test+"\" is not a valid parameter for langevin integrator.");
      throw runtime_error("Unknown parameter \""+param_test+"\" in langevin integrator.");
    }
    if (param.find("gamma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator for particle position. Friction constant not set. Using default value 1.");
      m_gamma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator for particle position. Setting friction constant to "+param["gamma"]+".");
      m_gamma = lexical_cast<double>(param["gamma"]);
    }
    m_msg->write_config("integrator.langevin.gamma",lexical_cast<string>(m_gamma));
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator for particle position. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
      m_msg->write_config("integrator.langevin.seed",lexical_cast<string>(0));
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator for particle position. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
      m_msg->write_config("integrator.langevin.seed",param["seed"]);
    }
    if (param.find("method") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Langevin dynamics integrator for particle position. No integration method specified. Using default BAOAB.");
      m_method = "baoab";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Langevin dynamics integrator for particle position. Using "+param["method"]+" integration method.");
      m_method = to_lower_copy(param["method"]);
      if (m_method == "bbk")
      {
        int N = m_system->get_group(m_group_name)->get_size();
        for (int i = 0; i < N; i++)
        {
          m_Rx.push_back(m_rng->gauss_rng(1.0));
          m_Ry.push_back(m_rng->gauss_rng(1.0));
          m_Rz.push_back(m_rng->gauss_rng(1.0));
        }
      }
    }
    m_msg->write_config("integrator.langevin.method",m_method);
  }

  virtual ~IntegratorLangevin()
  {
    m_Rx.clear();  m_Ry.clear(); m_Rz.clear();
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:
  
  RNGPtr  m_rng;             //!< Random number generator 
  double  m_gamma;           //!< Friction constant 
  string  m_method;          //!< Specify integration method (BAOAB, SPV, or BBK)
  vector<double>  m_Rx;      //!< Normally distributed random numbers for BBK integrator
  vector<double>  m_Ry;      //!< Normally distributed random numbers for BBK integrator
  vector<double>  m_Rz;      //!< Normally distributed random numbers for BBK integrator

  //! Integrate using BAOAB method (defualt)
  void integrate_baoab();

  //! Integrate using SPV method
  void integrate_spv();

  //! Integrate using BBK method
  void integrate_bbk();


};

typedef shared_ptr<IntegratorLangevin> IntegratorLangevinPtr;

#endif

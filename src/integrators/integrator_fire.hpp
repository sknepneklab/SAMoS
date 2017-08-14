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
 * \file integrator_fire.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 09-Aug-2017
 * \brief Declaration of IntegratorFIRE class
 */ 

#ifndef __INTEGRATOR_FIRE_H__
#define __INTEGRATOR_FIRE_H__

#include <cmath>

#include "integrator.hpp"
#include "rng.hpp"

using std::sqrt;
using std::fabs;
using std::min;

/*! IntegratorFIRE class handles FIRE minimization. 
 *  \note No activity. Just minimization. 
*/
class IntegratorFIRE : public Integrator
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
  IntegratorFIRE(SystemPtr sys, MessengerPtr msg, PotentialPtr pot, AlignerPtr align, NeighbourListPtr nlist,  ConstrainerPtr cons, ValuePtr temp, pairs_type& param) : Integrator(sys, msg, pot, align, nlist, cons, temp, param),
                                                                                                                                                                        m_converged(false),
                                                                                                                                                                        m_old_energy(1e15),
                                                                                                                                                                        m_last_neg(0)
  { 
    m_msg->write_config("integrator.fire","");
    if (param.find("alpha") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator intital alpha set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator alpha not set. Using default value of 0.1.");
      m_alpha = 0.1;
    }
    m_alpha_init = m_alpha;
    m_msg->write_config("integrator.FIRE.alpha",lexical_cast<string>(m_alpha));
    if (param.find("alpha_scale") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator alpha_scale set to "+param["alpha_scale"]+".");
      m_f_alpha = lexical_cast<double>(param["alpha_scale"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator alpha_scale not set. Using default value of 0.99.");
      m_f_alpha = 0.99;
    }
    m_msg->write_config("integrator.FIRE.f_alpha",lexical_cast<string>(m_f_alpha));
    if (param.find("f_tolerance") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator f_tolerance set to "+param["f_tolerance"]+".");
      m_F_tol = lexical_cast<double>(param["f_tolerance"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator f_tolerance not set. Using default value of 1e-4.");
      m_F_tol = 1e-4;
    }
    m_msg->write_config("integrator.FIRE.F_tol",lexical_cast<string>(m_F_tol));
    if (param.find("e_tolerance") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator e_tolerance set to "+param["e_tolerance"]+".");
      m_E_tol = lexical_cast<double>(param["e_tolerance"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator e_tolerance not set. Using default value of 1e-6.");
      m_E_tol = 1e-6;
    }
    m_msg->write_config("integrator.FIRE.E_tol",lexical_cast<string>(m_E_tol));
    if (param.find("dt_inc") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator step size increase factor dt_inc set to "+param["dt_inc"]+".");
      m_f_inc = lexical_cast<double>(param["dt_inc"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator step size increase factor dt_inc not set. Using default value of 1.1.");
      m_f_inc = 1.1;
    }
    m_msg->write_config("integrator.FIRE.f_inc",lexical_cast<string>(m_f_inc));
    if (param.find("dt_dec") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator step size decrease factor dt_dec set to "+param["dt_dec"]+".");
      m_f_dec = lexical_cast<double>(param["dt_dec"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator step size decrease factor dt_dec not set. Using default value of 0.5.");
      m_f_dec = 0.5;
    }
    m_msg->write_config("integrator.FIRE.f_dec",lexical_cast<string>(m_f_dec));
    if (param.find("dt_max") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator maximum step size dt_max set to "+param["dt_max"]+".");
      m_dt_max = lexical_cast<double>(param["dt_max"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator maximum step size dt_max not set. Using default value of 1.0.");
      m_dt_max = 1.0;
    }
    m_msg->write_config("integrator.FIRE.dt_max",lexical_cast<string>(m_dt_max));
    if (param.find("min_neg_steps") != param.end())
    {
      m_msg->msg(Messenger::INFO,"FIRE integrator minimum number of steps with energy decrease before increasing step size min_neg_steps set to "+param["min_neg_steps"]+".");
      m_N_min = lexical_cast<int>(param["min_neg_steps"]);
    }
    else
    {
      m_msg->msg(Messenger::WARNING,"FIRE integrator minimum number of steps with energy decrease before increasing step size min_neg_steps not set. Using default value of 5.");
      m_N_min = 5;
    }
    m_msg->write_config("integrator.FIRE.N_min",lexical_cast<string>(m_N_min));
    
  }
  
  //! Propagate system for a time step
  void integrate();
  
private:

  double m_alpha;                                   //!< alpha factor in the FIRE minimizer
  double m_F_tol;                                   //!< Force tolerance for checking convergence 
  double m_E_tol;                                   //!< Energy tolerance for checking convergence 
  int    m_N_min;                                   //!< Number of steps energy change is negative before allowing m_alpha and m_dt to adapt
  double m_f_inc;                                   //!< Increase factor for m_dt
  double m_f_dec;                                   //!< Decrease factor for m_dt
  double m_alpha_init;                              //!< Initial value for m_alpha
  double m_f_alpha;                                 //!< Decrease alpha by this much 
  int    m_last_neg;                                //!< Last time (in terms of step number when force.velocity was negative
  double m_old_energy;                              //!< Old value of energy
  bool   m_converged;                               //!< Flag which tests if the method has converged
  double m_dt_max;                                  //!< Maximum time step

  
};

typedef shared_ptr<IntegratorFIRE> IntegratorFIREPtr;

#endif

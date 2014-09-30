/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_coulomb_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Declaration of PairCoulombPotential class
 */ 

#ifndef __PAIR_COULOMB_POTENTIAL_HPP__
#define __PAIR_COULOMB_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the Coulomb pair potential
struct CoulombParameters
{
  double alpha;
  double sigma;
};

/*! PairCoulombPotential implements standard Coulomb potential with a LJ repulsive part
 *  to avoid diverging forces and interaction strengths for collapsing particles. Potential is given as
 *  \f$ U_{Coul}\left(r_{ij}\right) = \frac{\alpha}{r_{ij}} + 4\left|\alpha\right|\left(\frac \sigma r_{ij}\right)^{12} \f$,
 *  where \f$ \alpha \f$ is the potential strength, \f$ \sigma \f$ is the particle diameter and \f$ r_{ij} \f$ is the 
 *  interparticle distance.
 */
class PairCoulombPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (alpha and sigma)
  PairCoulombPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairPotential(sys, msg, nlist, param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("alpha") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (alpha) specified for Coulomb pair potential. Setting it to 1.");
      m_alpha = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (alpha) for Coulomb pair potential is set to "+param["alpha"]+".");
      m_alpha = lexical_cast<double>(param["alpha"]);
    }
    
    if (param.find("sigma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No particle diameter (sigma) specified for Coulomb pair potential. Setting it to 1.");
      m_sigma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global particle diameter (sigma) for Coulomb pair potential is set to "+param["sigma"]+".");
      m_sigma = lexical_cast<double>(param["sigma"]);
    }
    m_pair_params = new CoulombParameters*[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_pair_params[i] = new CoulombParameters[ntypes];
      for (int j = 0; j < ntypes; j++)
      {
        m_pair_params[i][j].alpha = m_alpha;
        m_pair_params[i][j].sigma = m_sigma;
      }
    }
  }
  
  virtual ~PairCoulombPotential()
  {
    for (int i = 0; i < m_system->get_ntypes(); i++)
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in Coulomb potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in Coulomb potential.");
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
    if (pair_param.find("sigma") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Coulomb pair potential. Setting sigma to "+pair_param["sigma"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = lexical_cast<double>(pair_param["sigma"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Coulomb pair potential. Using default sigma ("+lexical_cast<string>(m_sigma)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = m_sigma;
    }
    
    m_pair_params[type_1-1][type_2-1].alpha = param["alpha"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].alpha = param["alpha"];
    m_pair_params[type_1-1][type_2-1].sigma = param["sigma"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].sigma = param["sigma"];
    
    m_has_pair_params = true;
  }
  
  //! Returns false since Coulomb potential does not need neighbour list
  bool need_nlist() { return false; }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_alpha;       //!< potential strength
  double m_sigma;       //!< particle diameter
  CoulombParameters** m_pair_params;   //!< type specific pair parameters 
    
};

typedef shared_ptr<PairCoulombPotential> PairCoulombPotentialPtr;

#endif

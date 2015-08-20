/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file bond_harmonic_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 05-Nov-2014
 * \brief Declaration of BondHarmonicPotential class
 */ 

#ifndef __BOND_HARMONIC_POTENTIAL_HPP__
#define __BOND_HARMONIC_POTENTIAL_HPP__

#include <cmath>

#include "bond_potential.hpp"

using std::make_pair;
using std::sqrt;

//! Structure that handles parameters for the harmonic bond
struct BondHarmoicParameters
{
  double k;    // Spring constant
  double l0;   // Native length
};

/*! BondHarmonicPotential implements standard harmonic bond potential 
 *  Potential is given as \f$ U_{harm}\left(r_{ij}\right) = \frac{1}{2}k\left(\left|\vec r_{ij}\right|-l_0\right)^2 \f$,
 *  where \f$ k \f$ is the spring constant, \f$ l_0 \f$ is the the rest length and \f$ \vec r_{ij} = \vec r_i - \vec r_j \f$.
 */
class BondHarmonicPotential : public BondPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param param Contains information about all parameters (k and l0)
  BondHarmonicPotential(SystemPtr sys, MessengerPtr msg,  pairs_type& param) : BondPotential(sys, msg, param)
  {
    int ntypes = m_system->get_n_bond_types();
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No spring constant (k) specified for harmonic bond potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global spring constant (k) for harmonic bond potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    m_msg->write_config("potential.bond.harmonic.k",lexical_cast<string>(m_k));
    if (param.find("l_eq") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No rest length (l0) specified for harmonic bond potential. Setting it to 1.");
      m_l0 = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global rest length (l0) for harmonic bond potential is set to "+param["l_eq"]+".");
      m_l0 = lexical_cast<double>(param["l_eq"]);
    }
    m_msg->write_config("potential.bond.harmonic.l_eq",lexical_cast<string>(m_l0));
    m_bond_params = new BondHarmoicParameters[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_bond_params[i].k = m_k;
      m_bond_params[i].l0 = m_l0;
    }
  }
  
  virtual ~BondHarmonicPotential()
  {
    delete [] m_bond_params;
  }
                                                                                                                
  //! Set parameters data for specific bond type    
  void set_bond_parameters(pairs_type& bond_param)
  {
    map<string,double> param;
    
    int type;
    
    if (bond_param.find("type") == bond_param.end())
    {
      m_msg->msg(Messenger::ERROR,"Bond type has not been defined for harmonic bond potential.");
      throw runtime_error("Missing key for harmonic bond parameters.");
    }
    type = lexical_cast<int>(bond_param["type"]);
        
    if (bond_param.find("k") != bond_param.end())
    {
      m_msg->msg(Messenger::INFO,"Harmonic bond potential. Setting spring constant to "+bond_param["k"]+" for bond of types "+lexical_cast<string>(type)+".");
      param["k"] = lexical_cast<double>(bond_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Harmonic bond potential. Using default spring constant ("+lexical_cast<string>(m_k)+") for bonds of types "+lexical_cast<string>(type)+".");
      param["k"] = m_k;
    }
    if (bond_param.find("l_eq") != bond_param.end())
    {
      m_msg->msg(Messenger::INFO,"Harmonic bond potential. Setting rest length to "+bond_param["l_eq"]+" for bond of types "+lexical_cast<string>(type)+".");
      param["l_eq"] = lexical_cast<double>(bond_param["l_eq"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Harmonic bond potential. Using default rest length ("+lexical_cast<string>(m_l0)+") for bonds of types "+lexical_cast<string>(type)+".");
      param["l_eq"] = m_l0;
    }
    m_msg->write_config("potential.bond.harmonic.type_"+bond_param["type"]+".k",lexical_cast<string>(param["k"]));
    m_msg->write_config("potential.bond.harmonic.type_"+bond_param["type"]+".l_eq",lexical_cast<string>(param["l_eq"]));
    
    m_bond_params[type-1].k = param["k"];
    m_bond_params[type-1].l0 = param["l_eq"];
         
    m_has_bond_params = true;
  }
  
   
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
        
  double m_k;       //!< spring constant
  double m_l0;      //!< rest length
  BondHarmoicParameters* m_bond_params;   //!< type specific bond parameters 
    
};

typedef shared_ptr<BondHarmonicPotential> BondHarmonicPotentialPtr;

#endif

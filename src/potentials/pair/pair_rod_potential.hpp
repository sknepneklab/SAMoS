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
 *   (c) 2013, 2014, 2015
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file pair_rod_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 23-Mar-2015
 * \brief Declaration of PairRodPotential class
 */ 

#ifndef __PAIR_ROD_POTENTIAL_HPP__
#define __PAIR_ROD_POTENTIAL_HPP__

#include <cmath>
#include <string>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;
using std::fabs;
using std::string;

//! Structure that handles parameters for the rod pair potential
struct RodParameters
{
  double k;    // elastic constant
};

/*! PairRodPotential implements the anisotropic soft repulsion between rods of length \f$ l \f$ 
 *  and radius \f$ \sigma \f$. Details of the interaction potential, force and torques are straightforward
 *  to derive, but lengthy. They are given in a separate set of notes in the docs directory. 
 * 
 *  \note Rod length is encoded within the potential and not as a part of the Particle class as
 *  for the time being we would like to keep the core of the code work with spherical objects only. If the 
 *  anisotropy turns out to be important, this may change and the Particle class might start carrying
 *  quaternion describing its orientation. This would allow of a uniform interface for implementing 
 *  arbitrary anisotropic potentials. 
 * 
 *  \todo Add brief description of the model here.
 */
class PairRodPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param val Value control object (for phasing in)
  //! \param param Contains information about all parameters (k)
  PairRodPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, ValuePtr val, pairs_type& param) : PairPotential(sys, msg, nlist, val, param)
  {
    int ntypes = m_system->get_ntypes();
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (k) specified for rod pair potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (k) for rod pair potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
    if (param.find("phase_in") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod pair potential. Gradually phasing in the potential for new particles.");
      m_phase_in = true;
    }    
    if (param.find("model") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No contact model specified for rod pair potential. Assuming soft potential.");
      m_model = "soft";
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Contact model for rod pair potential is set to "+param["model"]+".");
      m_model = param["model"];
    }
    
    m_pair_params = new RodParameters*[ntypes];
    for (int i = 0; i < ntypes; i++)
    {
      m_pair_params[i] = new RodParameters[ntypes];
      for (int j = 0; j < ntypes; j++)
        m_pair_params[i][j].k = m_k;
    }
  }
  
  virtual ~PairRodPotential()
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
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in rod potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in rod potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("k") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Rod pair potential. Setting strength to "+pair_param["k"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = lexical_cast<double>(pair_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Rod pair potential. Using default strength ("+lexical_cast<string>(m_k)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = m_k;
    }
           
    m_pair_params[type_1-1][type_2-1].k = param["k"];
    if (type_1 != type_2)
      m_pair_params[type_2-1][type_1-1].k = param["k"];
        
    m_has_pair_params = true;
  }
  
  //! Returns true since rod potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute(double);
  
  
private:
       
  double m_k;                      //!< potential strength
  string m_model;                  //!< soft potential or Hertzian model
  RodParameters** m_pair_params;   //!< type specific pair parameters 
     
};

typedef shared_ptr<PairRodPotential> PairRodPotentialPtr;

#endif

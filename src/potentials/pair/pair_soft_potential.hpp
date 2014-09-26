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
 * \file pair_soft_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 31-Oct-2013
 * \brief Declaration of PairSoftPotential class
 */ 

#ifndef __PAIR_SOFT_POTENTIAL_HPP__
#define __PAIR_SOFT_POTENTIAL_HPP__

#include <cmath>

#include "pair_potential.hpp"

using std::make_pair;
using std::sqrt;

/*! PairSoftPotential implements the soft repulsive force between particles.
 *  Potential is given as \f$ U_{soft}\left(r_{ij}\right) = \frac{k}{2} \left(a_i + a_j - r_{ij}\right)^2 \f$, if
 *  \f$ r_{ij} \le a_i + a_j \f$ or \f$ U_{soft} = 0 \f$ if \f$ r_{ij} > a_i + a_j \f$.
 *  \f$ k \f$ is the potential strength, \f$ a_i \f$ and \f$ a_j \f$ are radii of the two particles and \f$ r_{ij} \f$ is the 
 *  interparticle distance.
 */
class PairSoftPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (k)
  PairSoftPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairPotential(sys, msg, nlist, param)
  {
    if (param.find("k") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential strength (k) specified for soft pair potential. Setting it to 1.");
      m_k = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential strength (k) for soft pair potential is set to "+param["k"]+".");
      m_k = lexical_cast<double>(param["k"]);
    }
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(pairs_type& pair_param)
  {
    map<string,double> param;
    int type_1, type_2;
    
    if (pair_param.find("type_1") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_1 has not been defined for pair potential parameters in soft potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    if (pair_param.find("type_2") == pair_param.end())
    {
      m_msg->msg(Messenger::ERROR,"type_2 has not been defined for pair potential parameters in soft potential.");
      throw runtime_error("Missing key for pair potential parameters.");
    }
    type_1 = lexical_cast<int>(pair_param["type_1"]);
    type_2 = lexical_cast<int>(pair_param["type_2"]);
    
    if (pair_param.find("k") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Soft pair potential. Setting strength to "+pair_param["k"]+" for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = lexical_cast<double>(pair_param["k"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Soft pair potential. Using default strength ("+lexical_cast<string>(m_k)+") for particle pair of types ("+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["k"] = m_k;
    }
        
    m_pair_params[make_pair(type_1,type_2)]["k"] = param["k"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["k"] = param["k"];
    
    m_has_pair_params = true;
  }
  
  //! Returns true since soft potential needs neighbour list
  bool need_nlist() { return true; }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_k;       //!< potential strength
     
};

typedef shared_ptr<PairSoftPotential> PairSoftPotentialPtr;

#endif

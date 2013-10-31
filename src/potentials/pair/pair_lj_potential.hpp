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
 * \file pair_lj_potential.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 22-Oct-2013
 * \brief Declaration of PairLJPotential class
 */ 

#ifndef __PAIR_LJ_POTENTIAL_HPP__
#define __PAIR_LJ_POTENTIAL_HPP__

#include "pair_potential.hpp"

using std::make_pair;

/*! PairLJPotential implements standard Lennard Jones potential 
 * \f$ U_{LJ}\left(r_{ij}\right) = 4\varepsilon \left[\left(\frac \sigma r_{ij}\right)^{12}-\left(\frac \sigma r_{ij}\right)^6)\right] \f$,
 *  where \f$ \varepsilon \f$ is the potential strength, \f$ \sigma \f$ is the particle diameter and \f$ r_{ij} \f$ is the 
 *  interparticle distance.
 */
class PairLJPotential : public PairPotential
{
public:
  
  //! Constructor
  //! \param sys Pointer to the System object
  //! \param msg Pointer to the internal state messenger
  //! \param nlist Pointer to the global neighbour list
  //! \param param Contains information about all parameters (epsilon, sigma, and rcut)
  PairLJPotential(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, pairs_type& param) : PairPotential(sys, msg, nlist, param)
  {
    if (!m_nlist)
    {
      m_msg->msg(Messenger::ERROR,"Lennard Jones pair potential requires neighbour list. None given.");
      throw runtime_error("Neighbour list required by Lennard Jones potential, but non specified.");
    }
    if (param.find("epsilon") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No potential depth (epsilon) specified for the Lennard Jones pair potential. Setting it to 1.");
      m_eps = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global potential depth (epsilon) for Lennard Jones pair potential is set to "+param["epsilon"]+".");
      m_eps = lexical_cast<double>(param["epsilon"]);
    }
    
    if (param.find("sigma") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No particle diameter (sigma) specified for the Lennard Jones pair potential. Setting it to 1.");
      m_sigma = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global particle diameter (sigma) for Lennard Jones pair potential is set to "+param["sigma"]+".");
      m_sigma = lexical_cast<double>(param["sigma"]);
    }
    
    if (m_param.find("rcut") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance (rcut) specified for the Lennard Jones pair potential. Setting it to 3.0.");
      m_rcut = 3.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Global cutoff distance (rcut) for Lennard Jones pair potential is set to "+param["rcut"]+".");
      m_rcut = lexical_cast<double>(param["rcut"]);
    }
    
    if (m_rcut < m_nlist.get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<double>(m_nlist.get_cutoff())+
      "is smaller than the Lennard-Jones cuttof distance ("+lexical_cast<double>(m_rcut)+
      "). Results will not be reliable.");
    
    if (param.find("shifted") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Lennard-Jones potential shifted to zero at cutoff.");
      m_shifted = true;
    }
  }
                                                                                                                
  //! Set pair parameters data for pairwise interactions    
  void set_pair_parameters(int type_1, int type_2, pairs_type& pair_param)
  {
    map<string,double> param;
    
    if (pair_param.find("epsilon") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Setting epsilon to "+pair_param["epsilon"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = lexical_cast<double>(pair_param["epsilon"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Using default epsilon ("+lexical_cast<string>(m_eps)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["epsilon"] = m_eps;
    }
    if (pair_param.find("sigma") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Setting sigma to "+pair_param["sigma"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = lexical_cast<double>(pair_param["sigma"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Using default sigma ("+lexical_cast<string>(m_sigma)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["sigma"] = m_sigma;
    }
    if (pair_param.find("rcut") != pair_param.end())
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Setting rcut to "+pair_param["rcut"]+" for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = lexical_cast<double>(pair_param["rcut"]);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Lennard Jones pair potential. Using default rcut ("+lexical_cast<string>(m_rcut)+") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+").");
      param["rcut"] = m_sigma;
    }
    
    m_pair_params[make_pair(type_1,type_2)]["epsilon"] = param["epsilon"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["epsilon"] = param["epsilon"];
    m_pair_params[make_pair(type_1,type_2)]["sigma"] = param["sigma"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["sigma"] = param["sigma"];
    m_pair_params[make_pair(type_1,type_2)]["rcut"] = param["rcut"];
    if (type_1 != type_2)
      m_pair_params[make_pair(type_2,type_1)]["rcut"] = param["rcut"];
    
    if (param["rcut"] < m_nlist.get_cutoff())
      m_msg->msg(Messenger::WARNING,"Neighbour list cutoff distance (" + lexical_cast<double>(m_nlist.get_cutoff())+
      "is smaller than the Lennard-Jones cuttof distance ("+lexical_cast<double>(m_rcut)+
      ") for particle pair of types "+lexical_cast<string>(type_1)+" and "+lexical_cast<string>(type_2)+"). Results will not be reliable.");
    
    m_has_pair_params = true;
  }
  
  //! Computes potentials and forces for all particles
  void compute();
  
  
private:
       
  double m_eps;    //!< potential depth 
  double m_sigma;  //!< particle diameter
  double m_rcut;   //!< cutoff distance
    
};


#endif

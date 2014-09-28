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
 *   (c) 2013, 2014
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file population_random.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Declaration of PopulationRandom class.
 */ 

#ifndef __POPULATION_RANDOM_HPP__
#define __POPULATION_RANDOM_HPP__

#include <list>
#include <stdexcept>

#include "population.hpp"
#include "rng.hpp"


using std::list;
using std::runtime_error;

/*! PopulationRandom class handles very simple particle division and death.
 *  We set division rate and death rate. Particles are divided and removed
 *  based on the their age, by a random event.
 * 
 *  Particles are always divided along the direction of the orientation vector 
 *  n.
 *
*/
class PopulationRandom : public Population
{
public:
  
  //! Construct PopulationRandom object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  PopulationRandom(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : Population(sys,msg,param)
  { 
    if (param.find("seed") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No random number generator seed specified. Using default 0.");
      m_rng = make_shared<RNG>(0);
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting random number generator seed to "+param["seed"]+".");
      m_rng = make_shared<RNG>(lexical_cast<int>(param["seed"]));
    }
    if (param.find("division_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No division rate set. Using default 1.");
      m_div_rate = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting division rate "+param["division_rate"]+".");
      m_div_rate =  lexical_cast<double>(param["division_rate"]);
    }
    if (param.find("death_rate") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Random population control. No death rate set. Using default 1.");
      m_death_rate = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Random population control. Setting death rate "+param["death_rate"]+".");
      m_death_rate =  lexical_cast<double>(param["death_rate"]);
    }
  }
  
  //! Particle division (emulates cell division)
  void divide(int);
  
  //! Remove particle (emulates cell death)
  void remove(int);
  
  //! Add particle (has no direct biological application)
  void add(int t) { }
  
  
  
private:
  
  RNGPtr  m_rng;          //!< Random number generator
  double  m_div_rate;     //!< Rate of division
  double  m_death_rate;   //!< Rate of death
   
};

typedef shared_ptr<PopulationRandom> PopulationRandomPtr;

#endif
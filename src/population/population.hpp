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
 * \file population.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 28-Sept-2014
 * \brief Declaration of Population class.
 */ 

#ifndef __POPULATION_HPP__
#define __POPULATION_HPP__

#include <string>
#include <vector>

#include <boost/make_shared.hpp>

#include "messenger.hpp"
#include "system.hpp"

#include "parse_parameters.hpp"

using std::string;
using std::vector;

using boost::make_shared;

/*! Population class is the abstract parent class that handles all population control 
 *  mechanism (division and death/removal of particles).
 *
*/
class Population
{
public:
  
  //! Construct Population object
  //! \param sys Reference to the System object
  //! \param msg Reference to the system wide messenger
  //! \param param Reference to the parameters that control given population control
  Population(SystemPtr sys, const MessengerPtr msg, pairs_type& param) : m_system(sys), m_msg(msg)
  { 
    if (param.find("group") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Population. No group has been set. Assuming group 'all'.");
      m_group_name = "all";
    }
    else
    {
      m_group_name = param["group"];
      m_msg->msg(Messenger::INFO,"Population. Applying population to group "+m_group_name+".");
    }
    if (param.find("freq") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Population. No population control frequency has been set. Assuming default 100.");
      m_freq = 100;
    }
    else
    {
      m_freq = lexical_cast<int>(param["freq"]);
      m_msg->msg(Messenger::INFO,"Population. Population control frequency set to "+param["freq"]+".");     
    }
    m_msg->write_config("population.freq",lexical_cast<string>(m_freq));
  }
  
  //! Particle division (emulates cell division)
  virtual void divide(int) = 0;
  
  //! Remove particle (emulates cell death)
  virtual void remove(int) = 0;
  
  //! Add particle (has no direct biological application)
  virtual void add(int) = 0;
  
  //! Change particle radius
  virtual void grow(int) = 0;
  
  //! Change particle length
  virtual void elongate(int) = 0;
  
protected:
  
  SystemPtr m_system;            //!< Contains pointer to the System object
  MessengerPtr m_msg;            //!< Handles messages sent to output
  string m_group_name;           //!< Name of the group to apply this population to
  int m_freq;                    //!< Frequency (in the number of time steps) with which to attempt population control
  double m_l_rescale;            //!< Rescale particle length by total of this amount
  int m_l_rescale_steps;         //!< Rescale particle lenght over this many steps
  double m_l_scale;              //!< Rescale particle lenth by this much in each step (=m_l_rescale**(m_grow_l_freq/m_l_rescale_steps))
   
};

typedef shared_ptr<Population> PopulationPtr;

#endif
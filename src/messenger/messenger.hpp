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
 * \file messenger.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Declaration of Messenger class.
 */ 

#ifndef __MESSENGER_HPP__
#define __MESSENGER_HPP__

#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/shared_ptr.hpp>

using std::ostream;
using std::cout;
using std::ofstream;
using std::string;
using std::map;
using std::endl;
using std::ios_base;

using boost::to_lower_copy;
using namespace boost::posix_time;
using namespace boost::gregorian;
using boost::shared_ptr;

/*! Messenger class handles all sort of messages that can either be sent to 
 *  the terminal or to a file.
*/

class Messenger
{
public:
  
  //! Message types
  enum MSG_TYPE
  {
    INFO,      /*!< Basic information */
    WARNING,   /*!< Warning */
    ERROR      /*!< Error message */
  };
  
  //! Construct Messenger object
  Messenger(string&);
  
  //! Destructor
  ~Messenger();
  
  //! Output message
  void msg(const Messenger::MSG_TYPE&, const string&) const;
  
private:
  
  ostream& m_out;       //!< Output stream to which to send output
  string m_file_name;   //!< Messenger file name
  map<const Messenger::MSG_TYPE,string> m_msg_type;  //!< maps message types to appropriate strings
  
};

typedef shared_ptr<Messenger> MessengerPtr;

#endif
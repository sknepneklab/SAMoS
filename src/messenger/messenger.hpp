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
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/xml_parser.hpp>

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
namespace pt = boost::property_tree;

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
  Messenger(const string&);
  
  //! Destructor
  ~Messenger();
  
  //! Add configuration file handling
  //! \param filename configuration file
  //! \param type type of the configuration file (XML or JSON)
  void add_config(const string& filename, const string& type)
  {
    m_config_type = type;
    m_config_file = filename;
    m_has_config = true;
  }
  
  //! Writes a note to the configuration file
  //! \param key name of the node to write
  //! \param val value to write
  void write_config(const string& key, const string& val)
  {
    if (m_has_config)
      m_config.put(key,val);
  }
  
  //! Output message
  void msg(const Messenger::MSG_TYPE&, const string&); //const;
  
private:
  
  ofstream m_out;      //!< Output stream to which to send output
  string m_file_name;   //!< Messenger file name
  map<const Messenger::MSG_TYPE,string> m_msg_type;  //!< maps message types to appropriate strings
  bool m_to_terminal;   //!< If true send output to terminal
  bool m_has_config;    //!< If true store configuration parameters
  string m_config_type; //!< Type of the file to store configuration into
  string m_config_file; //!< File name to store configuration into
  pt::ptree m_config;       //!< Boost property tree storing configuration parameters
  
};

typedef shared_ptr<Messenger> MessengerPtr;

#endif
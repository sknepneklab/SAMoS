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
 * \file dump.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Declaration of Dump class
 */ 

#ifndef __DUMP_H__
#define __DUMP_H__

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <list>
#include <boost/format.hpp>

#include "system.hpp"
#include "parse_parameters.hpp"

using std::string;
using std::ofstream;
using std::map;
using std::endl;
using std::list;
using boost::format;

/*! Dump class handles output of system's state, such 
 *  as particle coordinates, velocities, forces, etc.
 *  It supports number of formats that can be used for
 *  subsequent data analysis and visualization.
 */
class Dump
{
public:
  //! Constructor
  Dump(SystemPtr, MessengerPtr, const string&, pairs_type&);
  
  //! Destructor
  ~Dump()
  {
    if (m_multi_print) m_out.close();
    m_type_ext.clear();
    m_to_print.clear();
  }
  
  //! Do actual dump 
  void dump(int);
  
private:
  
  SystemPtr m_system;           //!< Pointer to the System object
  MessengerPtr m_msg;           //!< Handles system wide messages
  string m_file_name;           //!< Base file name for the output file
  pairs_type& m_params;         //!< Control parameters
  string m_ext;                 //!< File name extension
  ofstream m_out;               //!< Output file stream with the file
  int m_start;                  //!< First time step of the dump
  int m_freq;                   //!< Frequency with which to dump data
  bool m_multi_print;           //!< Print to multiple files
  string m_type;                //!< File format for the dump 
  bool m_no_header;             //!< Do not print header in some files
  
  // Auxiliary data structures
  map<string, string> m_type_ext;  //!< Hold extension for a given data type
  vector<string> m_to_print;       //!< List of quantities to dump (e.g., particle type, id, coordinate, velocity, etc.)
  
  
  // private member methods that do actual dumping
  // these methods cannot be called directly 
  //! DCD dump
  void dump_dcd();         
  //! XYZ dump (XYZ file format)
  void dump_xyz();
  //! Dump a selected list of parameters
  void dump_data();
  //! Dump input format for restarts
  void dump_input();
  //! Dump velocities
  void dump_velocity();
  
};

typedef shared_ptr<Dump> DumpPtr;

#endif
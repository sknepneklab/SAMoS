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
#include <boost/algorithm/string.hpp>

// Handles data compression (requires zlib)
#include <boost/iostreams/filtering_stream.hpp>    
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bo = boost::iostreams;

#include "system.hpp"
#include "parse_parameters.hpp"
#include "neighbour_list.hpp"

using std::string;
using std::ofstream;
using std::map;
using std::endl;
using std::list;
using boost::format;
using boost::replace_all;

/*! Dump class handles output of system's state, such 
 *  as particle coordinates, velocities, forces, etc.
 *  It supports number of formats that can be used for
 *  subsequent data analysis and visualization.
 */
class Dump
{
public:
  //! Constructor
  Dump(SystemPtr, MessengerPtr, NeighbourListPtr, const string&, pairs_type&);
  
  //! Destructor
  ~Dump()
  {
    if (!m_multi_print) m_file.close();
    m_type_ext.clear();
    m_to_print.clear();
  }
  
  //! Do actual dump 
  void dump(int);
  
private:
  
  SystemPtr m_system;           //!< Pointer to the System object
  MessengerPtr m_msg;           //!< Handles system wide messages
  NeighbourListPtr m_nlist;     //!< Pointer to the global neighbour list
  string m_file_name;           //!< Base file name for the output file
  pairs_type m_params;          //!< Control parameters
  string m_ext;                 //!< File name extension
  ofstream m_file;              //!< Output file stream with the file
  bo::filtering_ostream m_out;  //!< Filters data through a compressor
  int m_start;                  //!< First time step of the dump
  int m_freq;                   //!< Frequency with which to dump data
  bool m_multi_print;           //!< Print to multiple files
  bool m_print_header;          //!< If true, print header line on some file types
  string m_type;                //!< File format for the dump 
  bool m_no_header;             //!< Do not print header in some files
  double m_r_cut;               //!< Cutoff distance for the contact network
  bool m_compress;              //!< Compress output data using gzip compression
  
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
  //! Dump director 
  void dump_director();
  //! Dump XYZV file format for visualization with SimRePlay
  void dump_xyzv();
  //! Dump XYZC file format for visualization with SimRePlay
  void dump_xyzc();
  //! Dump MOL2 file format for visualization with VMD
  void dump_mol2();
  //! Dump contact network for further analysis
  void dump_contact();
  //! Dump faces based on the contact network
  void dump_faces();
  
};

typedef shared_ptr<Dump> DumpPtr;

#endif
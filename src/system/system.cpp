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
 * \file system.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of System class.
 */ 

#include "system.hpp"

// Set of auxiliary functions that help parse lines of the input file and
// reduce code bloating of the System class constructor
// These are defined as static 

//! Splits a string at one of the characters " ", "\t", or ","
//! \param str string to split
//! \return vector of split parts
static vector<string> split_line(const string& str)
{
  vector<string> split_string;
  split_regex(split_string, str, regex("\\s+|,|\\t"));
  return split_string;  
}


/*! System constructor
 *  The constructor reads in coordinates from a file.
 *  If there is an error, an exception will be thrown.
 * 
 *  File format:
 *   1. Any line starting with # sign is a comment and is ignored
 *  2. Otherwise, information about each particle is stored in a one 
 *  3. Particle ids should be stored in the incremental order
 *  4. Particle line
 *    [id] [type] [radius] [x] [y] [z] [vx] [vy] [vz] [omega]
 *  
 *
 * \param input_filename name of coordinate file
 * \param msg Pointer to the Messenger object
 * \param box Pointer to the simulation box object
 */
System::System(const string& input_filename, MessengerPtr msg, BoxPtr box) : m_msg(msg), m_box(box), m_periodic(false)
{
  vector<string> s_line;
  ifstream inp;
  inp.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  try
  {
    inp.open(input_filename.c_str());
  }
  catch (exception& e)
  {
    msg->msg(Messenger::ERROR,"Problem opening file "+input_filename);
    throw e;
  }
  
  string line;
  msg->msg(Messenger::INFO,"Reading particle coordinates from file: "+input_filename);
  
  inp.exceptions ( std::ifstream::badbit ); // need to reset ios exceptions to avoid EOF failure of getline
  while ( getline(inp, line) )
  {
    trim(line);
    to_lower(line);
    if (line[0] != '#' && line.size() > 0)
    {
      s_line = split_line(line);
      if (s_line.size() < NUM_PART_ATTRIB)
      {
        msg->msg(Messenger::ERROR,"Insufficient number of parameters to define a particle. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_PART_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the input file.");
      }
      const int id = lexical_cast<int>(s_line[0]);  // read particle id
      const int tp = lexical_cast<int>(s_line[1]);  // read particle type
      const double r = lexical_cast<double>(s_line[2]);  // read particle radius
      Particle p(id, tp, r);
      p.x = lexical_cast<double>(s_line[3]);
      if (p.x < m_box->xlo || p.x > m_box->xhi)
      {
        m_msg->msg(Messenger::ERROR,"X coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
        throw runtime_error("Particle outside the box.");
      }
      p.y = lexical_cast<double>(s_line[4]);
      if (p.y < m_box->ylo || p.y > m_box->yhi)
      {
        m_msg->msg(Messenger::ERROR,"Y coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
        throw runtime_error("Particle outside the box.");
      }
      p.z = lexical_cast<double>(s_line[5]);
      if (p.z < m_box->zlo || p.z > m_box->zhi)
      {
        m_msg->msg(Messenger::ERROR,"Z coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
        throw runtime_error("Particle outside the box.");
      }
      p.vx = lexical_cast<double>(s_line[6]);
      p.vy = lexical_cast<double>(s_line[7]);
      p.vz = lexical_cast<double>(s_line[8]);
      p.nx = lexical_cast<double>(s_line[9]);
      p.ny = lexical_cast<double>(s_line[10]);
      p.nz = lexical_cast<double>(s_line[11]);
      p.omega = lexical_cast<double>(s_line[12]);
      m_particles.push_back(p);
    }
  }
  msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_particles.size())+" particles.");
  inp.close();
}
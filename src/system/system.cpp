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
  
  // Populate group 'all'
  m_group["all"] = make_shared<Group>(Group(0,"all"));
  for (unsigned int i = 0; i < m_particles.size(); i++)
    m_group["all"]->add_particle(i);
  m_num_groups = 1;
  
  msg->msg(Messenger::INFO,"Generated group 'all' containing all particles.");
  
  this->disable_per_particle_eng();
}

/*! Generate a group of particles
    \param name name of the group
    \param param Contains information about all parameters
    
    Groups can be generated based on a number of criteria
    For example base on the particle type, particle size (radius), position, etc.
    
*/
void System::make_group(const string name, pairs_type& param)
{
  m_group[name] = make_shared<Group>(Group(m_num_groups++, name));
  map<string, vector<bool> > to_add;
  if (param.find("type") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["type"].push_back(false);
    m_msg->msg(Messenger::INFO,"Adding particle by group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_type() == lexical_cast<int>(param["type"]))
        to_add["type"][i] = true;
  }
  if (param.find("r_eq") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["r_eq"].push_back(false);
    m_msg->msg(Messenger::INFO,"Adding particle of radius equal to "+param["r_eq"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_radius() == lexical_cast<double>(param["r_eq"]))
        to_add["r_eq"][i] = true;
  }
  if (param.find("r_le") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["r_le"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle of radius less or equal "+param["r_le"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_radius() <= lexical_cast<double>(param["r_le"]))
        to_add["r_le"][i] = true;
  }
  if (param.find("r_ge") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["r_ge"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle of radius greater or equal to "+param["r_ge"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_radius() >= lexical_cast<double>(param["r_ge"]))
        to_add["r_ge"][i] = true;
  }
  if (param.find("x_le") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["x_le"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with x coordinate less or equal to "+param["x_le"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].x <= lexical_cast<double>(param["x_le"]))
        to_add["x_le"][i] = true;
  }
  if (param.find("x_ge") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["x_ge"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with x coordinate greater or equal to "+param["x_ge"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].x >= lexical_cast<double>(param["x_ge"]))
        to_add["x_ge"][i] = true;
  }
  if (param.find("y_le") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["y_le"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with y coordinate less or equal to "+param["y_le"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].y <= lexical_cast<double>(param["y_le"]))
        to_add["y_le"][i] = true;
  }
  if (param.find("y_ge") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["y_ge"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with y coordinate greater or equal to "+param["y_ge"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].y >= lexical_cast<double>(param["y_ge"]))
        to_add["y_ge"][i] = true;
  }
  if (param.find("z_le") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["z_le"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with z coordinate less or equal to "+param["z_le"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].z <= lexical_cast<double>(param["z_le"]))
        to_add["z_le"][i] = true;
  }
  if (param.find("z_ge") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["z_ge"].push_back(false);
          m_msg->msg(Messenger::INFO,"Adding particle with z coordinate greater or equal to "+param["z_ge"]+" to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].z >= lexical_cast<double>(param["z_ge"]))
        to_add["z_ge"][i] = true;
  }
  
  for (unsigned int i = 0; i < m_particles.size(); i++)
  {
    bool add = true;
    for (map<string, vector<bool> >::iterator it = to_add.begin(); it != to_add.end(); it++)
      add = add && (*it).second[i];
    if (add)
      m_group[name]->add_particle(i);
  }
  
  m_msg->msg(Messenger::INFO,"Added "+lexical_cast<string>(m_group[name]->get_size())+" particles to group "+name+".");
      
}
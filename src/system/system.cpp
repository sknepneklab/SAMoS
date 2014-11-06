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
System::System(const string& input_filename, MessengerPtr msg, BoxPtr box) : m_msg(msg), m_box(box), m_periodic(false), m_force_nlist_rebuild(false), m_current_particle_flag(0)
{
  vector<int> types;
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
      if (find(types.begin(), types.end(), tp) == types.end())
        types.push_back(tp);
      p.vx = lexical_cast<double>(s_line[6]);
      p.vy = lexical_cast<double>(s_line[7]);
      p.vz = lexical_cast<double>(s_line[8]);
      p.nx = lexical_cast<double>(s_line[9]);
      p.ny = lexical_cast<double>(s_line[10]);
      p.nz = lexical_cast<double>(s_line[11]);
      p.omega = lexical_cast<double>(s_line[12]);
      p.set_flag(m_current_particle_flag);
      m_current_particle_flag++;
      m_particles.push_back(p);
    }
  }
  msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_particles.size())+" particles.");
  inp.close();
  
  // Populate group 'all'
  m_group["all"] = make_shared<Group>(Group(0,"all"));
  for (unsigned int i = 0; i < m_particles.size(); i++)
  {
    m_group["all"]->add_particle(i);
    m_particles[i].groups.push_back("all");
  }
  m_num_groups = 1;
  
  msg->msg(Messenger::INFO,"Generated group 'all' containing all particles.");
  
  m_n_types = types.size();
   
  msg->msg(Messenger::INFO,"There are " + lexical_cast<string>(m_n_types) + " distinct particle types in the system.");
  
  
  this->disable_per_particle_eng();
  
  m_has_exclusions = false;
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
    {
      Particle& p = m_particles[i];
      m_group[name]->add_particle(i);
      p.groups.push_back(name);
    }
  }
  
  m_msg->msg(Messenger::INFO,"Added "+lexical_cast<string>(m_group[name]->get_size())+" particles to group "+name+".");
      
}

/*! Add particle to the system
 *  \param p Particle to add
 */ 
void System::add_particle(Particle& p)
{
  p.set_flag(m_current_particle_flag);
  m_particles.push_back(p);
  for (list<string>::iterator it = p.groups.begin(); it != p.groups.end(); it++)
    m_group[*it]->add_particle(p.get_id());
  // We need to force neighbour list rebuild
  m_force_nlist_rebuild = true;
  m_current_particle_flag++;
}

/*! Remove particle from the system
 *  \param p Particle to remove
 */ 
void System::remove_particle(Particle& p)
{
  int id = p.get_id();
  m_particles.erase(m_particles.begin() + id);
  // Shift down all ids
  for (unsigned int i = 0; i < m_particles.size(); i++)
  {
    Particle& p = m_particles[i];
    if (p.get_id() > id)
      p.set_id(p.get_id() - 1);
  }
  // Update all groups
  for(map<string, GroupPtr>::iterator it_g = m_group.begin(); it_g != m_group.end(); it_g++)
    (*it_g).second->shift(id);
  // We need to force neighbour list rebuild
  m_force_nlist_rebuild = true;
}

/*! Change group of the particle
 *  \param p Particle whose group we want to change
 *  \param old_group Old group
 *  \param new_group New group 
*/
void System::change_group(Particle& p, const string& old_group, const string& new_group)
{
  if (old_group == "all" && new_group != "all")
  {
    m_msg->msg(Messenger::ERROR,"Particle"+lexical_cast<string>(p.get_id())+": Cannot change from \"all\" to other group ("+new_group+").");
    throw runtime_error("Can't change out of \"all\" group.");
  }
  list<string>::iterator it_g = find(p.groups.begin(),p.groups.end(),old_group);
  if (it_g == p.groups.end())
  {
    m_msg->msg(Messenger::ERROR,"Particle"+lexical_cast<string>(p.get_id())+" does not belong to group "+old_group+".");
    throw runtime_error("Unknown particle group.");
  }
  else
    p.groups.erase(it_g);
  
  if (new_group != "all") // It's already in "all"
    p.groups.push_back(new_group);
  
  m_group[old_group]->remove_particle(p.get_id());
  if (new_group != "all") // It's already in "all"
    m_group[new_group]->add_particle(p.get_id());
}

/*! Kill off any non-zero momentum and angular momentum that a group of particles might have pick up
 * if the 3d-Newton's law breaking m_limit was enabled.
 * \param group Group name
*/
void System::zero_cm_momentum(const string& group)
{
  if (m_group.find(group) == m_group.end())
  {
    m_msg->msg(Messenger::ERROR,"Group"+group+" does not exist, so its momentum cannot be zeroed out.");
    throw runtime_error("Unknown group.");
  }
  
  int N = m_group[group]->get_size();
  vector<int> particles = m_group[group]->get_particles();
  double vcm_x = 0.0, vcm_y = 0.0, vcm_z = 0.0;
  double tau_cm_x = 0.0, tau_cm_y = 0.0, tau_cm_z = 0.0;
  N = this->size();
  
  m_msg->msg(Messenger::INFO,"Removing centre of mass momentum and angular momentum of the group "+group+".");
  
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_particles[pi];
    vcm_x += p.vx;  vcm_y += p.vy;   vcm_z += p.vz;
    tau_cm_x += p.tau_x;  tau_cm_y += p.tau_y;   tau_cm_z += p.tau_z;
  }
  vcm_x /= N; vcm_y /= N; vcm_z /= N;
  tau_cm_x /= N; tau_cm_y /= N; tau_cm_z /= N;
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_particles[pi];
    p.vx -= vcm_x;  p.vy -= vcm_y;  p.vz -= vcm_z;
    p.tau_x -= tau_cm_x; p.tau_y -= tau_cm_y; p.tau_z -= tau_cm_z;
  }
}

/*! Read in bond information from the bonds file.
 *  Bonds file has the following structure
 *  # - comments
 *  id  type i j
 *  where id is unique bond id (starts with 0)
 *  type is bond type (positive integer)
 *  i is the index of first particle in the bond (also indexed from zero)
 *  j is the index of second particle in the bond (indexed from zero)
 *  
 *  \param bond_file name of the file containing bonds information
*/
void System::read_bonds(const string& bond_file)
{
  vector<int> types;
  vector<string> s_line;
  ifstream inp;
  inp.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  try
  {
    inp.open(bond_file.c_str());
  }
  catch (exception& e)
  {
    m_msg->msg(Messenger::ERROR,"Problem opening (bond) file "+bond_file);
    throw e;
  }
  
  string line;
  m_msg->msg(Messenger::INFO,"Reading bond information from file: "+bond_file);
  
  // Handle exclusions
  if (m_exclusions.size() == 0)
    for (int i = 0; i < this->size(); i++)
      m_exclusions.push_back(vector<int>());
  
  inp.exceptions ( std::ifstream::badbit ); // need to reset ios exceptions to avoid EOF failure of getline
  while ( getline(inp, line) )
  {
    trim(line);
    to_lower(line);
    if (line[0] != '#' && line.size() > 0)
    {
      s_line = split_line(line);
      if (s_line.size() < NUM_BOND_ATTRIB)
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define a bond. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_BOND_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the bond file.");
      }
      int id = lexical_cast<int>(s_line[0]);  // read bond id
      int tp = lexical_cast<int>(s_line[1]);  // read particle type// Handle exclusions
      int i = lexical_cast<int>(s_line[2]);   // read id of first particle
      int j = lexical_cast<int>(s_line[3]);   // read id of second particle
      if (i >= m_particles.size())
      {
        m_msg->msg(Messenger::ERROR,"Index of first particle "+lexical_cast<string>(i)+" is larger than the total number of particles.");
        throw runtime_error("Wrong particle index.");
      }
      if (j >= m_particles.size())
      {
        m_msg->msg(Messenger::ERROR,"Index of second particle "+lexical_cast<string>(j)+" is larger than the total number of particles.");
        throw runtime_error("Wrong particle index.");
      }
      m_bonds.push_back(Bond(id, tp, i, j));
      m_particles[i].bonds.push_back(id);
      m_particles[j].bonds.push_back(id);
      if (!this->in_exclusion(i,j))
        m_exclusions[i].push_back(j);
      if (!this->in_exclusion(j,i))
        m_exclusions[j].push_back(i);
      if (find(types.begin(), types.end(), tp) == types.end())
        types.push_back(tp);
    }
  }
  m_msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_bonds.size())+" bonds.");
  inp.close();
  
  m_n_bond_types = types.size();
  
  m_has_exclusions = true;
}


/*! Read in angle information from the angles file.
 *  Angles file has the following structure
 *  # - comments
 *  id  type i j k
 *  where id is unique angle id (starts with 0)
 *  type is angle type (positive integer)
 *  i is the index of first (left) particle
 *  j is the index of second (middle) particle 
 *  k is the index of the third (right) particle
 *  
 *  \param angle_file name of the file containing angles information
*/
void System::read_angles(const string& angle_file)
{
  vector<int> types;
  vector<string> s_line;
  ifstream inp;
  inp.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  try
  {
    inp.open(angle_file.c_str());
  }
  catch (exception& e)
  {
    m_msg->msg(Messenger::ERROR,"Problem opening (angle) file "+angle_file);
    throw e;
  }
  
  string line;
  m_msg->msg(Messenger::INFO,"Reading angle information from file: "+angle_file);
  
  // Handle exclusions
  if (m_exclusions.size() == 0)
    for (int i = 0; i < this->size(); i++)
      m_exclusions.push_back(vector<int>());
  
  inp.exceptions ( std::ifstream::badbit ); // need to reset ios exceptions to avoid EOF failure of getline
  while ( getline(inp, line) )
  {
    trim(line);
    to_lower(line);
    if (line[0] != '#' && line.size() > 0)
    {
      s_line = split_line(line);
      if (s_line.size() < NUM_ANGLE_ATTRIB)
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define an angle. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_ANGLE_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the angle file.");
      }
      int id = lexical_cast<int>(s_line[0]);  // read bond id
      int tp = lexical_cast<int>(s_line[1]);  // read particle type
      int i = lexical_cast<int>(s_line[2]);   // read id of first (left) particle
      int j = lexical_cast<int>(s_line[3]);   // read id of second (middle) particle
      int k = lexical_cast<int>(s_line[4]);   // read id of third (right) particle
      if (i >= m_particles.size())
      {
        m_msg->msg(Messenger::ERROR,"Index of first (m_systemleft) particle "+lexical_cast<string>(i)+" is larger than the total number of particles.");
        throw runtime_error("Wrong particle index.");
      }
      if (j >= m_particles.size())
      {
        m_msg->msg(Messenger::ERROR,"Index of second (middle) particle "+lexical_cast<string>(j)+" is larger than the total number of particles.");
        throw runtime_error("Wrong particle index.");
      }
      if (k >= m_particles.size())
      {
        m_msg->msg(Messenger::ERROR,"Index of third (right) particle "+lexical_cast<string>(k)+" is larger than the total number of particles.");
        throw runtime_error("Wrong particle index.");
      }
      m_angles.push_back(Angle(id, tp, i, j, k));
      m_particles[i].angles.push_back(id);
      m_particles[j].angles.push_back(id);
      m_particles[k].angles.push_back(id);
      if (!this->in_exclusion(i,j))
        m_exclusions[i].push_back(j);
      if (!this->in_exclusion(i,k))
        m_exclusions[i].push_back(k);
      if (!this->in_exclusion(j,i))
        m_exclusions[j].push_back(i);
      if (!this->in_exclusion(j,k))
        m_exclusions[j].push_back(k);
      if (!this->in_exclusion(k,i))
        m_exclusions[k].push_back(i);
      if (!this->in_exclusion(k,j))
        m_exclusions[k].push_back(j);
      if (find(types.begin(), types.end(), tp) == types.end())
        types.push_back(tp);
    }
  }
  m_msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_angles.size())+" angles.");
  inp.close();
  m_has_exclusions = true;
  m_n_angle_types = types.size();
}

//! Compute tangent in the direction of the neighbour with the smallest index
//! Used for simulation of actomyosin
//! \param i particle at which we want to compute tangent vector
//! \param tx x component of the tangent vector (returned)
//! \param ty y component of the tangent vector (returned)
//! \param tz z component of the tangent vector (returned)
void System::compute_tangent(int i, double& tx, double& ty, double& tz)
{
  tx = ty = tz = 0.0;
  Particle& pi = m_particles[i];
  if (pi.bonds.size() > 0)
  {
    int j = *(min_element(pi.bonds.begin(), pi.bonds.end()));
    if (j < pi.get_id())
    {
      Particle& pj = m_particles[j];
      double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
      if (m_periodic)
      {
        if (dx > m_box->xhi) dx -= m_box->Lx;
        else if (dx < m_box->xlo) dx += m_box->Lx;
        if (dy > m_box->yhi) dy -= m_box->Ly;
        else if (dy < m_box->ylo) dy += m_box->Ly;
        if (dz > m_box->zhi) dz -= m_box->Lz;
        else if (dz < m_box->zlo) dz += m_box->Lz;
      }
      double l = sqrt(dx*dx + dy*dy + dz*dz);
      tx = dx/l;
      ty = dy/l;
      tz = dz/l;
    }
  }
}
 
//! Apply period boundary conditions on three coordinate
//! \param dx x value too apply periodic boundary to (will be overwritten)
//! \param dy y value too apply periodic boundary to (will be overwritten)
//! \param dz z value too apply periodic boundary to (will be overwritten)
void System::apply_periodic(double& dx, double& dy, double& dz)
{
  if (m_periodic)
  {
    if (dx > m_box->xhi) dx -= m_box->Lx;
    else if (dx < m_box->xlo) dx += m_box->Lx;
    if (dy > m_box->yhi) dy -= m_box->Ly;
    else if (dy < m_box->ylo) dy += m_box->Ly;
    if (dz > m_box->zhi) dz -= m_box->Lz;
    else if (dz < m_box->zlo) dz += m_box->Lz;
  }
}
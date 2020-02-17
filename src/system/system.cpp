/* ***************************************************************************
 *
 *  Copyright (C) 2013-2016 University of Dundee
 *  All rights reserved. 
 *
 *  This file is part of SAMoS (Soft Active Matter on Surfaces) program.
 *
 *  SAMoS is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  SAMoS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * ****************************************************************************/

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
System::System(const string& input_filename, MessengerPtr msg, BoxPtr box) : m_msg(msg), 
                                                                             m_box(box), 
                                                                             m_mesh(Mesh()),
                                                                             m_periodic(false),
                                                                             m_force_nlist_rebuild(false),
                                                                             m_nlist_rescale(1.0),
                                                                             m_current_particle_flag(0),
                                                                             m_dt(0.0),
                                                                             m_max_mesh_iter(100),
                                                                             m_boundary_type(1),
                                                                             m_has_boundary_neighbours(false),
                                                                             m_record_force_type(false)
{
  vector<int> types;
  vector<string> s_line;
  map<string, int> column_key;
  bool has_keys = false;
  bool write_keys = false;
  bool has_boundary = false;   // flag that determines is boundary particles are present
  int id = -1, tp;
  double r;
  // Here we list keys that appear in the old input file format
  vector<string> old_keys;
  old_keys.push_back("id");  old_keys.push_back("type");   old_keys.push_back("radius");
  old_keys.push_back("x");  old_keys.push_back("y");   old_keys.push_back("z");    
  old_keys.push_back("vx");  old_keys.push_back("vy");   old_keys.push_back("vz");    
  old_keys.push_back("nx");  old_keys.push_back("ny");   old_keys.push_back("nz");
  old_keys.push_back("omega");  old_keys.push_back("length");
  old_keys.push_back("ix");  old_keys.push_back("iy");   old_keys.push_back("iz");    
  // done with this odd setup
  
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

  // Create group "all"
  m_group["all"] = make_shared<Group>(Group(0,"all"));
  m_num_groups = 1;
  
  m_msg->msg(Messenger::INFO,"Generated group 'all' containing all particles.");

  string line;
  m_msg->msg(Messenger::INFO,"Reading particle coordinates from file: "+input_filename);
  m_msg->write_config("system.input_file",input_filename);
  
  inp.exceptions ( std::ifstream::badbit ); // need to reset ios exceptions to avoid EOF failure of getline
  while ( getline(inp, line) )
  {
    trim(line);
    to_lower(line);
    if (line.size() > 0)
    {
      if (has_keys && write_keys)
      {
        for (unsigned int col_idx = 1; col_idx < s_line.size(); col_idx++)
          column_key[s_line[col_idx]] = col_idx-1;
        for (unsigned int col_idx = 1; col_idx < s_line.size(); col_idx++) 
        {
          m_msg->msg(Messenger::INFO,"Column " + lexical_cast<string>(col_idx) + " of input file is : " + s_line[col_idx] + ".");
          cout << "Column " << lexical_cast<string>(col_idx) <<  " of input file is : " << s_line[col_idx] << endl;
        }
        if (column_key.find("boundary") != column_key.end())  // Create two groups for boundary and internal particles
        {
          m_group["tissue"] = make_shared<Group>(Group(1,"tissue"));
          m_group["environment"] = make_shared<Group>(Group(2,"environment"));
          m_group["boundary"] = make_shared<Group>(Group(3,"boundary"));
          m_group["internal"] = make_shared<Group>(Group(4,"internal"));
          m_num_groups = 5;
          m_msg->msg(Messenger::INFO,"Generated groups 'tissue', 'environment', 'boundary' and 'internal' to distinguish particles beloging to the tissue or not.");
          m_msg->msg(Messenger::INFO,"Generated groups 'boundary' and 'internal' to distinguish if tissue particles are inside or on the boundary.");
          has_boundary = true;
        }
        write_keys = false;
      }
      s_line = split_line(line);
      // Some variability in input files: a lot of them have the syntax # id type ... etc. This should also be read as keys!
      if (s_line[0] == "#")
      {
	      //if (s_line.size() > 1 && (s_line[1] == "id" || s_line[1] == "type" || s_line[1] == "radius" || s_line[1] == "x"))
        if (s_line.size() > 1 && find(old_keys.begin(), old_keys.end(), s_line[1]) != old_keys.end())
        {
          m_msg->msg(Messenger::WARNING,"Input format style '# id type radius' is deprecated and will be removed in a future version.");
          //cout << "Warning! Input format style '# id type radius' is deprecated and will be removed in a future version" << endl;
          has_keys = true;
          write_keys = true;
        }
      }
      else if (s_line[0] == "keys:")
      {
	      has_keys = true;
        write_keys = true;
      }
      else if (s_line.size() < NUM_PART_ATTRIB && (!has_keys))
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define a particle. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_PART_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the input file.");
      }
      else
      {
        // read particle id
        if (has_keys)
        {
          if (column_key.find("id") != column_key.end()) 
          {
            try
            {
              id = lexical_cast<int>(s_line[column_key["id"]]);
            }
            catch(bad_lexical_cast &e)
            {
              std::cerr << "Particle id has to be zero or a positive integer." << std::endl;
              throw;
            }
          }
          else
            id++;            
        }
        else
          id = lexical_cast<int>(s_line[0]);
        
        // read particle type
        if (has_keys)
        {
          if (column_key.find("type") != column_key.end()) 
          {
            try
            {
              tp = lexical_cast<int>(s_line[column_key["type"]]);
            }
            catch(bad_lexical_cast &e)
            {
              std::cerr << "Particle type has to be a positive integer." << std::endl;
              throw;
            }
          }
          else 
            tp = 1;
        }
        else
          tp = lexical_cast<int>(s_line[1]);
        if (tp < 1)
        {
          m_msg->msg(Messenger::ERROR,"Particle type has to be positive integer.");
          throw runtime_error("Wrong particle type.");
        }
        // read radius
        if (has_keys)
        {
          if (column_key.find("radius") != column_key.end()) 
          {
            try
            {
              r = lexical_cast<double>(s_line[column_key["radius"]]);
            }
            catch (bad_lexical_cast &e)
            {
              std::cerr << "Particle radius has to be a positive floating point number." << std::endl;
              throw;
            }
          }
          else
            r = 1.0;            
        }
        else
          r = lexical_cast<double>(s_line[2]);
        Particle p(id, tp, r);
        // read x coordinate
        if (has_keys)
        {
          p.x = 0.0;
          if (column_key.find("x") != column_key.end())     p.x = lexical_cast<double>(s_line[column_key["x"]]);            
        }
        else
          p.x = lexical_cast<double>(s_line[3]);
        if (p.x < m_box->xlo || p.x > m_box->xhi)
        {
          m_msg->msg(Messenger::ERROR,"X coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
          throw runtime_error("Particle outside the box.");
        }
        // read y coordinate
        if (has_keys)
        {
          p.y = 0.0;
          if (column_key.find("y") != column_key.end())     p.y = lexical_cast<double>(s_line[column_key["y"]]);            
        }
        else
          p.y = lexical_cast<double>(s_line[4]);
        if (p.y < m_box->ylo || p.y > m_box->yhi)
        {
          m_msg->msg(Messenger::ERROR,"Y coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
          throw runtime_error("Particle outside the box.");
        }
        // read z coordinate
        if (has_keys)
        {
          p.z = 0.0;
          if (column_key.find("z") != column_key.end())     p.z = lexical_cast<double>(s_line[column_key["z"]]);            
        }
        else
          p.z = lexical_cast<double>(s_line[5]);
        if (p.z < m_box->zlo || p.z > m_box->zhi)
        {
          m_msg->msg(Messenger::ERROR,"Z coordinate of particle "+lexical_cast<string>(p.get_id())+" is outside simulation box. Please update box size.");
          throw runtime_error("Particle outside the box.");
        }
        if (find(types.begin(), types.end(), tp) == types.end())
          types.push_back(tp);
        // read v
        if (has_keys)
        {
          p.vx = 0.0; p.vy = 0.0; p.vz = 0.0;
          if (column_key.find("vx") != column_key.end())   p.vx = lexical_cast<double>(s_line[column_key["vx"]]);  
          if (column_key.find("vy") != column_key.end())   p.vy = lexical_cast<double>(s_line[column_key["vy"]]);
          if (column_key.find("vz") != column_key.end())   p.vz = lexical_cast<double>(s_line[column_key["vz"]]); 
        }
        else
        {
          p.vx = lexical_cast<double>(s_line[6]);
          p.vy = lexical_cast<double>(s_line[7]);
          p.vz = lexical_cast<double>(s_line[8]);
        }
        // read 
        if (has_keys)
        {
          p.nx = 1.0;  p.ny = 0.0;  p.nz = 0.0;
          if (column_key.find("nx") != column_key.end())    p.nx = lexical_cast<double>(s_line[column_key["nx"]]);
          if (column_key.find("ny") != column_key.end())    p.ny = lexical_cast<double>(s_line[column_key["ny"]]);
          if (column_key.find("nz") != column_key.end())    p.nz = lexical_cast<double>(s_line[column_key["nz"]]);
        }
        else
        {
          p.nx = lexical_cast<double>(s_line[9]);
          p.ny = lexical_cast<double>(s_line[10]);
          p.nz = lexical_cast<double>(s_line[11]);
        }
        // read omega
        if (has_keys)
        {
          p.omega = 0.0;
          if (column_key.find("omega") != column_key.end())   p.omega = lexical_cast<double>(s_line[column_key["omega"]]);            
        }
        else if (s_line.size() > 12)
          p.omega = lexical_cast<double>(s_line[12]);
        // read length
        p.set_length(1.0);
        if (has_keys)
        {
          if (column_key.find("length") != column_key.end()) 
            p.set_length(lexical_cast<double>(s_line[column_key["length"]]));            
        }
        else
        {
          if (s_line.size() > 13)
            p.set_length(lexical_cast<double>(s_line[13]));
        }
        // read image flags
        p.ix = 0;  p.iy = 0;  p.iz = 0;
        if (has_keys)
        {
          if (column_key.find("ix") != column_key.end())      p.ix = lexical_cast<int>(s_line[column_key["ix"]]);
          if (column_key.find("iy") != column_key.end())      p.iy = lexical_cast<int>(s_line[column_key["iy"]]);
          if (column_key.find("iz") != column_key.end())      p.iz = lexical_cast<int>(s_line[column_key["iz"]]);
        }
        else
        {
          if (s_line.size() > 16)
          {
            p.ix = lexical_cast<int>(s_line[14]);
            p.iy = lexical_cast<int>(s_line[15]);
            p.iz = lexical_cast<int>(s_line[16]);
          }
        }
        if (has_keys)
          if (column_key.find("parent") != column_key.end())   p.set_parent(lexical_cast<int>(s_line[column_key["parent"]]));
        if (has_keys)
        {
          if (column_key.find("nvx") != column_key.end())      p.Nx = lexical_cast<double>(s_line[column_key["nvx"]]);
          if (column_key.find("nvy") != column_key.end())      p.Ny = lexical_cast<double>(s_line[column_key["nvy"]]);
          if (column_key.find("nvz") != column_key.end())      p.Nz = lexical_cast<double>(s_line[column_key["nvz"]]);
        }
        if (has_keys && (column_key.find("area") != column_key.end()))  
        {
          double A0 = lexical_cast<double>(s_line[column_key["area"]]);
          p.set_default_area(A0);
          p.A0 = A0;
        }
        if (has_keys)
          if (column_key.find("mass") != column_key.end())  p.mass = lexical_cast<double>(s_line[column_key["mass"]]);
        if (has_keys)
          if (column_key.find("molecule") != column_key.end())  
          {
            p.molecule = lexical_cast<double>(s_line[column_key["molecule"]]);
            if (p.molecule < m_molecules.size())
              m_molecules[p.molecule].push_back(p.get_id());
            else 
              m_molecules.push_back(vector<int>(1,p.get_id()));
          }
        if (has_keys && (column_key.find("boundary") != column_key.end()))  
        {
          if (lexical_cast<int>(s_line[column_key["boundary"]]) != 0)  
          {
            p.boundary = true;
            m_boundary.push_back(p.get_id());
          }
          else p.boundary = false;
          if (column_key.find("in_tissue") == column_key.end())  // boundary flags given but no in_tissue flag set, assume all in tissue
            p.in_tissue = true; 
        }
        if (has_keys && (column_key.find("in_tissue") != column_key.end()))
        {
          if (lexical_cast<int>(s_line[column_key["in_tissue"]]) != 0)
          {
            p.in_tissue = true;
            if (column_key.find("boundary") == column_key.end())
            {
              m_msg->msg(Messenger::ERROR,"For tissue simulations (in_tissue flag) boundary flag need to be specified.");
              throw runtime_error("No boundary flag specified for tissue simulations.");
            }
          }
        }
        p.set_flag(m_current_particle_flag);
        m_current_particle_flag++;
        m_particles.push_back(p);
      }
    }
  }
  if (!has_keys)
    m_msg->msg(Messenger::WARNING,"Data file uses obsolete format. Some code features will not function properly. Please consider using the new input format.");
  m_msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_particles.size())+" particles.");
  m_msg->write_config("system.n_particles",lexical_cast<string>(m_particles.size()));
  inp.close();
  
  // Populate group 'all'
  // As well as boundary and internal groups, if those are present
  for (unsigned int i = 0; i < m_particles.size(); i++)
  {
    m_group["all"]->add_particle(i);
    m_particles[i].add_group("all");
    if (has_boundary)
    {
      Particle& p = m_particles[i];
      if (p.in_tissue)
      {
        m_group["tissue"]->add_particle(i);
        m_particles[i].add_group("tissue");
        if (p.boundary)
        {
          m_group["boundary"]->add_particle(i);
          m_particles[i].add_group("boundary");
        }
        else
        {
          m_group["internal"]->add_particle(i);
          m_particles[i].add_group("internal");
        }
      }
      else
      {
        m_group["environment"]->add_particle(i);
        m_particles[i].add_group("environment");
      }
    }
  }
  // Make sure molecules exist
  if (m_molecules.size() == 0)
    for (unsigned int i = 0; i < m_particles.size(); i++)
    {
      Particle& p = m_particles[i];
      m_molecules.push_back(vector<int>(1,p.get_id()));
    }
      
  
  m_n_types = types.size();
   
  m_msg->msg(Messenger::INFO,"There are " + lexical_cast<string>(m_n_types) + " distinct particle types in the system.");
  m_msg->write_config("system.n_types",lexical_cast<string>(m_n_types));
  
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
  if (m_group.find(name) == m_group.end())
    m_group[name] = make_shared<Group>(Group(m_num_groups++, name));
  map<string, vector<bool> > to_add;
  if (param.find("type") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["type"].push_back(false);
    m_msg->msg(Messenger::INFO,"Adding particle to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_type() == lexical_cast<int>(param["type"]))
        to_add["type"][i] = true;
  }
  // Make group of particles that are not of a given type
  if (param.find("not_type") != param.end())
  {
    for (unsigned int i = 0; i < m_particles.size(); i++) to_add["not_type"].push_back(false);
    m_msg->msg(Messenger::INFO,"Adding particle to group "+name+".");
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i].get_type() != lexical_cast<int>(param["not_type"]))
        to_add["not_type"][i] = true;
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
      p.add_group(name);
    }
  }
  
  m_msg->msg(Messenger::INFO,"Added "+lexical_cast<string>(m_group[name]->get_size())+" particles to group "+name+".");
  m_msg->add_config("system.group.name",name);
      
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
  if (p.boundary) m_boundary.push_back(p.get_id());
  if (p.molecule < m_molecules.size())
    m_molecules[p.molecule].push_back(p.get_id());
  else 
    m_molecules.push_back(vector<int>(1,p.get_id()));
  // We need to force neighbour list rebuild
  m_force_nlist_rebuild = true;
  m_current_particle_flag++;
}

/*! Remove particle from the system
 *  \param id Id of particle to remove
 */ 
void System::remove_particle(int id)
{
  Mesh& mesh = this->get_mesh();
  Particle& pi = m_particles[id];
  int mol_id = pi.molecule;
  if (pi.boundary)
  {
    Particle& pj = m_particles[pi.boundary_neigh[0]];
    Particle& pk = m_particles[pi.boundary_neigh[1]];
    pj.boundary_neigh[(pj.boundary_neigh[0] == id) ? 0 : 1] = pk.get_id();
    pk.boundary_neigh[(pk.boundary_neigh[0] == id) ? 0 : 1] = pj.get_id();  
    // Vertex vi = mesh.get_vertices()[id];
    // if (vi.neigh.size() > 0)
    // {
    //   //if (!(vi.neigh[0] == pj.get_id() || vi.neigh[0] == pk.get_id()))
    //   //  throw runtime_error("Particle removal. Something is wrong with the ordering of one of the boundary vertices.");
    //   int first_neigh = vi.neigh[0] == pj.get_id() ? pj.get_id() : pk.get_id();
    //   int last_neigh  = vi.neigh[0] == pj.get_id() ? pk.get_id() : pj.get_id();
    //   m_particles[first_neigh].boundary_neigh[m_particles[first_neigh].boundary_neigh[0] == id ? 0 : 1] = vi.neigh[1];
    //   for (int i = 1; i < vi.neigh.size() - 1; i++)
    //   {
    //     m_particles[vi.neigh[i]].boundary = true;
    //     m_particles[vi.neigh[i]].boundary_neigh[0] = vi.neigh[i - 1];
    //     m_particles[vi.neigh[i]].boundary_neigh[1] = vi.neigh[i + 1];
    //   }
    //   m_particles[last_neigh].boundary_neigh[m_particles[last_neigh].boundary_neigh[0] == id ? 0 : 1] = vi.neigh[vi.neigh.size()-1];
    // }
    // else
    // {
    //   pj.boundary_neigh[(pj.boundary_neigh[0] == id) ? 0 : 1] = pk.get_id();
    //   pk.boundary_neigh[(pk.boundary_neigh[0] == id) ? 0 : 1] = pj.get_id();  
    // }
  }
  vector<int>::iterator it_m = find(m_molecules[mol_id].begin(), m_molecules[mol_id].end(),id);
  m_molecules[mol_id].erase(it_m);
  if (m_molecules[mol_id].size() == 0)
    m_molecules.erase(m_molecules.begin() + mol_id);
  m_particles.erase(m_particles.begin() + id);
  // Shift down all ids
  for (unsigned int i = 0; i < m_particles.size(); i++)
  {
    Particle& p = m_particles[i];
    if (p.get_id() > id)
      p.set_id(p.get_id() - 1);
    if (p.molecule > mol_id)
      p.molecule--;
    if (p.boundary)
    {
      if (p.boundary_neigh[0] > id) p.boundary_neigh[0]--;
      if (p.boundary_neigh[1] > id) p.boundary_neigh[1]--;
    }
  }
  for (unsigned int i = 0; i < m_molecules.size(); i++)
  {
    for (unsigned int j = 0; j < m_molecules[i].size(); j++)
      if (m_molecules[i][j] > mol_id) 
       m_molecules[i][j]--; 
  }
  // Update all groups
  for(map<string, GroupPtr>::iterator it_g = m_group.begin(); it_g != m_group.end(); it_g++)
    (*it_g).second->shift(id);
  
  vector<int>::iterator it = find(m_boundary.begin(), m_boundary.end(),id);
  if (it != m_boundary.end()) m_boundary.erase(it);
  for (unsigned int i = 0; i < m_boundary.size(); i++)  
    if (m_boundary[i] > id) m_boundary[i]--;
  
  m_force_nlist_rebuild = true;
}

/*! Change group of the particle
 *  \param p Particle whose group we want to change
 *  \param old_group Old group
 *  \param new_group New group 
*/
void System::change_group(int id, const string& old_group, const string& new_group)
{
  Particle& p = m_particles[id];
  
  if (old_group == "all" && new_group != "all")
  {
    m_msg->msg(Messenger::ERROR,"Particle "+lexical_cast<string>(p.get_id())+": Cannot change from \"all\" to other group ("+new_group+").");
    throw runtime_error("Can't change out of \"all\" group.");
  }
  list<string>::iterator it_g = find(p.groups.begin(),p.groups.end(),old_group);
  if (it_g == p.groups.end())
  {
    m_msg->msg(Messenger::ERROR,"Particle "+lexical_cast<string>(p.get_id())+" does not belong to group "+old_group+".");
    cout << p;
    throw runtime_error("Unknown particle group.");
  }
  else
    p.groups.erase(it_g);
  
  if (new_group != "all") // It's already in "all"
    p.add_group(new_group);
  
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
  m_msg->write_config("system.bond_file",bond_file);
  
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
      bool exclude = true;
      s_line = split_line(line);
      if (s_line.size() < NUM_BOND_ATTRIB)
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define a bond. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_BOND_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the bond file.");
      }
      unsigned int id = lexical_cast<int>(s_line[0]);  // read bond id
      unsigned int tp = lexical_cast<int>(s_line[1]);  // read particle type// Handle exclusions
      unsigned int i = lexical_cast<int>(s_line[2]);   // read id of first particle
      unsigned int j = lexical_cast<int>(s_line[3]);   // read id of second particle
      if (s_line.size() > 4)
        if (s_line[4] == "include_nb")
          exclude = false;
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
      if (exclude)
      {
        if (!this->in_exclusion(i,j))
          m_exclusions[i].push_back(j);
        if (!this->in_exclusion(j,i))
          m_exclusions[j].push_back(i);
      }
      if (find(types.begin(), types.end(), tp) == types.end())
        types.push_back(tp);
    }
  }
  m_msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_bonds.size())+" bonds.");
  m_msg->write_config("system.n_bonds",lexical_cast<string>(m_bonds.size()));
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
  m_msg->write_config("system.angle_file",angle_file);
  
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
      bool exclude = true;
      s_line = split_line(line);
      if (s_line.size() < NUM_ANGLE_ATTRIB)
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define an angle. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(NUM_ANGLE_ATTRIB) + " expected.");
        throw runtime_error("Insufficient number of parameters in the angle file.");
      }
      unsigned int id = lexical_cast<int>(s_line[0]);  // read bond id
      unsigned int tp = lexical_cast<int>(s_line[1]);  // read particle type
      unsigned int i = lexical_cast<int>(s_line[2]);   // read id of first (left) particle
      unsigned int j = lexical_cast<int>(s_line[3]);   // read id of second (middle) particle
      unsigned int k = lexical_cast<int>(s_line[4]);   // read id of third (right) particle
      if (s_line.size() > 5)
        if (s_line[5] == "include_nb")
          exclude = false;
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
      if (exclude)
      {
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
      }
      if (find(types.begin(), types.end(), tp) == types.end())
        types.push_back(tp);
    }
  }
  m_msg->msg(Messenger::INFO,"Read data for "+lexical_cast<string>(m_angles.size())+" angles.");
  m_msg->write_config("system.n_angles",lexical_cast<string>(m_angles.size()));
  inp.close();
  m_has_exclusions = true;
  m_n_angle_types = types.size();
}

/*! Read in connectivity information for boundary particles in tissue simulations
 *  Boundary neighbours file has the following structure
 *  # - comments
 *  id i j
 *  where id is unique bond id (starts with 0)
 *  i is the index of boundary particle
 *  j is the index of its neighbour
 *  \note only one pair has to be specified, the other one is set atomatically
 *  
 *  \param bound_neigh_file name of the file containing boundary connectivity information 
*/
void System::read_boundary_neighbours(const string& bound_neigh_file)
{
  vector<string> s_line;
  ifstream inp;
  inp.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
  try
  {
    inp.open(bound_neigh_file.c_str());
  }
  catch (exception& e)
  {
    m_msg->msg(Messenger::ERROR,"Problem opening (boundary connectivity) file "+bound_neigh_file);
    throw e;
  }
  
  string line;
  m_msg->msg(Messenger::INFO,"Reading boundary connectivity information from file: "+bound_neigh_file);
  m_msg->write_config("system.bound_neigh_file",bound_neigh_file);
  
  
  inp.exceptions ( std::ifstream::badbit ); // need to reset ios exceptions to avoid EOF failure of getline
  while ( getline(inp, line) )
  {
    trim(line);
    to_lower(line);
    if (line[0] != '#' && line.size() > 0)
    {
      s_line = split_line(line);
      if (s_line.size() < 3)
      {
        m_msg->msg(Messenger::ERROR,"Insufficient number of parameters to define a bondary neighbours. " + lexical_cast<string>(s_line.size()) + " given, but " + lexical_cast<string>(3) + " expected.");
        throw runtime_error("Insufficient number of parameters in the boundary connectivity file.");
      }
      unsigned int i = lexical_cast<int>(s_line[1]);   // read id of first particle
      unsigned int j = lexical_cast<int>(s_line[2]);   // read id of second particle
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
      Particle& pi = m_particles[i];
      Particle& pj = m_particles[j];
      if (!(pi.boundary && pj.boundary))
      {
        m_msg->msg(Messenger::ERROR,"Both particles have to be boundary particles in order to define them as boundary neighbours.");
        throw runtime_error("Inconsistent boundary connectivity information.");
      }
      pi.boundary_neigh.push_back(j);
      pj.boundary_neigh.push_back(i);
    }
  }
  m_has_boundary_neighbours = true;
  inp.close();
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
    int min_j = 10000000;   // some really large number 
    for (vector<int>::iterator it_b = pi.bonds.begin(); it_b != pi.bonds.end(); it_b++)
    {
      Bond& bond = m_bonds[*it_b];
      int j = (bond.i == pi.get_id()) ? bond.j : bond.i;
      if (j < min_j) min_j = j;
    }
    Particle& pj = m_particles[min_j];
    double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
    this->apply_periodic(dx,dy,dz);
    double l = sqrt(dx*dx + dy*dy + dz*dz);
    double fact = 1.0;
    if (min_j > pi.get_id()) fact = -1.0;
    tx = fact*dx/l;
    ty = fact*dy/l;
    tz = fact*dz/l;
  }
}


//! Compute system area by adding up areas of all cells (makes sense only for cell systems)
double System::compute_area()
{
  Mesh& mesh = this->get_mesh();
  double area = 0.0;
  int num = 0;
  for (int i = 0; i < mesh.size(); i++)
  {
    Vertex& V = mesh.get_vertices()[i];
    if (!V.boundary && V.attached)
    {
      area += V.area;
      num++;
    }
  }
  if (num > 0)
    return area/num;
  else
    return 0.0;
}

//! Compute average perimeter of cells in a cell system
double System::compute_average_perimeter()
{
  Mesh& mesh = this->get_mesh();
  double  perim = 0.0;
  int N = mesh.size();
  int num = 0;
  for (int i = 0; i < N; i++)
  {
    Vertex& V = mesh.get_vertices()[i];
    if (!V.boundary && V.attached)
    {
      perim += V.perim;
      num++;
    }
  }
  if (N > 0)
    return perim/num;
  else
    return 0.0;
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

//! Make sure that all group information on particles matches group information in the lists
//! \param group group name
bool System::group_ok(const string& group)
{
  int N = m_group[group]->get_size();
  vector<int> particles = m_group[group]->get_particles();
  for (int i = 0; i < N; i++)
  {
    int pi = particles[i];
    Particle& p = m_particles[pi];
    list<string>::iterator it_g = find(p.groups.begin(),p.groups.end(),group);
    if (it_g == p.groups.end())
    {
      cout << "For group : " << group << endl;
      cout << p;
      return false;
    }
  }
  return true;
}

//! Enforce period boundary conditions onto a particle
//! \param p reference to the particle
void System::enforce_periodic(Particle& p)
{
  if (m_periodic)
  {
    if (p.x <= m_box->xlo) 
    {
      p.x += m_box->Lx;
      p.ix--;
    }
    else if (p.x > m_box->xhi)
    {
      p.x -= m_box->Lx;
      p.ix++;
    }
    if (p.y <= m_box->ylo)
    {
      p.y += m_box->Ly;
      p.iy--;
    }
    else if (p.y > m_box->yhi)
    {
      p.y -= m_box->Ly;
      p.iy++;
    }
    if (p.z <= m_box->zlo)
    {
      p.z += m_box->Lz;
      p.iz--;
    }
    else if (p.z > m_box->zhi)
    {
      p.z -= m_box->Lz;
      p.iz++;
    }
  }
}

//! Update mesh for tissue simulations
void System::update_mesh()
{
  if (m_mesh.size() > 0)
  {
    for (int i = 0; i < this->size(); i++)
    {
      Particle& p = m_particles[i];
      m_mesh.update(p);
    }
    bool converged = false;
    int iter = 0;
    while(!converged && iter <= m_max_mesh_iter)
    {
      //cout << "iteration : " << iter++ << endl;
      converged = true;
      //converged = converged && m_mesh.remove_obtuse_boundary();
      /*
      vector<Vector3d> vecs = m_mesh.fix_obtuse_boundary();
      if (vecs.size() > 0)
      {
        for (vector<Vector3d>::iterator it_v = vecs.begin(); it_v != vecs.end(); it_v++)
        {
          Vector3d& v = *it_v;
          Particle p(this->size(),m_boundary_type,1.0);  
          p.x = v.x; 
          p.y = v.y; 
          p.z = v.z; 
          this->add_particle(p);
        }
        converged = false;
      }
      */
      converged = converged && m_mesh.equiangulate();
      //converged = converged && m_mesh.remove_edge_triangles();
      m_mesh.update_dual_mesh();
      m_mesh.update_face_properties();
      if (m_mesh.has_obtuse_boundary())
        this->set_force_nlist_rebuild(true);
      if (m_mesh.has_dangling_vertices())
        this->set_force_nlist_rebuild(true);
      //converged = converged && m_mesh.equiangulate();
      iter++;
    }
    if (iter >= m_max_mesh_iter)
    {
      cout << "Exceeded maximum number of iterations in boundary build. Most likely something is wrong with input paramters. Results will not be reliable." << endl;
      throw runtime_error("Exceeded maximum number of iterations in boundary build.");
    }
    for (int i = 0; i < m_mesh.size(); i++)
    {
      m_mesh.order_dual(i);
      m_mesh.dual_perimeter(i);
      m_mesh.dual_area(i); 
    }
  }
}

//! Compute centre of mass for a molecule
//! \param mol_id id of the molecule
//! \param xcm x-coordinate of the centre of mass
//! \param ycm y-coordinate of the centre of mass
//! \param zcm z-coordinate of the centre of mass
void System::molecule_cm(int mol_id, double& xcm, double& ycm, double& zcm)
{
  vector<int>& mol = this->get_mol_particles(mol_id);
  Particle& p0 = this->get_particle(mol[0]);
  int M = mol.size();
  xcm = 0.0; ycm = 0.0; zcm = 0.0;
  for (int i = 1; i < M; i++)
  {
    Particle& pi = this->get_particle(mol[i]);
    double dx = pi.x - p0.x, dy = pi.y - p0.y, dz = pi.z - p0.z;
    this->apply_periodic(dx,dy,dz);
    xcm += dx;  ycm += dy;  zcm += dz;
  }
  xcm = xcm/M + p0.x; ycm = ycm/M + p0.y; zcm = zcm/M + p0.z;
  this->apply_periodic(xcm,ycm,zcm);
}
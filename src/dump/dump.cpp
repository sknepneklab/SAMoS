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
 * \file dump.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 24-Oct-2013
 * \brief Implementation of the functions for dumping system state
 */ 

#include "dump.hpp"

//! Simple helper function to write an integer.
/*! 
 *  \param file file to write to
 *  \param val integer to write
 *  \note This function is adopted from HOOMD code (file: DCDDumpWriter.cc)
 */
static void write_int(ofstream &file, unsigned int val)
{
  file.write((char *)&val, sizeof(unsigned int));
}

/*! Construct a specific dump object
 *  \param sys Pointer to a System object
 *  \param msg Handles system wide messages
 *  \param nlist Pointer to the global neighbour list
 *  \param filename Base file name for the output 
 *  \param params list of parameters that control dump (e.g., frequency, output type, etc.)
*/
Dump::Dump(SystemPtr sys, MessengerPtr msg, NeighbourListPtr nlist, const string& filename, pairs_type& params) : m_system(sys), 
                                                                                                                  m_msg(msg), 
                                                                                                                  m_file_name(filename), 
                                                                                                                  m_params(params),
                                                                                                                  m_time_step_offset(0),
                                                                                                                  m_print_header(false),
                                                                                                                  m_print_keys(false),
                                                                                                                  m_compress(false),
                                                                                                                  m_output_dual(false),
                                                                                                                  m_dual_boundary(false),
                                                                                                                  m_include_bonds(false),
                                                                                                                  m_include_mesh(false), 
                                                                                                                  m_group("all"),
                                                                                                                  m_directory(".")
{
  m_type_ext["velocity"] = "vel";
  m_type_ext["xyz"] = "xyz";
  m_type_ext["full"] = "dat";
  m_type_ext["input"] = "inp";
  m_type_ext["dcd"] = "dcd";
  m_type_ext["director"] = "dir";
  m_type_ext["xyzv"] = "xyzv";
  m_type_ext["xyzc"] = "xyzc";
  m_type_ext["mol2"] = "mol2";
  m_type_ext["contact"] = "con";
  m_type_ext["face"] = "fc";
  m_type_ext["mesh"] = "off";
  m_type_ext["vtp"] = "vtp";
  m_type_ext["boundary"] = "boundary";
  m_type_ext["ajm"] = "ajm";
  
  string fname = filename;
  
  replace_all(fname,".","_");
  
  if (nlist)
    m_nlist = nlist;
  else
    m_msg->msg(Messenger::WARNING,"Neighbour list has not been specified. Contact network won't be produced.");
    
  if (params.find("type") == params.end())
  {
    m_msg->msg(Messenger::ERROR,"No dump type specified.");
    throw runtime_error("No dump type specified.");
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Dump type set to "+params["type"]);
    m_type = params["type"];
  }
  m_msg->write_config("dump."+fname+".type",m_type);
  if (params.find("directory") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"No dump directory specified. Using current directory.");
    m_directory = ".";
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Dump directory set to "+params["directory"]);
    m_directory = params["directory"];
  }
  m_msg->write_config("dump."+fname+".directory",m_directory);
  if (m_type_ext.find(m_type) == m_type_ext.end())
  {
    m_msg->msg(Messenger::ERROR,"Unsupported dump type "+m_type+".");
    throw runtime_error("Unsupported dump type");
  }
  else
  {
    m_ext = m_type_ext[m_type];
    m_msg->msg(Messenger::INFO,"Dump will be sent to the file with base name "+m_file_name+" with extension "+m_ext+".");
  }
  m_msg->write_config("dump."+fname+".extension",m_ext);
  if (params.find("start") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"No starting step for dump specified. Using default 0, i.e., start dumping after at the beginning of the simulation.");
    m_start = 0;
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Dumping will start after "+params["start"]+" time steps.");
    m_start = lexical_cast<int>(params["start"]);
  }
  m_msg->write_config("dump."+fname+".start",lexical_cast<string>(m_start));
  if (params.find("freq") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"No dump frequency specified. Using default of dumping each 100 time steps.");
    m_freq = 100;
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Dump will be produced every "+params["freq"]+" time steps.");
    m_freq = lexical_cast<int>(params["freq"]);
  }
  m_msg->write_config("dump."+fname+".freq",lexical_cast<string>(m_freq));
  if (params.find("compress") != params.end())
  {
    m_msg->msg(Messenger::INFO,"Output data will be compressed.");
    m_msg->write_config("dump."+fname+".compress","true");
    m_out.push(bo::gzip_compressor(9)); 
    m_compress = true;
  }
  else
    m_msg->write_config("dump."+fname+".compress","false");  
  if (params.find("multi") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"All time steps will be concatenated to a single file.");
    m_multi_print = false;
    string file_name = m_directory+"/"+m_file_name+"."+m_ext;
    if (m_compress)
    {
      m_file.open(file_name.c_str(), std::ios_base::out | std::ios_base::binary);
      m_ext += ".gz";
    }
    else
      m_file.open(file_name.c_str()); 
    m_out.push(m_file);
  }
  else
  {
    m_msg->msg(Messenger::INFO,"Each dumped time step will be stored in a separate file labelled by the time step.");
    m_multi_print = true;
  }
  if (params.find("header") != params.end())
  {
    m_print_header = true;
    m_msg->msg(Messenger::INFO,"Include info header into the dump file.");
    m_msg->write_config("dump."+fname+".header","true");
  }
  if (params.find("keys") != params.end())
  {
    m_print_keys = true;
    m_msg->msg(Messenger::INFO,"Include keys line into the dump file.");
    m_msg->write_config("dump."+fname+".keys","true");
  }
  if (params.find("dual") != params.end())
  {
    m_output_dual = true;
    m_msg->msg(Messenger::INFO,"VTP file will contain dual mesh.");
    m_msg->write_config("dump."+fname+".dual","true");
  }
  if (params.find("dual_boundary") != params.end())
  {
    m_dual_boundary = true;
    m_msg->msg(Messenger::INFO,"VTP file will contain boundary faces for dual mesh.");
    m_msg->write_config("dump."+fname+".dual_boundary","true");
  }
  if (params.find("step_offset") != params.end())
  {
    m_time_step_offset = lexical_cast<int>(params["step_offset"]);;
    m_msg->msg(Messenger::INFO,"Dump steps will be offset by "+params["step_offset"]+".");
    m_msg->write_config("dump."+fname+".time_step_offset",params["step_offset"]);
  }
  if (params.find("include_bonds") != params.end())
  {
    m_include_bonds = true;
    m_msg->msg(Messenger::INFO,"Bonds will be included in the dump file.");
    m_msg->write_config("dump."+fname+".include_bonds","true");
  }
  if (params.find("include_mesh") != params.end())
  {
    m_include_mesh = true;
    m_msg->msg(Messenger::INFO,"Mesh will be included in the dump file.");
    m_msg->write_config("dump."+fname+".include_mesh","true");
  }
  if ((params.find("record_force_types") != params.end()) && (m_type == "vtp"))
  {
    m_system->set_record_force_type(true);
    m_msg->msg(Messenger::INFO,"Force component for each force type will be recorded in VTP files.");
    m_msg->write_config("dump."+fname+".record_force_type","true");
  }
  if (m_type == "vtp" && m_include_bonds && m_include_mesh)
  {
    m_msg->msg(Messenger::ERROR,"Including bonds and mesh in a VTP file would produce an inconsistent file. Therefore, it is not allowed.");
    runtime_error("Unsupported output parameters for VTP files.");
  }
  if (params.find("id") != params.end())
    m_msg->add_config("dump."+fname+".quantity","id");
  if (params.find("tp") != params.end())
    m_msg->add_config("dump."+fname+".quantity","tp");
  if (params.find("flag") != params.end())
    m_msg->add_config("dump."+fname+".quantity","flag");
  if (params.find("radius") != params.end())
    m_msg->add_config("dump."+fname+".quantity","radius");
  if (params.find("coordinate") != params.end())
    m_msg->add_config("dump."+fname+".quantity","coordinate");
  if (params.find("velocity") != params.end())
    m_msg->add_config("dump."+fname+".quantity","velocity");
  if (params.find("force") != params.end())
    m_msg->add_config("dump."+fname+".quantity","force");
  if (params.find("director") != params.end())
    m_msg->add_config("dump."+fname+".quantity","director");
  if (params.find("omega") != params.end())
    m_msg->add_config("dump."+fname+".quantity","omega");
  if (params.find("image_flags") != params.end())
    m_msg->add_config("dump."+fname+".quantity","image_flags");
  if (params.find("parent") != params.end())
    m_msg->add_config("dump."+fname+".quantity","parent");
  if (params.find("a0") != params.end())
    m_msg->add_config("dump."+fname+".quantity","a0");
  if (params.find("cell_area") != params.end())
    m_msg->add_config("dump."+fname+".quantity","cell_area");
  if (params.find("cell_perim") != params.end())
    m_msg->add_config("dump."+fname+".quantity","cell_perim");
  if (params.find("boundary") != params.end())
    m_msg->add_config("dump."+fname+".quantity","boundary");
  if (params.find("group") != params.end())
    m_group = params["group"];
  m_msg->add_config("dump."+fname+".group",m_group);
  if (!m_system->has_group(m_group))
  {
    m_msg->msg(Messenger::ERROR,"Dump: Group "+m_group+" is not defined.");
    throw runtime_error("Attempted to dump group that does not exist.");
  }
  if (m_type == "vtp" && m_group != "all")
    m_msg->msg(Messenger::WARNING,"Group is not set to \all\". No bond and mesh info will be included into the VTP file. ");
  if (m_type == "mol2" && m_group != "all")
    m_msg->msg(Messenger::WARNING,"MOL2 output does not support groups. All particles will be printed out.");
  if (m_output_dual && m_group != "all")
  {
    m_msg->msg(Messenger::ERROR,"Output of the dual mesh does not support groups.");
    throw runtime_error("Groups not suported for dual dump.");
  }
  
  if (m_type == "dcd")
  {
    if (m_group != "all")
    {
      m_msg->msg(Messenger::ERROR,"DCD files does not group specific dump.");
      throw runtime_error("Attempted to dump groups into a DCD output.");
    }
    if (m_compress)
    {
      m_msg->msg(Messenger::ERROR,"DCD files does not support data compression.");
      throw runtime_error("Attempted to compress DCD output.");
    }
    if (m_multi_print)
    {
      m_msg->msg(Messenger::ERROR,"DCD files cannot be split. mulifile flag is not permitted.");
      throw runtime_error("Attempted to spread DCD output into multiple files.");
    }
    else
    {
      write_int(m_file, 84);
      
      // the next 4 bytes in the file must be "CORD"
      char cord_data[] = "CORD";
      m_file.write(cord_data, 4);
      write_int(m_file, 0);      // Number of frames in file, none written yet
      write_int(m_file, m_start); // Starting timestep
      write_int(m_file, m_freq);  // Timesteps between frames written to the file
      write_int(m_file, 0);      // Number of timesteps in simulation
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);                  // timestep (unused)
      write_int(m_file, 1);                  // include unit cell
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 0);
      write_int(m_file, 24); // Pretend to be CHARMM version 24
      write_int(m_file, 84);
      write_int(m_file, 164);
      write_int(m_file, 2);
    
      char title_string[81];
      char remarks[] = "Created by SAMoS";
      strncpy(title_string, remarks, 80);
      title_string[79] = '\0';
      m_file.write(title_string, 80);
  
      char time_str[81];
      time_t cur_time = time(NULL);
      tm *tmbuf=localtime(&cur_time);
      strftime(time_str, 80, "REMARKS Created  %d %B, %Y at %H:%M", tmbuf);
      m_file.write(time_str, 80);
      
      write_int(m_file, 164);
      write_int(m_file, 4);
      write_int(m_file, m_system->size());
      write_int(m_file, 4);
      
      // check for errors
      if (!m_file.good())
      {
        m_msg->msg(Messenger::ERROR,"Error writing DCD header");
        throw runtime_error("Error writing DCD file");
      }
    }
  }
}

//! Do actual dump.
//! \param step current time step
void Dump::dump(int step)
{
  if (step < m_start)
    return;
  if (step % m_freq != 0)
    return;
  if (m_multi_print)
  {
    string file_name = m_directory+"/"+m_file_name+"_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+"."+m_ext;
    if (m_compress)
    {
      file_name += ".gz";
      if (m_ext != "vtp" && m_ext != "ajm")
        m_file.open(file_name.c_str(),std::ios_base::out | std::ios_base::binary);
    }
    else
    {  
      if (m_ext != "vtp" && m_ext != "ajm")
        m_file.open(file_name.c_str()); 
    }
    if (m_ext != "vtp" && m_ext != "ajm")
      m_out.push(m_file);
  }
  
  if (m_type == "xyz")
    this->dump_xyz();
  else if (m_type == "full")
    this->dump_data();
  else if (m_type == "input")
    this->dump_input();
  else if (m_type == "velocity")
    this->dump_velocity();
  else if (m_type == "dcd")
    this->dump_dcd();
  else if (m_type == "director")
    this->dump_director();
  else if (m_type == "xyzv")
    this->dump_xyzv();
  else if (m_type == "xyzc")
    this->dump_xyzc();
  else if (m_type == "mol2")
    this->dump_mol2();
  else if (m_type == "contact")
    this->dump_contact();
  else if (m_type == "face")
    this->dump_faces();
  else if (m_type == "mesh")
    this->dump_mesh();
  else if (m_type == "boundary")
    this->dump_boundary();
  else if (m_type == "ajm")
    this->dump_ajm(step);
#ifdef HAS_VTK  
  else if (m_type == "vtp")
    this->dump_vtp(step);
#endif
  
  if (m_multi_print)
    if (m_ext != "vtp" && m_ext != "ajm")
    {
      m_out.pop();
      m_file.close();
    }
}

// Private functions 

//! Dump coordinated in a DCD file
void Dump::dump_dcd()
{
  double Lx = m_system->get_box()->Lx;
  double Ly = m_system->get_box()->Ly;
  double Lz = m_system->get_box()->Lz;
  int N = m_system->size();
  float* buffer = new float[N];   // buffer for printing out coordinates; note: DCD requires floats
  double unitcell[6];
  // set box dimensions (not really important 
  unitcell[0] = 100.0f;
  unitcell[2] = 100.0f;
  unitcell[5] = 100.0f;
  // box angles are 90 degrees
  unitcell[1] = 0.0f;
  unitcell[3] = 0.0f;
  unitcell[4] = 0.0f;
  write_int(m_file, 48);
  m_file.write((char *)unitcell, 48);
  write_int(m_file, 48);
  // check for errors
  if (!m_file.good())
  {
    m_msg->msg(Messenger::ERROR,"Error writing into DCD file.");
    throw runtime_error("Error writing DCD file");
  } 
  // Prepare x coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_params.find("unwrap") != m_params.end())
      buffer[i] = static_cast<float>(p.x+p.ix*Lx);
    else
      buffer[i] = static_cast<float>(p.x);
  }
  // write x coords
  write_int(m_file, N * sizeof(float));
  m_file.write((char *) buffer, N * sizeof(float));
  write_int(m_file, N * sizeof(float));
  
  // Prepare y coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_params.find("unwrap") != m_params.end())
      buffer[i] = static_cast<float>(p.y+p.iy*Ly);
    else
      buffer[i] = static_cast<float>(p.y);
  }
  // write y coords
  write_int(m_file, N * sizeof(float));
  m_file.write((char *) buffer, N * sizeof(float));
  write_int(m_file, N * sizeof(float));
  
  // Prepare z coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_params.find("unwrap") != m_params.end())
      buffer[i] = static_cast<float>(p.z+p.iz*Lz);
    else
      buffer[i] = static_cast<float>(p.z);
  }
  // write z coords
  write_int(m_file, N * sizeof(float));
  m_file.write((char *) buffer, N * sizeof(float));
  write_int(m_file, N * sizeof(float));
  
  delete [] buffer;
}

//! Dump particle coordinates in the common XYZ format
void Dump::dump_xyz()
{
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  m_out << N << endl;
  m_out << "Generated by SAMoS code." << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    m_out << format("%3d\t%10.6f\t%10.6f\t%10.6f") % p.get_type() % p.x % p.y % p.z << endl;
  }
}

//! Dump selected set of data
void Dump::dump_data()
{
  double Lx = m_system->get_box()->Lx;
  double Ly = m_system->get_box()->Ly;
  double Lz = m_system->get_box()->Lz;
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  Mesh& mesh = m_system->get_mesh();
  if (m_print_header)
  {
    m_out << "# ";
    if (m_params.find("id") != m_params.end())
      m_out << " id ";
    if (m_params.find("tp") != m_params.end())
      m_out << " type ";
    if (m_params.find("flag") != m_params.end())
      m_out << " flag ";
    if (m_params.find("radius") != m_params.end())
      m_out << " radius ";
    if (m_params.find("coordinate") != m_params.end())
      m_out << " x  y  z ";
    if (m_params.find("velocity") != m_params.end())
      m_out << " vx  vy  vz ";
    if (m_params.find("force") != m_params.end())
      m_out << " fx  fy  fz ";
    if (m_params.find("director") != m_params.end())
      m_out << " nx  ny  nz ";
    if (m_params.find("omega") != m_params.end())
      m_out << " omega ";
    if (m_params.find("image_flags") != m_params.end())
      m_out << " ix  iy  iz ";
    if (m_params.find("normal") != m_params.end())
      m_out << " nvx  nvy  nvz ";
    if (m_params.find("parent") != m_params.end())
      m_out << " parent ";
    if (m_params.find("area") != m_params.end())
      m_out << " area ";
    if (m_params.find("cell_area") != m_params.end())
      m_out << " cell_area ";
    if (m_params.find("cell_perim") != m_params.end())
      m_out << " cell_perim ";
    if (m_params.find("cont_num") != m_params.end())
      m_out << " cont_num ";
    if (m_params.find("boundary") != m_params.end())
      m_out << " boundary ";
    if (m_params.find("in_tissue") != m_params.end())
      m_out << " in_tissue ";  
    if (m_params.find("molecule") != m_params.end())
      m_out << " molecule ";
    if (m_params.find("shape_param") != m_params.end())
      m_out << " shape_param";
    m_out << endl;
  }
  if (m_print_keys)
  {
    m_out << "keys: ";
    if (m_params.find("id") != m_params.end())
      m_out << " id ";
    if (m_params.find("tp") != m_params.end())
      m_out << " type ";
    if (m_params.find("flag") != m_params.end())
      m_out << " flag ";
    if (m_params.find("radius") != m_params.end())
      m_out << " radius ";
    if (m_params.find("coordinate") != m_params.end())
      m_out << " x  y  z ";
    if (m_params.find("velocity") != m_params.end())
      m_out << " vx  vy  vz ";
    if (m_params.find("force") != m_params.end())
      m_out << " fx  fy  fz ";
    if (m_params.find("director") != m_params.end())
      m_out << " nx  ny  nz ";
    if (m_params.find("omega") != m_params.end())
      m_out << " omega ";
    if (m_params.find("image_flags") != m_params.end())
      m_out << " ix  iy  iz ";
    if (m_params.find("normal") != m_params.end())
      m_out << " nvx  nvy  nvz ";
    if (m_params.find("parent") != m_params.end())
      m_out << " parent ";
    if (m_params.find("area") != m_params.end())
      m_out << " area ";
    if (m_params.find("cell_area") != m_params.end())
      m_out << " cell_area ";
    if (m_params.find("cell_perim") != m_params.end())
      m_out << " cell_perim ";
    if (m_params.find("cont_num") != m_params.end())
      m_out << " cont_num ";
    if (m_params.find("boundary") != m_params.end())
      m_out << " boundary ";
    if (m_params.find("in_tissue") != m_params.end())
      m_out << " in_tissue ";  
    if (m_params.find("molecule") != m_params.end())
      m_out << " molecule ";
    m_out << endl;
  }
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    if (m_params.find("id") != m_params.end())
      m_out << format("%5d ") % p.get_id();
    if (m_params.find("tp") != m_params.end())
      m_out << format("%2d ") % p.get_type();
    if (m_params.find("flag") != m_params.end())
      m_out << format("%5d ") % p.get_flag();
    if (m_params.find("radius") != m_params.end())
      m_out << format("%8.5f ") % p.get_radius();
    if (m_params.find("coordinate") != m_params.end())
    {
      if (m_params.find("unwrap") != m_params.end())
        m_out << format(" %16.10e  %16.10e  %16.10e") % (p.x + p.ix*Lx) % (p.y + p.iy*Ly) % (p.z + p.iz*Lz);
      else
        m_out << format(" %16.10e  %16.10e  %16.10e") % p.x % p.y % p.z;
    }
    if (m_params.find("velocity") != m_params.end())
      m_out << format(" %12.7e  %12.7e  %12.7e") % p.vx % p.vy % p.vz;
    if (m_params.find("force") != m_params.end())
      m_out << format(" %12.7e  %12.7e  %12.7e") % p.fx % p.fy % p.fz;
    if (m_params.find("director") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.nx % p.ny % p.nz;
    if (m_params.find("omega") != m_params.end())
      m_out << format("%8.5f ") % p.omega;
    if (m_params.find("image_flags") != m_params.end())
      m_out << format(" %3d  %3d  %3d ") % p.ix % p.iy % p.iz;
    if (m_params.find("normal") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.Nx % p.Ny % p.Nz;
    if (m_params.find("parent") != m_params.end())
      m_out << format(" %3d ") % p.get_parent();
    if (m_params.find("area") != m_params.end())
      m_out << format("%16.6f ") % p.A0;
    if (m_nlist->has_faces())
    {
      Vertex& V = mesh.get_vertices()[p.get_id()];
      if (m_params.find("cell_area") != m_params.end())
      {
        m_out << format("%16.6f ") % V.area;
      }
      if (m_params.find("cell_perim") != m_params.end())
      {
        m_out << format("%16.6f ") % V.perim;
      }
      if (m_params.find("shape_param") != m_params.end())
      {
        m_out << format("%10.6f ") % (V.perim/sqrt(V.area));
      }
    }
    if (m_params.find("cont_num") != m_params.end())
    {
        if (m_nlist->has_contacts())
          m_out << format("%2d ") % m_nlist->get_contacts(i).size();
    }
    if (m_params.find("boundary") != m_params.end())
    {
      if (p.boundary)
        m_out << " 1 ";
      else
        m_out << " 0 ";
    }   
    if (m_params.find("in_tissue") != m_params.end())
    {
      if (p.in_tissue)
        m_out << " 1 ";
      else
        m_out << " 0 ";
    }  
    if (m_params.find("molecule") != m_params.end())
      m_out << format("%2d ") % p.molecule;
    m_out << endl;
  }
}

//! Dump format suitable for input
void Dump::dump_input()
{
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  if (m_print_header)
  {
    m_out << "# " << format(" Lx = %10.6f, Ly = %10.6f, Lz = %10.6f") % m_system->get_box()->Lx % m_system->get_box()->Ly % m_system->get_box()->Lz << endl;
    m_out << "# id type  radius x y z vx vy vz nx ny nz  omega length ix iy iz" << endl;
  }
  m_out << "keys: id molecule type radius  x   y   z   vx   vy   vz   nx   ny   nz  length ix iy iz" << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    //m_out << format("%5d  %2d  %8.4f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %3d  %3d  %3d") % p.get_id() % p.get_type() % p.get_radius() % p.x % p.y % p.z % p.vx % p.vy % p.vz % p.nx % p.ny % p.nz % p.omega % p.get_length() % p.ix % p.iy % p.iz << endl;
    m_out << format("%5d  %3d  %2d  %8.4f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %3d  %3d  %3d") % p.get_id() % p.molecule % p.get_type() % p.get_radius() % p.x % p.y % p.z % p.vx % p.vy % p.vz % p.nx % p.ny % p.nz % p.omega % p.ix % p.iy % p.iz << endl;
  }
}

//! Dumps velocities in a format suitable for visualization with GNUPlot
//! If scale is present, all velocities are rescales by a given factor.
void Dump::dump_velocity()
{
  double scale = 1.0;
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  if (m_params.find("scale") != m_params.end())
  {
    m_msg->msg(Messenger::INFO,"Scaling all velocities by "+m_params["scale"]+".");
    scale = lexical_cast<double>(m_params["scale"]);
  }
  if (m_print_header)
  {
    m_out << "# scale a = " << scale << endl;
    m_out << "# x-0.5*a*vx  y-0.5*a*vy z-0.5*a*vz  a*vx  a*vy  a*vz" << endl;
  }
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    m_out << format("%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % (p.x-0.5*scale*p.vx) % (p.y-0.5*scale*p.vy) % (p.z-0.5*scale*p.vz) % (scale*p.vx) % (scale*p.vy) % (scale*p.vz) << endl;
  }
}

//! Dumps director vectors in a format suitable for visualization with GNUPlot
//! If scale is present, all vectors are rescales by a given factor.
void Dump::dump_director()
{
  double scale = 1.0;
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  if (m_params.find("scale") != m_params.end())
  {
    m_msg->msg(Messenger::INFO,"Scaling all director vectors by "+m_params["scale"]+".");
    scale = lexical_cast<double>(m_params["scale"]);
  }
  if (m_print_header)
  {
    m_out << "# scale a = " << scale << endl;
    m_out << "# x-0.5*a*nx  y-0.5*a*ny z-0.5*a*nz  a*nx  a*ny  a*nz" << endl;
  }
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    m_out << format("%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % (p.x-0.5*scale*p.nx) % (p.y-0.5*scale*p.ny) % (p.z-0.5*scale*p.nz) % (scale*p.nx) % (scale*p.ny) % (scale*p.nz) << endl;
  }
}

//! Dump particle coordinates in the common XYZV format
//! suitable for visualization with SimRePlay toolkit
void Dump::dump_xyzv()
{
  double scale = 1.0;
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  if (m_params.find("scale") != m_params.end())
  {
#ifndef NDEBUG
    m_msg->msg(Messenger::INFO,"Scaling all director vectors by "+m_params["scale"]+".");
#endif
    scale = lexical_cast<double>(m_params["scale"]);
  }
  m_out << N << endl;
  if (m_params.find("velocity") != m_params.end())
    m_out << "Generated by SAMoS code. Vectors are velocities scaled by " << scale << endl;
  else
    m_out << "Generated by SAMoS code. Vectors are directors scaled by " << scale << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    if (m_params.find("velocity") != m_params.end())
      m_out << format("%c\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % (scale*p.vx) % (scale*p.vy) % (scale*p.vz) << endl;
    else
      m_out << format("%c\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % (scale*p.nx) % (scale*p.ny) % (scale*p.nz) << endl;
  }
}

//! Dump particle coordinates in the common XYZC format
//! suitable for visualization with SimRePlay toolkit
void Dump::dump_xyzc()
{
  int N = m_system->get_group(m_group)->get_size();
  vector<int> particles = m_system->get_group(m_group)->get_particles();
  m_system->enable_per_particle_eng();
  m_msg->msg(Messenger::WARNING,"XYZC file format output enabled per particle energy tracking. There fill be a substantial performance penalty (using slow STL maps).");
  m_out << N << endl;
  if (m_params.find("potential") != m_params.end())
    m_out << "Generated by SAMoS code. Printing potential of type " << m_params["potential"] << " for each particle." << endl;
  else if (m_params.find("aligner") != m_params.end())
    m_out << "Generated by SAMoS code. Printing alignment potential of type " << m_params["aligner"] << " for each particle." << endl;
  else
  {
    m_msg->msg(Messenger::ERROR,"Either potential or alignment have to be specified for the XYZC file format.");
    throw runtime_error("Missing print type in XYZC dump.");
  }
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(particles[i]);
    if (m_params.find("potential") != m_params.end())
      m_out << format("%c\t%10.6f\t%10.6f\t%10.6f\t%10.6f") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % p.get_pot_energy(m_params["potential"]) << endl;
    else
      m_out << format("%c\t%10.6f\t%10.6f\t%10.6f\t%10.6f") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % p.get_align_energy(m_params["aligner"]) << endl;
  }
}

//! Dump particle coordinates and bonds in the common MOL2 format
//! suitable for visualization with VMD
void Dump::dump_mol2()
{
  int N = m_system->size();
  int Nbonds = m_system->num_bonds();
  m_out << "@<TRIPOS>MOLECULE" << endl;
  m_out << "Generated by SAMoS code" << endl;
  m_out << N << "  " << Nbonds << endl;
  m_out << "NO_CHARGES" << endl;
  m_out << "@<TRIPOS>ATOM" << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    m_out << format("%d\t%d\t%10.6f\t%10.6f\t%10.6f\t%d") % (p.get_id()+1) % p.get_type() % p.x % p.y % p.z % p.get_type() << endl;
  }
  m_out << "@<TRIPOS>BOND" << endl;
  for (int i = 0; i < Nbonds; i++)
  {
    Bond& b = m_system->get_bond(i);
    m_out << format("%d\t%d\t%d\t%d") % (b.id + 1) % (b.i + 1) % (b.j + 1) % b.type << endl;
  }
}

//! Dump contact network for particles
//! Format is a simple text file with 3 columns
//! Column 1: contact id (starting with 0)
//! Column 2: id of the first particle in the contact (lower of the two ids)
//! Column 3: id of the second particle in the contact (larger of the two ids)
void Dump::dump_contact()
{
  if (!m_nlist)
  {
    m_msg->msg(Messenger::ERROR,"In order to produce contact network you need to specify neighbour list. Please use nlist command.");
    throw runtime_error("No neighbour list specified for contact network dump.");
  }
  else
  {
    double rcut;
    bool include_flag = false;
    if (m_params.find("rcut") == m_params.end())
    {
      m_msg->msg(Messenger::WARNING,"No cutoff distance for the contact network set. Setting it to 1.");
      rcut = 1.0;
    }
    else
    {
      m_msg->msg(Messenger::INFO,"Cutoff distance for contact network set to "+m_params["rcut"]+".");
      rcut = lexical_cast<double>(m_params["rcut"]);
    }
    if (rcut > m_nlist->get_cutoff())
    {
      m_msg->msg(Messenger::WARNING,"Contact network cutoff distance larger than the neighbour list cutoff. Setting it to that of the neighbour list.");
      rcut = m_nlist->get_cutoff();
    }
    if (m_params.find("include_flag") != m_params.end())
      include_flag = true;
    if (m_print_header)
    {
       m_out << "#  contact_id   id_p1  id_p2";
       if (include_flag)
         m_out << "  flag_p1   flag_p2";
       m_out << endl;
    }
    bool periodic = m_system->get_periodic();
    int N = m_system->get_group(m_group)->get_size();
    vector<int> particles = m_system->get_group(m_group)->get_particles();
    int contact = 0;
    double rcut2 = rcut*rcut;
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(particles[i]);
      vector<int>& neigh = m_nlist->get_neighbours(pi.get_id());
      for (unsigned int j = 0; j < neigh.size(); j++)
      {
        Particle& pj = m_system->get_particle(neigh[j]);
        double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        if (periodic)
          m_system->apply_periodic(dx,dy,dz);
        double r_sq = dx*dx + dy*dy + dz*dz;
        if (r_sq <= rcut2)
        {
          m_out << format("%d\t%d\t%d") % (contact++) % pi.get_id() % pj.get_id();
          if (include_flag)
            m_out << format("\t%d\t%d") % pi.get_flag() % pj.get_flag();
          m_out << endl;
        }
      }
    }
  }
}

//! Dump faces based on the contact network for particles
//! Format is a simple text file with f+1 columns, where f in the number of vertices in a face
//! Column 1: contact id (starting with 0)
void Dump::dump_faces()
{
  if (!m_nlist)
  {
    m_msg->msg(Messenger::ERROR,"In order to produce faces based on the contact network you need to specify neighbour list. Please use nlist command.");
    throw runtime_error("No neighbour list specified for faces dump.");
  }
  else if (m_nlist->has_faces())
  {
    Mesh& mesh = m_system->get_mesh();
    if (m_print_header)
       m_out << "#  face_id   partice_ids" << endl;
    vector<Face>& faces = mesh.get_faces();
    for (unsigned int i = 0; i < faces.size(); i++)
    {
      m_out << format("%d  ") % faces[i].id;
      for (unsigned int j = 0; j < faces[i].vertices.size(); j++)
        m_out << format("%d ") % faces[i].vertices[j];
      m_out << endl;
    }
  }
}

//! Dump mesh, that is dual of the particle position network 
//! This will be done in the OFF format
void Dump::dump_mesh()
{
  Mesh& m = m_system->get_mesh();
  vector<Face>& faces = m.get_faces();
  vector<Vertex>& vertices = m.get_vertices();
  int Nface = 0;
  for (unsigned int i = 0; i < vertices.size(); i++)
    if (!vertices[i].boundary) Nface++;
  m_out << "OFF" << endl;
  m_out << faces.size() << " " << Nface << " 0 " << endl;
  for (unsigned int i = 0; i < faces.size(); i++)
    m_out << format("%f\t%f\t%f") % faces[i].rc.x % faces[i].rc.y % faces[i].rc.z << endl;
  for (unsigned int i = 0; i < vertices.size(); i++)
    if (!vertices[i].boundary)
    {
      m_out << vertices[i].faces.size() << " ";
      for (unsigned int k = 0; k < vertices[i].faces.size(); k++)
        m_out << vertices[i].faces[k] << " ";
      m_out << endl;
    }
}

//! Dump vertices that are on the boundary 
//! This is primarily used to for restarting tissue simulations. 
void Dump::dump_boundary()
{
  Mesh& mesh = m_system->get_mesh();
  int i = 0;
  if (m_print_header)
    m_out << "# id  v1  v2" << endl;
  for (vector<pair<int,int> >::iterator it = mesh.get_boundary().begin(); it != mesh.get_boundary().end(); it++)
    m_out << i++ << "  " << (*it).first << "  " << (*it).second << endl;
}

//! Dump three files (vertices, faces and boundary)  
//! This is to be used to for simulations with the Active Junction Model 
void Dump::dump_ajm(int step)
{
  Mesh& mesh = m_system->get_mesh();
  if (m_params.find("json") == m_params.end())
  {
    string vert_file_name = m_directory+"/"+m_file_name+"_vert_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+"."+m_ext;
    string cell_file_name = m_directory+"/"+m_file_name+"_cell_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+"."+m_ext;
    string bnd_file_name = m_directory+"/"+m_file_name+"_boundary_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+"."+m_ext;
    string tp_file_name = m_directory+"/"+m_file_name+"_types_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+"."+m_ext;
    ofstream vert_file, cell_file, bnd_file, tp_file;

    vert_file.open(vert_file_name.c_str());
    cell_file.open(cell_file_name.c_str());
    bnd_file.open(bnd_file_name.c_str());
    tp_file.open(tp_file_name.c_str());

    PlotArea& pa = mesh.plot_area(false);
    for (unsigned int f = 0; f < pa.points.size(); f++)
      vert_file << format("%5d   %.7f   %.7f") % f % pa.points[f].x % pa.points[f].y << endl;
    
    for (unsigned int f = 0; f < pa.sides.size(); f++)
    {
      for (unsigned int d = 0; d < pa.sides[f].size(); d++)
        cell_file << format("%5d  ") % pa.sides[f][d];
      cell_file << endl;
    }
    
    for (unsigned int f = 0; f < pa.type.size(); f++)
      tp_file << format("%2d") % (pa.type[f]-1) << endl;

    tp_file << endl;

    for (unsigned int f = 0; f < pa.boundary_faces.size(); f++)
    {
      cell_file << format("%5d  ") % pa.boundary_faces[f];  
      bnd_file << format("%5d") % pa.boundary_faces[f] << endl;
    }
    cell_file << endl;
    
    vert_file.close();
    cell_file.close();
    bnd_file.close();
    tp_file.close();
  }
  else
  {
    pt::ptree json; 
    string json_file_name = m_directory+"/"+m_file_name+"_"+lexical_cast<string>(format("%010d") % (step+m_time_step_offset))+".json";  
    PlotArea& pa = mesh.plot_area(false);

    // collect vertices
    pt::ptree vertices;
    for (unsigned int f = 0; f < pa.points.size(); f++)    
    {
      pt::ptree vertex;
      vertex.put("id", f);
      pt::ptree r; 
      pt::ptree pos;
      pos.put("",pa.points[f].x);
      r.push_back(std::make_pair("",pos));
      pos.put("",pa.points[f].y);
      r.push_back(std::make_pair("",pos));
      vertex.add_child("r", r);
      if (find(pa.boundary_faces.begin(),pa.boundary_faces.end(), f) == pa.boundary_faces.end())
        vertex.put("boundary", "false");
      else
        vertex.put("boundary", "true");
      vertex.put("constraint", "none");
      vertex.put("type", "regular");
      vertices.push_back(std::make_pair("", vertex));
    }
    json.add_child("mesh.vertices", vertices);
    // collect cells
    pt::ptree cells;
    for (unsigned int f = 0; f < pa.sides.size(); f++)
    {
      pt::ptree cell;
      cell.put("id",f);
      if (m_params.find("A0") != m_params.end())
        cell.put("A0", lexical_cast<double>(m_params["A0"]));
      else
        cell.put("A0", pa.area[f]);  
      if (m_params.find("P0") != m_params.end())
        cell.put("P0", lexical_cast<double>(m_params["P0"]));
      else
        cell.put("P0", pa.perim[f]);  
      pt::ptree vert; 
      for (unsigned int d = 0; d < pa.sides[f].size(); d++)
      {
        pt::ptree val; 
        val.put("", pa.sides[f][d]);
        vert.push_back(std::make_pair("", val));
      }
      cell.add_child("vertices",vert);
      cell.put("type", pa.type[f]);
      cell.put("outer", "false");
      cell.put("nsides", pa.sides[f].size());
      cells.push_back(std::make_pair("", cell));
    }
    // collect the outer face
    pt::ptree outer;
    outer.put("id",pa.sides.size());
    outer.put("A0", 0.0);
    outer.put("P0", 0.0);
    pt::ptree vert;
    for (unsigned int d = 0; d < pa.boundary_faces.size(); d++)
    {
      pt::ptree val; 
      val.put("", pa.boundary_faces[d]);
      vert.push_back(std::make_pair("", val));
    }
    outer.add_child("vertices",vert);
    outer.put("outer", "true");
    outer.put("type", "passive");
    outer.put("nsides", pa.boundary_faces.size());
    cells.push_back(std::make_pair("", outer));
    json.add_child("mesh.faces", cells);

    // Write everything to the file
    pt::write_json(json_file_name, json);
  }
}

#ifdef HAS_VTK
//! Dump meshes into VTK output 
void Dump::dump_vtp(int step)
{
  vector<pair<int,int> > visited_edges;
  string file_name = m_directory+"/"+m_file_name+"_"+lexical_cast<string>(format("%010d") % step)+"."+m_ext;
  vtkSmartPointer<vtkPolyData> polydata =  vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> lines =  vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkCellArray> faces = vtkSmartPointer<vtkCellArray>::New();
  map<int,int> index_map;   // Maps particle indices to the current index (for different groups)
  Mesh& mesh = m_system->get_mesh();
  
  if (!m_output_dual)
  {
    int N = m_system->get_group(m_group)->get_size();
    vector<int> particles = m_system->get_group(m_group)->get_particles();
    
    vtkSmartPointer<vtkIntArray> ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkIntArray> types =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> radii =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> a0 =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> press =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> vel =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> force =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> dir =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> ndir =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> dual_area =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> num_neigh =  vtkSmartPointer<vtkDoubleArray>::New();
    map<string,vtkSmartPointer<vtkDoubleArray> > force_components;
    
    ids->SetName("Id");
    ids->SetNumberOfComponents(1);
    types->SetName("Type");
    types->SetNumberOfComponents(1);
    radii->SetName("Radius");
    radii->SetNumberOfComponents(1);
    a0->SetName("NativeArea");
    a0->SetNumberOfComponents(1);
    press->SetName("Pressure");
    press->SetNumberOfComponents(1);
    vel->SetName("Velocity");
    vel->SetNumberOfComponents(3);
    force->SetName("Force");
    force->SetNumberOfComponents(3);
    dir->SetName("Director");
    dir->SetNumberOfComponents(3);
    ndir->SetName("NDirector");
    ndir->SetNumberOfComponents(3);
    dual_area->SetName("DualArea");
    dual_area->SetNumberOfComponents(1);
    num_neigh->SetName("NumNeigh");
    num_neigh->SetNumberOfComponents(1);

    if (m_system->record_force_type())
    {
      Particle& pi = m_system->get_particle(particles[0]);
      for (map<string,ForceType>::iterator it_t = pi.get_force_type().begin(); it_t != pi.get_force_type().end(); it_t++)
      {
        force_components[(*it_t).first] = vtkSmartPointer<vtkDoubleArray>::New();
        force_components[(*it_t).first]->SetName((*it_t).first.c_str());
        force_components[(*it_t).first]->SetNumberOfComponents(3);
      }
    }
      
    for (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(particles[i]);
      double v[3] = {pi.vx, pi.vy, pi.vz};
      double f[3] = {pi.fx, pi.fy, pi.fz};
      double n[3] = {pi.nx, pi.ny, pi.nz};
      double nn[3] = {-pi.nx, -pi.ny, -pi.nz};
      index_map[pi.get_id()] = i;
      points->InsertNextPoint ( pi.x, pi.y, pi.z );
      ids->InsertNextValue(pi.get_id());
      types->InsertNextValue(pi.get_type());
      radii->InsertNextValue(pi.get_radius());
      a0->InsertNextValue(pi.A0);
      press->InsertNextValue(pi.s_xx + pi.s_yy + pi.s_zz);
      vel->InsertNextTuple(v);
      force->InsertNextTuple(f);
      dir->InsertNextTuple(n);
      ndir->InsertNextTuple(nn);
      if (m_system->record_force_type())
      {
        for (map<string,ForceType>::iterator it_t = pi.get_force_type().begin(); it_t != pi.get_force_type().end(); it_t++)
        {
          double fc[3] = {(*it_t).second.fx, (*it_t).second.fy, (*it_t).second.fz};
          force_components[(*it_t).first]->InsertNextTuple(fc);
        }
      }
      if (mesh.size() > 0)
      {
        Vertex& V = mesh.get_vertices()[i];
        // This is just a workaround a bug in some versions of Paraview that do not properly parse numbers in scientific notation
        if (V.area > 1e-7)  dual_area->InsertNextValue(V.area);
         else  dual_area->InsertNextValue(0.0);
        num_neigh->InsertNextValue(V.n_edges);
      }
      else
        num_neigh->InsertNextValue(pi.coordination);

    }
    
    polydata->SetPoints(points);
    
    polydata->GetPointData()->AddArray(ids);
    polydata->GetPointData()->AddArray(types);
    polydata->GetPointData()->AddArray(radii);
    polydata->GetPointData()->AddArray(a0);
    polydata->GetPointData()->AddArray(press);
    polydata->GetPointData()->AddArray(vel);
    polydata->GetPointData()->AddArray(force);
    polydata->GetPointData()->AddArray(dir);
    polydata->GetPointData()->AddArray(ndir);
    polydata->GetPointData()->AddArray(num_neigh); // A header for this with no data gets put into the .vtp file even when there is no mesh
    if (mesh.size() > 0)
    {
      polydata->GetPointData()->AddArray(dual_area);
    }
    if (m_system->record_force_type())
    {
      Particle& pi = m_system->get_particle(particles[0]);
      for (map<string,ForceType>::iterator it_t = pi.get_force_type().begin(); it_t != pi.get_force_type().end(); it_t++)
        polydata->GetPointData()->AddArray(force_components[(*it_t).first]);
    }

        
    if (m_system->num_bonds() > 0 && m_include_bonds)
    {
      vtkSmartPointer<vtkLine> line =  vtkSmartPointer<vtkLine>::New();
      vtkSmartPointer<vtkDoubleArray> lens =  vtkSmartPointer<vtkDoubleArray>::New();
      lens->SetName("Length");
      lens->SetNumberOfComponents(1);
      for (int b = 0; b < m_system->num_bonds(); b++)
      {
        Bond& bond = m_system->get_bond(b);
        line->GetPointIds()->SetId(0, index_map[bond.i]); 
        line->GetPointIds()->SetId(1, index_map[bond.j]);
        lines->InsertNextCell(line);
        Particle& pi = m_system->get_particle(bond.i);
        Particle& pj = m_system->get_particle(bond.j);
        double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
        m_system->apply_periodic(dx,dy,dz);
        lens->InsertNextValue(sqrt(dx*dx + dy*dy + dz*dz));
      }
      polydata->SetLines(lines);
      polydata->GetCellData()->AddArray(lens);
    }
    
    if (mesh.size() > 0 && m_include_mesh)
    {
      vtkSmartPointer<vtkLine> edge =  vtkSmartPointer<vtkLine>::New();
      vtkSmartPointer<vtkDoubleArray> lens =  vtkSmartPointer<vtkDoubleArray>::New();
      vtkSmartPointer<vtkDoubleArray> boundary =  vtkSmartPointer<vtkDoubleArray>::New();
      vtkSmartPointer<vtkDoubleArray> boundary_part =  vtkSmartPointer<vtkDoubleArray>::New();
      vtkSmartPointer<vtkIntArray> boundary_edges =  vtkSmartPointer<vtkIntArray>::New();
      lens->SetName("Length");
      lens->SetNumberOfComponents(1);
      boundary->SetName("Boundary");
      boundary->SetNumberOfComponents(1);
      boundary_part->SetName("BoundaryPart");
      boundary_part->SetNumberOfComponents(1);
      boundary_edges->SetName("Boundary");
      boundary_edges->SetNumberOfComponents(1);
      // Edited for this to be backwards for better paraview plotting
      for (int i = 0; i < mesh.size(); i++)
      {
        Vertex& vert = mesh.get_vertices()[i];
        if (vert.boundary && vert.attached)
          boundary->InsertNextValue(0);
        else if (!vert.attached)
          boundary->InsertNextValue(1);
        else
          boundary->InsertNextValue(2);
      }
      polydata->GetPointData()->AddArray(boundary);
      for (int i = 0; i < N; i++)
      {
        Particle& pi = m_system->get_particle(particles[i]);
        if (pi.boundary)
          boundary_part->InsertNextValue(0);
        else if (!pi.boundary && pi.in_tissue)
          boundary_part->InsertNextValue(1);
        else
          boundary_part->InsertNextValue(2);    // non-tissue (environment) particles have flag 2
      }
      polydata->GetPointData()->AddArray(boundary_part);
      for (int e = 0; e < mesh.nedges(); e++)
      {
        Edge& ee = mesh.get_edges()[e];
        if ( find(visited_edges.begin(),visited_edges.end(),make_pair(ee.from,ee.to)) == visited_edges.end() )
        {
          edge->GetPointIds()->SetId(0, index_map[ee.from]); 
          edge->GetPointIds()->SetId(1, index_map[ee.to]);
          lines->InsertNextCell(edge);
          Particle& pi = m_system->get_particle(ee.from);
          Particle& pj = m_system->get_particle(ee.to);
          double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
          m_system->apply_periodic(dx,dy,dz);
          lens->InsertNextValue(sqrt(dx*dx + dy*dy + dz*dz));
          if (ee.boundary || mesh.get_edges()[ee.pair].boundary)
            boundary_edges->InsertNextValue(0);
          else
            boundary_edges->InsertNextValue(2);
          visited_edges.push_back(make_pair(ee.from,ee.to));
          visited_edges.push_back(make_pair(ee.to,ee.from));
        }
      }
      polydata->SetLines(lines);
      polydata->GetCellData()->AddArray(lens);
      polydata->GetCellData()->AddArray(boundary_edges);
      
      vtkSmartPointer<vtkPolygon> face =  vtkSmartPointer<vtkPolygon>::New();
      for (int f = 0; f < mesh.nfaces(); f++)
      {
        Face& ff = mesh.get_faces()[f];
        face->GetPointIds()->SetNumberOfIds(ff.n_sides);
        for (int fi = 0; fi < ff.n_sides; fi++)
          face->GetPointIds()->SetId(fi, index_map[ff.vertices[fi]]);
        faces->InsertNextCell(face);
      }
    }
  }
  else
  {
    vtkSmartPointer<vtkPolygon> face =  vtkSmartPointer<vtkPolygon>::New();
    vtkSmartPointer<vtkIntArray> ids =  vtkSmartPointer<vtkIntArray>::New();
    vtkSmartPointer<vtkDoubleArray> areas =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> perims =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> p0 =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkDoubleArray> circum_radius =  vtkSmartPointer<vtkDoubleArray>::New();
    vtkSmartPointer<vtkIntArray> types =  vtkSmartPointer<vtkIntArray>::New();
    ids->SetName("Id");
    areas->SetNumberOfComponents(1);
    areas->SetName("Area");
    areas->SetNumberOfComponents(1);
    perims->SetName("Perimeter");
    perims->SetNumberOfComponents(1);
    circum_radius->SetName("CircumRadius");
    circum_radius->SetNumberOfComponents(1);
    p0->SetName("p0");
    p0->SetNumberOfComponents(1);
    types->SetName("type");
    
    PlotArea& pa = mesh.plot_area(m_dual_boundary);
    for (unsigned int f = 0; f < pa.points.size(); f++)
    {
      points->InsertNextPoint(pa.points[f].x, pa.points[f].y, pa.points[f].z);
      circum_radius->InsertNextValue(pa.circum_radius[f]);
      ids->InsertNextValue(f);
    }
    polydata->SetPoints(points);
    polydata->GetPointData()->AddArray(ids);
    polydata->GetPointData()->AddArray(circum_radius);

    for (unsigned int f = 0; f < pa.sides.size(); f++)
    {
      face->GetPointIds()->SetNumberOfIds(pa.sides[f].size());
      for (unsigned int d = 0; d < pa.sides[f].size(); d++)
        face->GetPointIds()->SetId(d, pa.sides[f][d]);
      faces->InsertNextCell(face);
      areas->InsertNextValue(pa.area[f]);
      perims->InsertNextValue(pa.perim[f]);
      p0->InsertNextValue(pa.perim[f]/sqrt(pa.area[f]));
      types->InsertNextValue(pa.type[f]);
    }
    
    polydata->SetPolys(faces);
    polydata->GetCellData()->AddArray(ids);
    polydata->GetCellData()->AddArray(areas);
    polydata->GetCellData()->AddArray(perims);
    polydata->GetCellData()->AddArray(p0);
    polydata->GetCellData()->AddArray(types);
  }
  
  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(file_name.c_str());
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif
  writer->SetDataModeToAscii();
  writer->Write();
}
#endif

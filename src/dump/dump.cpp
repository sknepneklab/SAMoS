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
 *  \param fname Base file name for the output 
 *  \param params list of parameters that control dump (e.g., frequency, output type, etc.)
*/
Dump::Dump(SystemPtr sys, MessengerPtr msg, const string& fname, pairs_type& params) : m_system(sys), m_msg(msg), m_file_name(fname), m_params(params)
{
  m_type_ext["velocity"] = "vel";
  m_type_ext["xyz"] = "xyz";
  m_type_ext["full"] = "dat";
  m_type_ext["input"] = "inp";
  m_type_ext["dcd"] = "dcd";
  m_type_ext["director"] = "dir";
  m_type_ext["xyzv"] = "xyzv";

  m_print_header = false;
  
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
  if (params.find("multi") == params.end())
  {
    m_msg->msg(Messenger::WARNING,"All time steps will be concatenated to a single file.");
    m_multi_print = false;
    string file_name = m_file_name+"."+m_ext;
    m_out.open(file_name.c_str()); 
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
  }
  if (m_type == "dcd")
  {
    if (m_multi_print)
    {
      m_msg->msg(Messenger::ERROR,"DCD files cannot be split. mulifile flag is not permitted.");
      runtime_error("Attempted to spread DCD output into multiple files.");
    }
    else
    {
      write_int(m_out, 84);
      
      // the next 4 bytes in the file must be "CORD"
      char cord_data[] = "CORD";
      m_out.write(cord_data, 4);
      write_int(m_out, 0);      // Number of frames in file, none written yet
      write_int(m_out, m_start); // Starting timestep
      write_int(m_out, m_freq);  // Timesteps between frames written to the file
      write_int(m_out, 0);      // Number of timesteps in simulation
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);                  // timestep (unused)
      write_int(m_out, 1);                  // include unit cell
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 0);
      write_int(m_out, 24); // Pretend to be CHARMM version 24
      write_int(m_out, 84);
      write_int(m_out, 164);
      write_int(m_out, 2);
    
      char title_string[81];
      char remarks[] = "Created by APCS";
      strncpy(title_string, remarks, 80);
      title_string[79] = '\0';
      m_out.write(title_string, 80);
  
      char time_str[81];
      time_t cur_time = time(NULL);
      tm *tmbuf=localtime(&cur_time);
      strftime(time_str, 80, "REMARKS Created  %d %B, %Y at %H:%M", tmbuf);
      m_out.write(time_str, 80);
      
      write_int(m_out, 164);
      write_int(m_out, 4);
      write_int(m_out, m_system->size());
      write_int(m_out, 4);
      
      // check for errors
      if (!m_out.good())
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
    string file_name = m_file_name+"_"+lexical_cast<string>(format("%010d") % step)+"."+m_ext;
    m_out.open(file_name.c_str());
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
  
  if (m_multi_print)
    m_out.close();
}

// Private functions 

//! Dump coordinated in a DCD file
void Dump::dump_dcd()
{
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
  write_int(m_out, 48);
  m_out.write((char *)unitcell, 48);
  write_int(m_out, 48);
  // check for errors
  if (!m_out.good())
  {
    m_msg->msg(Messenger::ERROR,"Error writing into DCD file.");
    throw runtime_error("Error writing DCD file");
  } 
  // Prepare x coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    buffer[i] = static_cast<float>(p.x);
  }
  // write x coords
  write_int(m_out, N * sizeof(float));
  m_out.write((char *) buffer, N * sizeof(float));
  write_int(m_out, N * sizeof(float));
  
  // Prepare y coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    buffer[i] = static_cast<float>(p.y);
  }
  // write y coords
  write_int(m_out, N * sizeof(float));
  m_out.write((char *) buffer, N * sizeof(float));
  write_int(m_out, N * sizeof(float));
  
  // Prepare z coordinate
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    buffer[i] = static_cast<float>(p.z);
  }
  // write z coords
  write_int(m_out, N * sizeof(float));
  m_out.write((char *) buffer, N * sizeof(float));
  write_int(m_out, N * sizeof(float));
  
  delete [] buffer;
}

//! Dump particle coordinates in the common XYZ format
void Dump::dump_xyz()
{
  int N = m_system->size();
  m_out << N << endl;
  m_out << "Generated by APCS code." << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    m_out << format("%3d\t%10.6f\t%10.6f\t%10.6f") % p.get_type() % p.x % p.y % p.z << endl;
  }
}

//! Dump selected set of data
void Dump::dump_data()
{
  int N = m_system->size();
  if (m_print_header)
  {
    m_out << "# ";
    if (m_params.find("id") != m_params.end())
      m_out << " id ";
    if (m_params.find("tp") != m_params.end())
      m_out << " type ";
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
    m_out << endl;
  }
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_params.find("id") != m_params.end())
      m_out << format("%5d ") % p.get_id();
    if (m_params.find("tp") != m_params.end())
      m_out << format("%2d ") % p.get_type();
    if (m_params.find("radius") != m_params.end())
      m_out << format("%8.5f ") % p.get_radius();
    if (m_params.find("coordinate") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.x % p.y % p.z;
    if (m_params.find("velocity") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.vx % p.vy % p.vz;
    if (m_params.find("force") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.fx % p.fy % p.fz;
    if (m_params.find("director") != m_params.end())
      m_out << format(" %10.6f  %10.6f  %10.6f") % p.nx % p.ny % p.nz;
    if (m_params.find("omega") != m_params.end())
      m_out << format("%8.5f ") % p.omega;
    m_out << endl;
  }
}

//! Dump format suitable for input
void Dump::dump_input()
{
  int N = m_system->size();
  if (m_print_header)
    m_out << "# id type  radius x y z vx vy vz nx ny nz  omega" << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    m_out << format("%5d  %2d  %8.4f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f") % p.get_id() % p.get_type() % p.get_radius() % p.x % p.y % p.z % p.vx % p.vy % p.vz % p.nx % p.ny % p.nz % p.omega << endl;
  }
}

//! Dumps velocities in a format suitable for visualization with GNUPlot
//! If scale is present, all velocities are rescales by a given factor.
void Dump::dump_velocity()
{
  double scale = 1.0;
  int N = m_system->size();
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
    Particle& p = m_system->get_particle(i);
    m_out << format("%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % (p.x-0.5*scale*p.vx) % (p.y-0.5*scale*p.vy) % (p.z-0.5*scale*p.vz) % (scale*p.vx) % (scale*p.vy) % (scale*p.vz) << endl;
  }
}

//! Dumps director vectors in a format suitable for visualization with GNUPlot
//! If scale is present, all vectors are rescales by a given factor.
void Dump::dump_director()
{
  double scale = 1.0;
  int N = m_system->size();
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
    Particle& p = m_system->get_particle(i);
    m_out << format("%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e\t%12.7e") % (p.x-0.5*scale*p.nx) % (p.y-0.5*scale*p.ny) % (p.z-0.5*scale*p.nz) % (scale*p.nx) % (scale*p.ny) % (scale*p.nz) << endl;
  }
}

//! Dump particle coordinates in the common XYZV format
//! suitable for visualization with SimRePlay toolkit
void Dump::dump_xyzv()
{
  double scale = 1.0;
  int N = m_system->size();
  if (m_params.find("scale") != m_params.end())
  {
    m_msg->msg(Messenger::INFO,"Scaling all director vectors by "+m_params["scale"]+".");
    scale = lexical_cast<double>(m_params["scale"]);
  }
  m_out << N << endl;
  if (m_params.find("velocity") != m_params.end())
    m_out << "Generated by APCS code. Vectors are velocities scaled by " << scale << endl;
  else
    m_out << "Generated by APCS code. Vectors are directors scaled by " << scale << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (m_params.find("velocity") != m_params.end())
      m_out << format("%c\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % (scale*p.vx) % (scale*p.vy) % (scale*p.vz) << endl;
    else
      m_out << format("%c\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f") % static_cast<char>(64+p.get_type()) % p.x % p.y % p.z % (scale*p.nx) % (scale*p.ny) % (scale*p.nz) << endl;
  }
}
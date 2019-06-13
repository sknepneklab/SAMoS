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
 * \file neighbour_list.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Declaration of NeighbourList class.
 */ 

#ifndef __NEIGHBOUR_LIST_HPP__
#define __NEIGHBOUR_LIST_HPP__

#include <vector>
#include <stdexcept>
#include <utility>
#include <string>
#include <fstream>
#include <list>
#include <memory>


#ifdef HAS_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#endif

#include "messenger.hpp"
#include "system.hpp"
#include "cell_list.hpp"
#include "parse_parameters.hpp"

using std::vector;
using std::pair;
using std::string;
using std::ofstream;
using std::list;
using std::shared_ptr;

#ifdef HAS_CGAL

//! Auxiliary data structure for CGAL-enebled cleanup of the vertices outside the boundary.
struct FaceInfo2
{
  FaceInfo2(){}
  bool in_domain()
  { 
    return nesting_level%2 == 1;
  }
  int nesting_level;
};

/*! Typdefs for CGAL library */
typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,Kernel>       Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel,Fbb>           Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>                       TDS;
typedef CGAL::Exact_predicates_tag                                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>     Delaunay;
typedef Delaunay::Point                                                   Point;

//typedef CGAL::Triangulation_data_structure_2<Vb>                          Tds;
//typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                       Delaunay;
//typedef Delaunay::Vertex_circulator                                       Vertex_circulator;
//typedef Kernel::Point_2                                                   Point;
using std::make_pair;
using std::pair;
#endif

/*! Auxiliary structure to keep coordinates of the first state 
 *  after the build so we can check if it is necessary to rebuild.
*/
struct PartPos
{
  //! Constructor
  //! \param x x coordinate of the stored position 
  //! \param y y coordinate of the stored position 
  //! \param z z coordinate of the stored position 
  PartPos(double x, double y, double z) : x(x), y(y), z(z) { }
  //@{
  double x, y, z;  //!< Stores actual position 
  //@}
};


/*! This class handles neighbour lists for fast potential and force 
 *  calculations. It is implemented as a vector of vectors (for performance). 
 *  Since the system is on a curved surface
 *  it would be hard to make a general cell list. Therefore, we relay on
 *  the slower but generic N^2 list generation 
*/
class NeighbourList
{
public:
  
  //! Construct NeighbourList object
  //! \param sys Reference to the System object
  //! \param msg Constant reference to the Messenger object
  //! \param cutoff Cutoff distance (should be set to potential cutoff distance + padding distance)
  //! \param pad Padding distance
  NeighbourList(SystemPtr sys, MessengerPtr msg, double cutoff, double pad, pairs_type& param) : m_system(sys), 
                                                                                                 m_msg(msg),
                                                                                                 m_cut(cutoff), 
                                                                                                 m_pad(pad), 
                                                                                                 m_triangulation(false),
                                                                                                 m_max_perim(20.0),
                                                                                                 m_circumcenter(true),
                                                                                                 m_disable_nlist(false),
                                                                                                 m_remove_detached(true),
                                                                                                 m_static_boundary(false)
  {
    m_msg->write_config("nlist.cut",lexical_cast<string>(m_cut));
    m_msg->write_config("nlist.pad",lexical_cast<string>(m_pad));
    // Check if box is large enough for cell list
    if (m_system->get_box()->Lx > 2.0*(cutoff+pad) && m_system->get_box()->Ly > 2.0*(cutoff+pad) && m_system->get_box()->Lz > 2.0*(cutoff+pad))
    {
      m_use_cell_list = true;
      m_cell_list = shared_ptr<CellList>(new CellList(m_system,m_msg,cutoff+pad));
      m_msg->msg(Messenger::INFO,"Using cell lists for neighbour list builds.");
      m_msg->write_config("nlist.build_type","cell");
    }
    else
    {
      m_use_cell_list = false;
      m_msg->msg(Messenger::INFO,"Box dimensions are too small to be able to use cell lists. Neighbour list will be built using N^2 algorithm.");
      m_msg->write_config("nlist.build_type","n_square");
    }
    if (param.find("triangulation") != param.end())
    {
      if (m_system->get_periodic())
      {
        m_msg->msg(Messenger::ERROR,"Delaunay triangulation is not supported for periodic systems.");
        throw runtime_error("Delaunay triangulation not supported in periodic systems.");
      }
      m_triangulation = true;
      m_msg->msg(Messenger::INFO,"Faces will be build using Delaunay triangulation.");
      m_msg->write_config("nlist.triangulation","true"); 
    }
    if (param.find("max_perimeter") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list. No maximum face perimeter set. Assuming default value of 20.");
      m_msg->write_config("nlist.max_perimeter","20.0");
      m_max_perim = 20.0;
    }
    else    
    {
      m_msg->msg(Messenger::INFO,"Neighbour list.  Setting maximum face perimeter to "+param["max_perimeter"]+".");
      m_msg->write_config("nlist.max_perimeter",param["max_perimeter"]);
      m_max_perim =  lexical_cast<double>(param["max_perimeter"]);
    }
    if (param.find("circumcenter") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list. Using circumcenters for mesh dual.");
      m_msg->write_config("nlist.circumcenter","true");
      m_circumcenter = true;
    }
    else    
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Using geometric centers for mesh dual.");
      m_msg->write_config("nlist.circumcenter","false");
      m_circumcenter = false;
    }
    if (param.find("disable_nlist") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list will not be built. This should be used in pair with triangulations when particle connections are build based on tringulations.");
      m_msg->write_config("nlist.disable_nlist","true");
      m_disable_nlist = true;
    }
    if (param.find("keep_detached") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list. Particles detached from the mesh (tissue) will be kept.");
      m_msg->write_config("nlist.remove_detached","false");
      m_remove_detached = false;
    }
    if (param.find("max_iter") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Setting maximum number of iterations for boundary builds in tissue simulations to "+param["max_iter"]+".");
      m_msg->write_config("nlist.max_iter",param["max_iter"]);
      m_system->set_max_mesh_iterations(lexical_cast<int>(param["max_iter"]));
    }
    if (param.find("boundary_type") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Setting type of boundary particles in tissue simulations to "+param["boundary_type"]+".");
      m_msg->write_config("nlist.boundary_type",param["boundary_type"]);
      m_system->set_boundary_type(lexical_cast<int>(param["boundary_type"]));
    }
    if (param.find("static_boundary") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Using static boundaries in the tissue simulation.");
      m_msg->write_config("nlist.static_boundary","true");
      m_static_boundary = true;
    }
    this->build();
  }
  
  //! Destructor
  ~NeighbourList()
  {
    for(unsigned int i = 0; i < m_list.size(); i++)
      m_list[i].clear();
    m_list.clear();
    m_old_state.clear();
  }
  
  //! Check is neighbour list of the given particle needs update
  bool need_update(Particle&);
  
  //! Returns true is faces list exists
  bool has_faces() { return m_triangulation; }
  
  //! Get neighbour list for a give particle
  //! \param id Particle id
  //! \return Reference to the particle's neighbour list
  vector<int>& get_neighbours(int id) { return m_list[id]; }
  
  //! Get contacts for a given particle
  //! \param id Particle id
  //! \return Reference to the particle's contact list
  vector<int>& get_contacts(int id) { return m_contact_list[id]; }
  
  //! Check if contact list exists
  bool has_contacts() { return (m_contact_list.size() > 0); }
  
  //! Get neighbour list cutoff distance
  double get_cutoff() { return m_cut;  }  //!< \return neighbour list cutoff distance
  
  //! Rescales neigbour list cutoff
  //! \param scale scale factor
  void rescale_cutoff(double scale)
  {
    m_cut *= scale;
    if (m_use_cell_list && m_system->get_box()->Lx > 2.0*(m_cut+m_pad) && m_system->get_box()->Ly > 2.0*(m_cut+m_pad) && m_system->get_box()->Lz > 2.0*(m_cut+m_pad))
    {
      m_cell_list = shared_ptr<CellList>(new CellList(m_system,m_msg,m_cut+m_pad));
      m_msg->msg(Messenger::INFO,"Rescaling neighbour list cutoff.");
      m_msg->msg(Messenger::INFO,"Still using cell lists for neighbour list builds.");
    }
    else
    {
      m_use_cell_list = false;
      m_msg->msg(Messenger::INFO,"Rescaling neighbour list cutoff.");
      m_msg->msg(Messenger::INFO,"No longer possible to use cell lists for neighbour list builds.");
    }
    this->build();
  }
  
  
  //! Build neighbour list
  void build();
  
  // Does actual contact and face building 
  void build_mesh();
    
  
private:
  
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  CellListPtr m_cell_list;         //!< Pointer to the cell list
  vector<vector<int> >  m_list;    //!< Holds the list for each particle 
  vector<PartPos> m_old_state;     //!< Coordinates of particles right after the build
  double m_cut;                    //!< List build cutoff distance 
  double m_pad;                    //!< Padding distance (m_cut should be set to potential cutoff + m_pad)
  bool m_use_cell_list;            //!< If true, use cell list to speed up neighbour list builds
  bool m_triangulation;            //!< If true, build Delaunay triangulation for faces
  double m_max_perim;              //!< Maximum value of the perimeter beyond which face becomes a hole.
  bool m_circumcenter;             //!< If true, use cell circumcenters when computing duals. 
  bool m_disable_nlist;            //!< If true, neigbour list is not built (only used for cell simulations)
  bool m_remove_detached;          //!< If true, remove detached particles (vertices) before rebuilding neighbour list (for cell simulations)
  bool m_static_boundary;          //!< If true, treat tissue boundary as static, i.e., do not add new boundary particles 
  vector<vector<int> >  m_contact_list;    //!< Holds the contact list for each particle
    
  // Actual neighbour list builds
  void build_nsq(int);    //!< Build with N^2 algorithm
  void build_cell();      //!< Build using cells list
  
   //! Build faces
  void build_faces(bool);

  // Remove dangling edges
  void remove_dangling();
  
  // Remove detached particles
  void remove_detached();
 
  
#ifdef HAS_CGAL
  // Build Delaunay triangulation
  bool build_triangulation();
#endif
    
  //! Dump particles and connectivity into a MOL2 file for debugging
  void debug_dump(const string&);

};

typedef shared_ptr<NeighbourList> NeighbourListPtr;

//! Dot product between two vectors given by three particles
double dot(const Particle&, const Particle&, const Particle&);

//! Mirror image of a point with respect to an edge
void mirror(const Particle&, const Particle&, const Particle&, double&, double&, double&);


#endif

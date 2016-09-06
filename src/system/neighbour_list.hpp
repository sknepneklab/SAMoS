/* *************************************************************
 *  
 *   Soft Active Mater on Surfaces (SAMoS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013, 2014
 * 
 *   School of Science and Engineering
 *   School of Life Sciences 
 *   University of Dundee
 * 
 *   (c) 2015
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

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

//#include <boost/property_map/property_map.hpp>
//#include <boost/ref.hpp>

#ifdef HAS_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#endif

#include "messenger.hpp"
#include "system.hpp"
#include "cell_list.hpp"
#include "parse_parameters.hpp"

using std::vector;
using std::pair;
using std::string;
using std::ofstream;

#ifdef HAS_CGAL
/*! Typdefs for CGAL library */
typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                          Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>                       Delaunay;
typedef Delaunay::Vertex_circulator                                       Vertex_circulator;
typedef Kernel::Point_2                                                   Point;
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
                                                                                                 m_build_contacts(false),
                                                                                                 m_build_faces(false),
                                                                                                 m_triangulation(false),
                                                                                                 m_contact_dist(0.0),
                                                                                                 m_max_perim(20.0),
                                                                                                 m_max_edge_len(7.0),
                                                                                                 m_circumcenter(true),
                                                                                                 m_disable_nlist(false),
                                                                                                 m_remove_detached(false)
  {
    m_msg->write_config("nlist.cut",lexical_cast<string>(m_cut));
    m_msg->write_config("nlist.pad",lexical_cast<string>(m_pad));
    // Check if box is large enough for cell list
    if (m_system->get_box()->Lx > 2.0*(cutoff+pad) && m_system->get_box()->Ly > 2.0*(cutoff+pad) && m_system->get_box()->Lz > 2.0*(cutoff+pad))
    {
      m_use_cell_list = true;
      m_cell_list = boost::shared_ptr<CellList>(new CellList(m_system,m_msg,cutoff+pad));
      m_msg->msg(Messenger::INFO,"Using cell lists for neighbour list builds.");
      m_msg->write_config("nlist.build_type","cell");
    }
    else
    {
      m_use_cell_list = false;
      m_msg->msg(Messenger::INFO,"Box dimensions are too small to be able to use cell lists. Neighbour list will be built using N^2 algorithm.");
      m_msg->write_config("nlist.build_type","n_square");
    }
    if (param.find("build_contacts") != param.end())
    {
      m_build_contacts = true;
      m_msg->msg(Messenger::INFO,"Neighbour list will also build contact network.");
      m_msg->write_config("nlist.contact_network","true"); 
    }
    if (param.find("triangulation") != param.end())
    {
      m_triangulation = true;
      m_msg->msg(Messenger::INFO,"Faces will be build using Delaunay triangulation.");
      m_msg->write_config("nlist.triangulation","true"); 
      if (param.find("max_edge_len") == param.end())
      {
        m_msg->msg(Messenger::WARNING,"Neighbour list. No maximum edge lenght set. Assuming default value of 7.");
        m_msg->write_config("nlist.max_edge_len","7.0");
        m_max_edge_len = 7.0;
      }
      else    
      {
        m_msg->msg(Messenger::INFO,"Neighbour list. Setting maximum edge length "+param["max_edge_len"]+".");
        m_msg->write_config("nlist.max_edge_len",param["max_edge_len"]);
        m_max_edge_len =  lexical_cast<double>(param["max_edge_len"]);
      }
    }
    if (param.find("build_faces") != param.end())
    {
      if (m_system->get_periodic())
      {
        m_msg->msg(Messenger::ERROR,"Building faces is not supported for periodic boundary conditions.");
        throw runtime_error("Faces not supported in periodic sytems.");
      }
      m_build_contacts = true;
      m_build_faces = true;
      m_msg->msg(Messenger::INFO,"Neighbour list will also build faces.");
      m_msg->write_config("nlist.faces","true"); 
    }
    if (param.find("contact_distance") == param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list. No contact distance set. Assuming default sum of particle radii.");
      m_msg->write_config("nlist.contact_distance","0.0"); 
    }
    else
    {
      
      m_contact_dist = lexical_cast<double>(param["contact_distance"]);
      if (m_contact_dist > m_cut)
      {
        m_msg->msg(Messenger::WARNING,"Neighbour list. Contact distance "+param["contact_distance"]+" is larger than the neighbour list cuttoff distance. Using neigbour list cuttoff instead.");
        m_msg->write_config("nlist.contact_distance",lexical_cast<string>(m_cut)); 
        m_contact_dist = m_cut;
      }
      else
      {
        m_msg->msg(Messenger::INFO,"Neighbour list.  Setting contact distance to "+param["contact_distance"]+".");
        m_msg->write_config("nlist.contact_distance",param["contact_distance"]); 
      }
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
    if (param.find("remove_detached") != param.end())
    {
      m_msg->msg(Messenger::WARNING,"Neighbour list. Particles detached from the mesh (tissue) will be erased.");
      m_msg->write_config("nlist.remove_detached","true");
      m_remove_detached = true;
    }
    if (param.find("max_iter") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Setting maximum number of iterations for boundary builds in tissue simulations to "+param["max_iter"]+".");
      m_msg->write_config("nlist.max_iter",param["max_iter"]);
      m_system->set_max_mesh_iterations(lexical_cast<int>(param["max_iter"]));
    }
    if (param.find("boundary_type") != param.end())
    {
      m_msg->msg(Messenger::INFO,"Neighbour list. Setting typo of boundary particles in tissue simulations to "+param["boundary_type"]+".");
      m_msg->write_config("nlist.boundary_type",param["boundary_type"]);
      m_system->set_boundary_type(lexical_cast<int>(param["boundary_type"]));
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
  bool has_faces() { return m_build_faces; }
  
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
    m_contact_dist *= scale;
    if (m_use_cell_list && m_system->get_box()->Lx > 2.0*(m_cut+m_pad) && m_system->get_box()->Ly > 2.0*(m_cut+m_pad) && m_system->get_box()->Lz > 2.0*(m_cut+m_pad))
    {
      m_cell_list = boost::shared_ptr<CellList>(new CellList(m_system,m_msg,m_cut+m_pad));
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
  bool m_build_contacts;           //!< If true, build list of contacts
  bool m_build_faces;              //!< If true, build list of faces based on contact network
  bool m_triangulation;            //!< If true, build Delaunay triangulation for faces
  double m_contact_dist;           //!< Distance over which to assume particles to be in contact 
  double m_max_perim;              //!< Maximum value of the perimeter beyond which face becomes a hole.
  double m_max_edge_len;           //!< Maximum value of the edge beyond which we drop it (for triangulations)
  bool m_circumcenter;             //!< If true, use cell circumcenters when computing duals. 
  bool m_disable_nlist;            //!< If true, neigbour list is not built (only used for cell simulations)
  bool m_remove_detached;          //!< If true, remove detached particles (vertices) before rebuilding neighbour list (for cell simulations)
  vector<vector<int> >  m_contact_list;    //!< Holds the contact list for each particle
    
  // Actual neighbour list builds
  void build_nsq(int);    //!< Build with N^2 algorithm
  void build_cell();      //!< Build using cells list
  
  // Does contact list build
  void build_contacts();
  
  // Checks if the contact is intersecting with other contacts
  bool contact_intersects(int, int);
  
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
   
  // Check if all particles are on the same side
  bool same_side_line(Particle&, Particle&, vector<int>&);
  
  //! Dump particles and connectivity into a MOL2 file for debugging
  void debug_dump(const string&);

};

typedef shared_ptr<NeighbourList> NeighbourListPtr;

//! Dot product between two vectors given by three particles
double dot(const Particle&, const Particle&, const Particle&);

//! Mirror image of a point with respect to an edge
void mirror(const Particle&, const Particle&, const Particle&, double&, double&, double&);


#endif

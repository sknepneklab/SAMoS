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
 * \file vertex.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Declaration of Vertex class.
 */ 

#ifndef __VERTEX_HPP__
#define __VERTEX_HPP__

#include "particle.hpp"
#include "vector3d.hpp"

#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

#include <boost/format.hpp>

using boost::format;
using std::string;
using std::vector;
using std::endl;
using std::sin;
using std::cos;
using std::runtime_error;

/*! Vertex class is a light-weight class that handles vertices in a mesh.
 *  A vertex can simply be a particle position or the position of the dual lattice
 *  in tissue models. 
 * 
 *  This may apper as doubling of the data strucures, but gives us fexibilty to 
 *  work with meshes that do not coincide with the lattice positions.
 *
 */
struct Vertex
{
  //! Construct a Vertex object
  //! \param id vertex id
  //! \param x x coordinate
  //! \param y x coordinate
  //! \param z x coordinate
  Vertex(int id, double x, double y, double z) : id(id), 
                                                 type(1),
                                                 r(x,y,z),
                                                 N(0,0,1),
                                                 z(0),
                                                 n_edges(0), n_faces(0), 
                                                 boundary(false),
                                                 ordered(false),
                                                 attached(true) 
                                                 {   }
  
  //! Constract from particle position 
  //! \param p particle
  Vertex(Particle& p) 
  {
   id = p.get_id();
   type = p.get_type();
   r = Vector3d(p.x,p.y,p.z);
   if (p.Nx == 0.0 && p.Ny == 0.0 && p.Nz == 0.0)
     throw runtime_error("Surface normal for each mesh vertex has to non-zero."); 
   N = Vector3d(p.Nx,p.Ny,p.Nz);
   z = 0;
   n_edges = 0;
   n_faces = 0;
   boundary = false;
   ordered = false;
   attached = true;
  }
  
  ~Vertex()
  {
    neigh.clear();           
    edges.clear();           
    faces.clear();           
    dual.clear();            
  }
  
  //! Add neighbour
  //! \param v neighbour index
  void add_neighbour(int v)
  {
    neigh.push_back(v);
    z++;
  }
  
  //! Remove neighbour
  //! \param v neighbour index
  void remove_neighbour(int v)
  {
    vector<int>::iterator it = find(neigh.begin(), neigh.end(), v);
    if (it != neigh.end())
    {
      neigh.erase(it);
      z--;
    }
  }
  
  //! Add edge
  //! \param e edge index
  void add_edge(int e)
  {
    edges.push_back(e);
    n_edges++;
  }
   
  //! Remove edge
  //! \param e edge index
  void remove_edge(int e)
  {
    vector<int>::iterator it = find(edges.begin(), edges.end(), e);
    if (it != edges.end())
    {
      edges.erase(it);
      n_edges--;
    }
  }
  
  //! Add face
  //! \param f face index
  void add_face(int f)
  {
    faces.push_back(f);
    n_faces++;
  } 
  
  //! Remove face
  //! \param f face index
  void remove_face(int f)
  {
    vector<int>::iterator it = find(faces.begin(), faces.end(), f);
    if (it != faces.end())
    {
      faces.erase(it);
      n_faces--;
    }
  }
  
  //! Compare if the id of the vertex is equal to an integer value
  bool operator==(int val)
  {
     return (id == val);
  }
  
  //! Get the corresponding angle deficit derivative
  Vector3d& get_angle_def(int v)
  {
    if (v == id) return angle_def[0];  // note that vertices in the boundary star are ordered in a way to allow this assumption on indices
    else if (v == neigh[0]) return angle_def[1];
    else return angle_def[2]; 
  }
   
  int id;                      //!< Vertex id
  int type;                    //!< Vertex type 
  Vector3d r;                  //!< Position in the embedding 3d flat space
  Vector3d N;                  //!< Normal to the surface 
  
  int z;                       //!< Coordination number (number of neighbours)
  int n_edges;                 //!< Number of neighours this vertex has
  int n_faces;                 //!< Number of faces this vertex belongs to
  
  double area;                 //!< Area of dual associated with the vertex
  double perim;                //!< Perimeter of dual associated with the vertex
  
  bool boundary;               //!< If true, vertex is a boundary vertex
  bool ordered;                //!< If true, vertex star is ordered
  bool attached;               //!< If true, vertex is attached to a mesh
    
  vector<int> neigh;           //!< Contains indices of all neighbours
  vector<int> edges;           //!< Contains indices of all edges that originate at this vertex
  vector<int> faces;           //!< Contains indices of faces this vertex belongs to
  vector<int> dual;            //!< Centres of all faces surrounding it. Boundary vertices are special.
  vector<Vector3d> angle_def;  //!< Contains derivatives with respect to three coordinates at the boundary 
    
};

ostream& operator<<(ostream&, const Vertex&);

#endif
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
 * \file mesh.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Declaration of Mesh class.
 */ 

#ifndef __MESH_HPP__
#define __MESH_HPP__

#include "vertex.hpp"
#include "edge.hpp"
#include "face.hpp"

#include <vector>
#include <map>
#include <string>
#include <utility>
#include <algorithm>
#include <exception>
#include <cassert>

#include <boost/format.hpp>

using boost::format;
using std::string;
using std::vector;
using std::map;
using std::pair;
using std::make_pair;
using std::endl;
using std::copy;
using std::sort;
using std::reverse;
using std::rotate;
using std::runtime_error;

typedef pair<int,int> VertexPair;

//!< Data structure that holds data for ploting poygons
typedef struct
{
  vector<Vector3d> points;
  vector<vector<int> > sides;
  vector<double> area;
  vector<double> perim;
  vector<double> circum_radius;
} PlotArea;

/*! Mesh class handles basic manipuations with mesesh
 *
 */
class Mesh
{
public:
  //! Construct a Mesh object
  Mesh() : m_size(0), 
           m_nedge(0), 
           m_nface(0), 
           m_is_triangulation(true), 
           m_max_face_perim(20.0),
           m_circumcenter(true)
  {   }
  
  //! Get mesh size
  int size() { return m_size; }
  
  //! Get number of edges 
  int nedges() { return m_nedge; }
  
  //! Get number of faces
  int nfaces() { return m_nface; }
  
  //! Get list of all vertices
  vector<Vertex>& get_vertices() { return m_vertices; }
  
  //! Get edge list 
  vector<Edge>& get_edges() { return m_edges; }
  
  //! Get list of faces
  vector<Face>& get_faces() { return m_faces; }
  
  //! Get edge-face data structure
  map<pair<int,int>, int>& get_edge_face() { return m_edge_face; }
  
  //! Get the information about boundary vertex pairs
  vector<pair<int,int> >& get_boundary() { return m_boundary; }
  
  //! Get maximum face perimeter
  double get_max_face_perim() {  return m_max_face_perim; }
  
  //! Set value of the maximum face perimeter
  //! \param val new value
  void set_max_face_perim(double val) { m_max_face_perim = val; }
  
  //! Resets the mesh
  void reset();
  
  //! Sets the circumcenter flag
  //! \param val value of the circumcenter flag
  void set_circumcenter(bool val) { m_circumcenter = val; }
  
  //! Add a vertex
  //! \param p particle
  //add_vertex(Particle& p)
  //{
  //  vertices.push_back(Vertex(p));
  //  size++;
  //}
  
  //! Add a vertex
  //! \param vid vertex id
  //! \param x x-coordinate
  //! \param y y-coordinate
  //! \param z z-coordinate
  void add_vertex(int vid, double x, double y, double z)
  {
    m_vertices.push_back(Vertex(vid,x,y,z));
    m_size++;
  }
  
  //! Add from particle 
  //! \param p particle
  void add_vertex(Particle& p)
  {
    m_vertices.push_back(Vertex(p));
    m_size++;
  }
  
  //! Add an edge
  void add_edge(int,int);
    
  //! Generate faces of the mesh
  void generate_faces();
  
  //! Generate dual mesh
  void generate_dual_mesh();
  
  //! Update dual mesh
  void update_dual_mesh();
  
  //! Updates vertex positions 
  //! \param p particle
  void update(Particle& p)
  {
    m_vertices[p.get_id()].r = Vector3d(p.x, p.y, p.z);
    m_vertices[p.get_id()].type = p.get_type();
  }
  
  //! Post-processes the mesh
  void postprocess(bool);
  
  //! Compute face centre
  void compute_centre(int);
  
  //! Compute face angles
  void compute_angles(int);
  
  //! Order vertex star
  void order_star(int);
  
  //! Order dual star
  void order_dual(int);
  
  //! Compute dual area
  double dual_area(int);
  
  //! Dual perimeter
  double dual_perimeter(int);
  
  //! Opposite vertex
  int opposite_vertex(int);
  
  //! Flip edge
  void edge_flip(int);
  
  //! Mesh equiangulation
  bool equiangulate();
  
  //! Face centre Jacobian
  void fc_jacobian(int);
  
  //! Update face properties
  void update_face_properties();
  
  //! Remove obtuse boundary faces
  bool remove_obtuse_boundary();

  //! Fix obtuse boundary 
  vector<Vector3d> fix_obtuse_boundary();
  
  //! Remove edge triangles
  bool remove_edge_triangles();
  
  //! Return true if the vertex is a boundary vertex
  //! \param v index of the vertex
  bool is_boundary_vertex(int v)
  {
    return m_vertices[v].boundary;
  }
  
  //! Compute angle area scaling factor for boundary vertices
  double angle_factor(int);
  
  //! Compute derivatives of the angle factor for boudnary vertices
  void angle_factor_deriv(int);
  
  //! Compute radius of a circumscribed circle
  double circum_radius(int);
  
  //! Compute data for ploting polyons
  PlotArea& plot_area(bool);
     
private:  
  
  int m_size;    //!< Mesh size
  int m_nedge;   //!< Number of edges
  int m_nface;   //!< Number of faces
  bool m_is_triangulation;    //!< If true, all faces are triangles (allows more assumptions)
  double m_max_face_perim;    //!< If face perimeter is greater than this value, reject face and treat it as a hole.
  bool m_circumcenter;        //!< If true, compute face circumcenters. Otherwise compute geometric centre. 
    
  vector<Vertex> m_vertices;           //!< Contains all vertices
  vector<Edge> m_edges;                //!< Contains all edge
  vector<Face> m_faces;                //!< Contains all faces
  map<pair<int,int>, int> m_edge_map;  //!< Relates vertex indices to edge ids
  map<pair<int,int>, int> m_edge_face; //!< Relates pairs of faces to edges
  vector<pair<int,int> > m_boundary;   //!< List of vertex pair that are on the boundary
  vector<int> m_boundary_edges;        //!< List of all edges that are at the boundary
  vector<int> m_obtuse_boundary;       //!< List of all boundary edges that have obtuse angle opposite to them  
  PlotArea m_plot_area;                //!< Used to preapre polygonal data for plotting
  
  //! Compute face circumcentre
  void compute_circumcentre(int);
  
  //! Compute face geometric centre
  void compute_geometric_centre(int);
  
  //! Remove pair of edges
  bool remove_edge_pair(int);
  
  //! Compute face area
  double face_area(int);
  
  //! Order boundary star
  void order_boundary_star(int);
  
  //! Remove edge
  void remove_edge(int);
  
  //! Remove face
  void remove_face(int);
  
  //! Remove edge face
  bool remove_edge_face(int);

  //! Returns coordinates of the mirror image of a vertex opposite to a boundary edge
  Vector3d mirror_vertex(int);
  
  //! Functor used to compare lengths of two edges
  struct CompareEdgeLens
  {
    CompareEdgeLens(const Mesh& mesh) : m_mesh(mesh) { }
    bool operator()(const int e1, const int e2)
    {
      const Edge& E1 = m_mesh.m_edges[e1];
      const Edge& E2 = m_mesh.m_edges[e2];
      
      const Vertex& V1 = m_mesh.m_vertices[E1.from];
      const Vertex& V2 = m_mesh.m_vertices[E1.to];
      
      const Vertex& V3 = m_mesh.m_vertices[E2.from];
      const Vertex& V4 = m_mesh.m_vertices[E2.to];
      
      double r12 = (V1.r - V2.r).len();
      double r34 = (V3.r - V4.r).len();
      
      return (r12 > r34);
    }
    const Mesh& m_mesh;
  };
  
  //! Functor used to compare radii of circumscribed circles
  struct CompareRadii
  {
    CompareRadii(const Mesh& mesh) : m_mesh(mesh) { }
    bool operator()(const int e1, const int e2)
    {
      const Edge& E1 = m_mesh.m_edges[e1];
      const Edge& E2 = m_mesh.m_edges[e2];
      
      const Face& f1 = m_mesh.m_faces[m_mesh.m_edges[E1.pair].face];
      const Face& f2 = m_mesh.m_faces[m_mesh.m_edges[E2.pair].face];
      
      double r1 = (m_mesh.m_vertices[f1.vertices[0]].r-f1.rc).len();
      double r2 = (m_mesh.m_vertices[f2.vertices[0]].r-f2.rc).len();
        
      return (r1 > r2);
    }
    const Mesh& m_mesh;
  };
 
};

#endif

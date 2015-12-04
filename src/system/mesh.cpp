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
 * \file mesh.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Implementation of Mesh class memebr functions
 */ 

#include "mesh.hpp"

#include <iostream>

using std::cerr;
using std::cout;
using std::endl;
using std::unique;
using std::sort;

typedef pair<int,double> vert_angle;  //!< Used to sort angles

/*! Comprison used to sort edges in the star */
static bool comp(vert_angle v1, vert_angle v2)
{
  return v1.second < v2.second;
}


// Class members defined below

/*! Clean up the entire mesh data structure. */
void Mesh::reset()
{
  m_size = 0;  m_nedge = 0;  m_nface = 0;
  m_vertices.clear();           
  m_edges.clear();              
  m_faces.clear();              
  m_edge_map.clear();  
}

/*! Add and edge to the list of edges. Edge is defined
 *  by the indices of two vertices that belong to it.
 *  We also populate the auxiliary data structure edge_map
 *  that allows us to quickly search over vertex pairs.
 *  \param ei index of the 1st vertex
 *  \param ej index of the 2nd vertex
*/
void Mesh::add_edge(int ei, int ej)
{
  m_edges.push_back(Edge(m_nedge,ei,ej));
  m_vertices[ei].add_edge(m_nedge);
  m_vertices[ei].add_neighbour(ej);
  m_edge_map[make_pair(ei,ej)] = m_nedge;
  if (m_edge_map.find(make_pair(ej,ei)) != m_edge_map.end())
    m_edges[m_edge_map[make_pair(ej,ei)]].pair = m_nedge;
  m_nedge++;
}

/*! Generates faces from the edge information
*/
void Mesh::generate_faces()
{
  for (int v = 0; v < m_size; v++)
    this->order_star(v);
  for (int i = 0; i < m_nedge; i++)
  {
    Face face = Face(m_nface);
    Edge& E = m_edges[i];
    if (!E.visited)
    {
      E.visited = true;
      int seed = E.from;
      face.add_vertex(seed);
      face.add_vertex(E.to);
      face.add_edge(E.id);
      int v = E.to;
      cout << "Face : " << m_nface << " --> " << seed << " " << v << " ";
      while (v != seed)
      {
        Vertex& V = m_vertices[v];
        for (unsigned int e = 0; e < V.edges.size(); e++)
        {
          Edge& Ej = m_edges[V.edges[e]];
          if (!Ej.visited)
          {
            Ej.visited = true;
            if (Ej.to != seed) face.add_vertex(Ej.to);
            face.add_edge(Ej.id);
            v = Ej.to;
            cout << v << " ";
            break;
          }
        }
      }
      cout << endl;
    }
    double perim = 0.0;
    for (unsigned int f = 0; f < face.vertices.size(); f++)
    {
      int f_n = ( f == face.vertices.size() - 1) ? 0 : f + 1;
      Vector3d r1 = m_vertices[face.vertices[f]].r;
      Vector3d r2 = m_vertices[face.vertices[f_n]].r;
      perim += (r1-r2).len();
    }
    if (perim < m_max_face_perim)
    {
      m_faces.push_back(face);
      for (unsigned int f = 0; f < face.vertices.size(); f++)
        m_vertices[face.vertices[f]].add_face(m_nface);
      for (unsigned int f = 0; f < face.edges.size(); f++)
        m_edges[face.edges[f]].face = m_nface;
      m_nface++;
    }
  }
}


/*! Once the mesh is read in, we need to set things like
 *  boundary flags.
*/ 
void Mesh::postprocess()
{
  for (int e = 0; e < m_nedge; e++)
  {
    Edge& E = m_edges[e];
    if (E.face == NO_FACE)
    {
      m_vertices[E.from].boundary = true;
      m_vertices[E.to].boundary = true;
    }
  }
}

/*! Computes geometric centre of a face. Coordinates are stored in 
 *  the Face object.
 *  \param id of the face 
*/
void Mesh::compute_centre(int f)
{
  Face& face = m_faces[f];
  double xc = 0.0, yc = 0.0, zc = 0.0;
  for (int i = 0; i < face.n_sides; i++)
  {
    xc += m_vertices[face.vertices[i]].r.x;
    yc += m_vertices[face.vertices[i]].r.y;
    zc += m_vertices[face.vertices[i]].r.z;
  }
  face.rc = Vector3d(xc/face.n_sides,yc/face.n_sides,zc/face.n_sides);
}

/*! Order faces, edges and neighbours in the vertex star. At this point it is not possible
 *  to determine if the order is clockwise or counterclockwise. This 
 *  will be corrected for once the normal to the vertex is known.
 *  \param v vertex index
*/
void Mesh::order_star(int v)
{
  Vertex& V = m_vertices[v];
  int size = V.edges.size();
  if (size > 1)
  {
    int e = V.edges[0];
    Edge& Ei = m_edges[e];
    Vector3d ri = m_vertices[Ei.to].r - V.r;
    vector<vert_angle> angles;
    angles.push_back(make_pair(e,0.0));
    for (int i = 1; i < size; i++)
    {
      e = V.edges[i];
      Edge& Ej = m_edges[e];
      Vector3d rj = m_vertices[Ej.to].r - V.r;
      angles.push_back(make_pair(e,angle(ri,rj,V.N)));
    }
    sort(angles.begin(),angles.end(),comp);
    for (int i = 0; i < size; i++)
      V.edges[i] = angles[i].first;
  }
  V.ordered = true;
}

/*! Compute dual area by using expression 
 *  \f$ A_i = \frac{1}{2}\sum_{\mu}\left(\vec r_\mu\times\vec r_{\mu+1}\right)\cdot\vec n_i \f$
 *  were \f$ \vec r_\mu \f$ is the coodinate of the cetre of face \f$ \mu \f$ and \f$ \vec n_i \f$ 
 *  is the normal vector to the vertex.
 *  \note We assume that faces are ordered, otherwise the result will be 
 *  wrong.
 *  \param v verex index
*/
double Mesh::dual_area(int v)
{
  if (!m_vertices[v].ordered)
    throw runtime_error("Vertex star has to be ordered before dual area can be computed.");
  Vertex& V = m_vertices[v];
  if (V.boundary) return 0.0;
  
  V.area = 0.0;
  for (int i = 0; i < V.n_edges; i++)
  {
    int j = ( i == V.n_edges-1) ? 0 : i + 1;
    Vector3d& r1 = m_faces[m_edges[V.edges[i]].face].rc;
    Vector3d& r2 = m_faces[m_edges[V.edges[j]].face].rc;
    Vector3d  rr = cross(r1,r2); 
    V.area += dot(rr,V.N.unit());
  }
   
  V.area *= 0.5;
  if (V.area < 0.0) V.area = -V.area;
  
  return V.area;
}

/*! Compute lenght of duals perimeter
 *  \f$ l_i = \sum_{\mu}\left|\vec r_\mu-\vec r_{\mu+1}\right| \f$
 *  were \f$ \vec r_\mu \f$ is the coodinate of the cetre of face \f$ \mu \f$.
 *  \note We assume that faces are ordered, otherwise the result will be 
 *  wrong.
 *  \param v verex index
*/
double Mesh::dual_perimeter(int v)
{
  if (!m_vertices[v].ordered)
    throw runtime_error("Vertecx star has to be ordered before dual premeter can be computed.");
  Vertex& V = m_vertices[v];
  if (V.boundary) return 0.0;
  
  V.perim  = 0.0;
  for (int i = 0; i < V.n_edges; i++)
  {
    int j = ( i == V.n_edges-1) ? 0 : i + 1;
    Vector3d& r1 = m_faces[m_edges[V.edges[i]].face].rc;
    Vector3d& r2 = m_faces[m_edges[V.edges[j]].face].rc;
    V.perim += (r1-r2).len();
  }
  return V.perim;
}



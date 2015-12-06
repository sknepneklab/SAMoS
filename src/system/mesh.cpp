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
  if (m_edge_map.find(make_pair(ei,ej)) == m_edge_map.end())
  {
    m_edges.push_back(Edge(m_nedge,ei,ej));
    m_vertices[ei].add_edge(m_nedge);
    m_vertices[ei].add_neighbour(ej);
    m_edge_map[make_pair(ei,ej)] = m_nedge;
    m_nedge++;
  }
}

/*! Generates faces from the edge information
*/
void Mesh::generate_faces()
{
  vector<vert_angle> angles;
  for (int i = 0; i < m_nedge; i++)
  {
    Face face = Face(m_nface);
    Edge& E = m_edges[i];
    if (!E.visited)
    {
      E.visited = true;
      int seed = E.from;
      int vn = E.to;
      face.add_vertex(seed);
      face.add_vertex(vn);
      face.add_edge(E.id);
      int vp = seed;
      while (vn != seed)
      {
        Vertex& Vp = m_vertices[vp];
        Vertex& Vn = m_vertices[vn];
        Vector3d ri = Vn.r - Vp.r;
        angles.clear();
        for (unsigned int e = 0; e < Vn.edges.size(); e++)
        {
          Edge& Ej = m_edges[Vn.edges[e]];
          if (!Ej.visited && Ej.to != vp)
          {
            Vector3d rj = m_vertices[Ej.to].r - Vn.r;
            angles.push_back(make_pair(e,M_PI - angle(ri,rj,Vn.N)));
          }
        }
        sort(angles.begin(),angles.end(),comp);
        Edge& Ej = m_edges[Vn.edges[angles[0].first]];
        Ej.visited = true;
        if (Ej.to != seed) face.add_vertex(Ej.to);
        face.add_edge(Ej.id);
        vp = vn;
        vn = Ej.to;
      }      
    }
    double perim = 0.0;
    for (unsigned int f = 0; f < face.vertices.size(); f++)
    {
      int f_n = ( f == face.vertices.size() - 1) ? 0 : f + 1;
      Vector3d r1 = m_vertices[face.vertices[f]].r;
      Vector3d r2 = m_vertices[face.vertices[f_n]].r;
      perim += (r1-r2).len();
    }
    if (face.vertices.size() > 0)
    {
      if (perim >= m_max_face_perim) face.is_hole = true;
      m_faces.push_back(face);
      for (unsigned int f = 0; f < face.vertices.size(); f++)
        m_vertices[face.vertices[f]].add_face(m_nface);
      for (unsigned int f = 0; f < face.edges.size(); f++)
        m_edges[face.edges[f]].face = m_nface;
      m_nface++;
    }
  }
}

/*! Genererate position of the dual vertices */
void Mesh::generate_dual_mesh()
{
  m_dual.clear();
  m_ndual = 0;
  for (int v = 0; v < m_size; v++)
    m_vertices[v].dual.clear();
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (!face.is_hole)
    {
      this->compute_centre(f);
      m_dual.push_back(face.rc);
      for (int e = 0; e < face.n_sides; e++)
      {
        Edge& E = m_edges[face.edges[e]];
        E.dual = m_ndual;
      }
      m_ndual++;
    }
    else
    {
      for (int e = 0; e < face.n_sides; e++)
      {
        Edge& E = m_edges[face.edges[e]];
        Vector3d rc = m_faces[m_edges[E.pair].face].rc;
        Vector3d r = m_vertices[E.to].r - m_vertices[E.from].r;
        Vector3d rn = mirror(m_vertices[E.from].r,r,rc);
        m_dual.push_back(rn);
        m_vertices[E.from].dual.push_back(m_ndual);
        m_vertices[E.to].dual.push_back(m_ndual);
        E.dual = m_ndual;
        m_ndual++;
      }
    }
  }
  m_ndual = m_dual.size();    
}

/*! Update position of the dual vertices */
void Mesh::update_dual_mesh()
{
  int i = 0;
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (!face.is_hole)
    {
      this->compute_centre(f);
      m_dual[i++] = face.rc;
    }
    else
    {
      for (int e = 0; e < face.n_sides; e++)
      {
        Edge& E = m_edges[face.edges[e]];
        Vector3d rc = m_faces[m_edges[E.pair].face].rc;
        Vector3d r = m_vertices[E.to].r - m_vertices[E.from].r;
        Vector3d rn = mirror(m_vertices[E.from].r,r,rc);
        m_dual[i++] = rn;
      }
    }
  }   
}


/*! Once the mesh is read in, we need to set things like
 *  boundary flags.
*/ 
void Mesh::postprocess()
{
  m_size = m_vertices.size();
  m_nedge = m_edges.size();
  m_nface = m_faces.size();
  //for (int i = 0; i < m_size; i++)
  //  if (m_vertices[i].edges.size() == 0) m_vertices[i].boundary = true;
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (face.is_hole)
    {
      for (unsigned int v = 0; v < face.vertices.size(); v++)
        m_vertices[face.vertices[v]].boundary = true;
      for (unsigned int e = 0; e < face.edges.size(); e++)
        m_edges[face.edges[e]].boundary = true;
    }
  }
  for (int e = 0; e < m_nedge; e++)
  {
    Edge& E = m_edges[e];
    Edge& Epair = m_edges[m_edge_map[make_pair(E.to,E.from)]];
    E.pair = Epair.id;
    Epair.pair = E.id;
  }
  for (int v = 0; v < m_size; v++)
    this->order_star(v);
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
  V.dual.clear();
  vector<int> edges;
  
  if (V.edges.size() > 0)
  {
    edges.push_back(V.edges[0]);
    while (edges.size() < V.edges.size())
    {
      int idx = edges.size() - 1;
      int face = m_edges[edges[idx]].face;
      for (unsigned int e = 0; e < V.edges.size(); e++)
      {
        Edge& Ej = m_edges[V.edges[e]];
        if (m_edges[Ej.pair].face == face)
        {
          edges.push_back(V.edges[e]);
          break;
        }
      }
    }
    copy(edges.begin(),edges.end(),V.edges.begin());
    for (int e = 0; e < V.n_edges; e++)
    {
      Edge& E = m_edges[V.edges[e]];
      if (!m_faces[E.face].is_hole)
      {
        V.dual.push_back(E.dual);
      }
      else
      {
        int prev = (e == 0) ? V.n_edges - 1 : e-1;
        int next = (e == V.n_edges - 1) ? 0 : e+1;
        Edge& Ep = m_edges[m_edges[V.edges[prev]].pair];
        Edge& En = m_edges[m_edges[V.edges[next]].pair];
        if (m_faces[Ep.face].is_hole)
        { 
          V.dual.push_back(Ep.dual);
          V.dual.push_back(E.dual);
        }
        else 
        {
          V.dual.push_back(E.dual);
          V.dual.push_back(En.dual);
        }
      }
    }
    //cout << V << endl;
  }
  else 
    V.attached = false;
    
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
  if (!V.attached) return 0.0;
  //if (V.boundary) return 0.0;
  
  V.area = 0.0;
  for (unsigned int i = 0; i < V.dual.size(); i++)
  {
    int j = ( i == V.dual.size()-1) ? 0 : i + 1;
    Vector3d& r1 = m_dual[V.dual[i]];
    Vector3d& r2 = m_dual[V.dual[j]];
    //Vector3d& r1 = m_faces[m_edges[V.edges[i]].face].rc;
    //Vector3d& r2 = m_faces[m_edges[V.edges[j]].face].rc;
    
    Vector3d  rr = cross(r1,r2); 
    V.area += dot(rr,V.N.unit());
  }
   
  V.area *= 0.5;
  if (V.area < 0.0) 
  {
    V.area = -V.area;
    reverse(V.dual.begin(),V.dual.end());
    reverse(V.edges.begin(),V.edges.end());
    reverse(V.neigh.begin(),V.neigh.end());
  }
  
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
  Vertex& V = m_vertices[v];
  if (!V.attached) return 0.0;
  if (!V.ordered)
    throw runtime_error("Vertex star has to be ordered before dual premeter can be computed.");
  //if (V.boundary) return 0.0;
  
  V.perim  = 0.0;
  for (int i = 0; i < V.dual.size(); i++)
  {
    int j = ( i == V.dual.size()-1) ? 0 : i + 1;
    //Vector3d& r1 = m_faces[m_edges[V.edges[i]].face].rc;
    //Vector3d& r2 = m_faces[m_edges[V.edges[j]].face].rc;
    Vector3d& r1 = m_dual[V.dual[i]];
    Vector3d& r2 = m_dual[V.dual[j]];
    V.perim += (r1-r2).len();
  }
  return V.perim;
}



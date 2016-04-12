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
  m_edge_face.clear();
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
  m_is_triangulation = true;
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
      int prev_edge = E.id;
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
        m_edges[prev_edge].next = Ej.id;
        prev_edge = Ej.id;
        vp = vn;
        vn = Ej.to;
        if (vn == seed) m_edges[prev_edge].next = E.id;   // original edge is the "next" edge for last edge in the loop 
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
    if (face.vertices.size() > 3 && (!face.is_hole))
      m_is_triangulation = false;
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
      this->compute_angles(f);
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
        //Vector3d rc = m_faces[m_edges[E.pair].face].rc;
        Vector3d r = m_vertices[E.to].r - m_vertices[E.from].r;
        Vector3d rn = m_vertices[E.from].r + r.scaled(0.5);//mirror(m_vertices[E.from].r,r,rc);
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

/*! Update position of the dual vertices as well as the cell centre Jacobian */
void Mesh::update_dual_mesh()
{
  int i = 0;
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (!face.is_hole)
    {
      this->compute_angles(f);
      this->compute_centre(f);
      m_dual[i++] = face.rc;
    }
    else
    {
      for (int e = 0; e < face.n_sides; e++)
      {
        Edge& E = m_edges[face.edges[e]];
        //Vector3d rc = m_faces[m_edges[E.pair].face].rc;
        Vector3d r = m_vertices[E.to].r - m_vertices[E.from].r;
        Vector3d rn = m_vertices[E.from].r + r.scaled(0.5); //mirror(m_vertices[E.from].r,r,rc);
        m_dual[i++] = rn;
      }
    }
    this->fc_jacobian(f);
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
  m_boundary.clear();
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
      {
        Edge& E = m_edges[face.edges[e]];
        E.boundary = true;
        m_boundary.push_back(make_pair(E.from,E.to));
        m_boundary.push_back(make_pair(E.to,E.from));
      }
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

/*! Computes centre of a face. If the face is not triangle, compute geometric centre.
 *  If the face is a triangle, compute circumcenter or geometric center depending on the
 *  value of the m_circumcenter flag. Always compute geometric centre for faces at the boundary.
 *  the Face object.
 *  \param id of the face 
*/
void Mesh::compute_centre(int f)
{
  Face& face = m_faces[f];
  if (!m_circumcenter && face.n_sides > 3)
    this->compute_geometric_centre(f);
  else
    this->compute_circumcentre(f); 
}

/*! Computes cosines of interior angles at each vertex of the face. 
 *  It assumes that vertices are ordered and computes the vertex angle as the 
 *  angle between vectors along two edges meeting at the vertex.
 *  \param f id of the face 
*/
void Mesh::compute_angles(int f)
{
  Face& face = m_faces[f];
  face.angles.clear();
  int i_m, i_p;
  for (int i = 0; i < face.n_sides; i++)
  {
    i_m = (i == 0) ? face.n_sides - 1 : i - 1;
    i_p = (i == face.n_sides-1) ? 0 : i + 1;
    Vector3d& ri = m_vertices[face.vertices[i]].r;
    Vector3d& ri_m = m_vertices[face.vertices[i_m]].r;
    Vector3d& ri_p = m_vertices[face.vertices[i_p]].r;
    Vector3d dr1 = (ri_p - ri).unit();
    Vector3d dr2 = (ri_m - ri).unit();
    face.angles.push_back(dr1.dot(dr2));
    //face.angles.push_back(std::abs(angle(dr1,dr2,m_vertices[face.vertices[i]].N)));
  }
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
    // Here we handle edges along the hole.
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
    
  V.area = 0.0;
  if (!V.boundary)
  {
    for (unsigned int i = 0; i < V.dual.size(); i++)
    {
      int j = ( i == V.dual.size()-1) ? 0 : i + 1;
      Vector3d& r1 = m_dual[V.dual[i]];
      Vector3d& r2 = m_dual[V.dual[j]];
            
      Vector3d  rr = cross(r1,r2); 
      V.area += dot(rr,V.N.unit());
    }
  }
  else
  {
    for (int f = 0; f < V.n_faces; f++)
    {
      Face& face = m_faces[V.faces[f]];
      int fn = (f == V.n_faces-1) ? 0 : f + 1;
      Face& face_n = m_faces[V.faces[fn]];
      if (!(face.is_hole || face_n.is_hole))
      {
        Vector3d rr = cross(V.r-face.rc,V.r-face_n.rc);
        V.area += fabs(dot(rr,V.N.unit()));
      }
      if (face.boundary && !face.obtuse)
      {
        for (int i = 0; i < face.n_sides; i++)
          if (m_edges[m_edges[face.edges[i]].pair].boundary) 
          {
            Edge& E = m_edges[face.edges[i]];
            Vector3d& rj = m_vertices[E.from].r;
            Vector3d& rk = m_vertices[E.to].r;
            Vector3d r = 0.5*(rj+rk);
            Vector3d  rr = cross(V.r-r,V.r-face.rc);
            V.area += fabs(dot(rr,V.N.unit()));
            break;
          }
      }
    }
  }
   
  V.area *= 0.5;
  if (V.area < 0.0) 
  {
    V.area = -V.area;
    reverse(V.dual.begin(),V.dual.end());
    reverse(V.edges.begin(),V.edges.end());
    reverse(V.neigh.begin(),V.neigh.end());
  }
  
  // Now we make sure that the start around boundary vertex is
  // ordered such that fist and last edges are boundary edges
  if (V.boundary)
  {
    int pos = 0;
    for (int i = 0; i < V.n_edges; i++)
    {
      Edge& E = m_edges[V.edges[i]];
      if (m_edges[E.pair].boundary)
      {
        pos = i;
        break;
      }
    }
    rotate(V.edges.begin(), V.edges.begin()+pos, V.edges.end());
    rotate(V.dual.begin(), V.dual.begin()+pos, V.dual.end());
    rotate(V.neigh.begin(), V.neigh.begin()+pos, V.neigh.end());
  }
  
  return V.area;
}

/*! Compute length of dual's perimeter
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
  if (!V.boundary)
  {
    for (unsigned int i = 0; i < V.dual.size(); i++)
    {
      int j = ( i == V.dual.size()-1) ? 0 : i + 1;
      //Vector3d& r1 = m_faces[m_edges[V.edges[i]].face].rc;
      //Vector3d& r2 = m_faces[m_edges[V.edges[j]].face].rc;
      Vector3d& r1 = m_dual[V.dual[i]];
      Vector3d& r2 = m_dual[V.dual[j]];
      if (!(V.boundary && j == 0))
        V.perim += (r1-r2).len();
    }
  }
  else
  {
    for (int f = 0; f < V.n_faces; f++)
    {
      Face& face = m_faces[V.faces[f]];
      int fn = (f == V.n_faces-1) ? 0 : f + 1;
      Face& face_n = m_faces[V.faces[fn]];
      if (!(face.is_hole || face_n.is_hole))
      {
        Vector3d rr = face.rc-face_n.rc;
        V.perim += rr.len();
      }
      if (face.boundary && !face.obtuse)
      {
        for (int i = 0; i < face.n_sides; i++)
          if (m_edges[m_edges[face.edges[i]].pair].boundary) 
          {
            Edge& E = m_edges[face.edges[i]];
            Vector3d& rj = m_vertices[E.from].r;
            Vector3d& rk = m_vertices[E.to].r;
            Vector3d r = 0.5*(rj+rk);
            Vector3d  rr = r-face.rc;
            V.perim += rr.len();
            break;
          }
      }
    }
  }
  return V.perim;
}

/*! Return index of the vertex oposite to an edge. 
 *  This only makes sense if the face is a triangle.
 *  If it is not, we throw a runtime error message.
 *  \param e edge index
*/
int Mesh::opposite_vertex(int e)
{
  Edge& edge = m_edges[e];
  if (!edge.boundary)
  {
    Face& f = m_faces[edge.face];
    if (f.n_sides > 3)
      throw runtime_error("Vertex opposite to an edge is only defined for triangular faces.");
    for (int i = 0; i < f.n_sides; i++)
      if ((f.vertices[i] != edge.from) && (f.vertices[i] != edge.to))
        return f.vertices[i];
    cout << edge << endl;
    cout << f << endl;
    throw runtime_error("Vertex opposite to an edge: Mesh is not consistent. There is likely a bug in how edges and faces are updated.");
  }
  return -1;  // If edge is boundary, return -1.
}

/*! This is one of the simplest mesh moves that changes it topology. 
 *  Namely we flip an edge (pair of half-edges) shred by two trangles. Needless to say,
 *  this move is only defiened for triangulations.
 *  \param e id of the half-edge to flip (its pair is also flipped)
*/
void Mesh::edge_flip(int e)
{
  if (!m_is_triangulation)
    return;  // Edge flip is only defined for triangulations
  
  Edge& E = m_edges[e];
  Edge& Ep = m_edges[E.pair];
    
  if (E.boundary || Ep.boundary)
    return;   // We cannot flip a boundary edge.
  
  Face& F  = m_faces[E.face];
  Face& Fp = m_faces[Ep.face];
  
//   cout << "########## before ##########" << endl;
//   cout << "------ F -------" << endl;
//   cout << F << endl;
//   cout << "------ Fp -------" << endl;
//   cout << Fp << endl;
  
  // Get four edges surrounding the edge to be flipped
  Edge& E1 = m_edges[E.next];
  Edge& E2 = m_edges[E1.next];
  Edge& E3 = m_edges[Ep.next];
  Edge& E4 = m_edges[E3.next];
  
  assert(E2.next == E.id && E4.next == Ep.id);
  
  // Get four vertices that are involved in the flip
  Vertex& V1 = m_vertices[E.from];
  Vertex& V2 = m_vertices[Ep.from];
  Vertex& V3 = m_vertices[this->opposite_vertex(E.id)];
  Vertex& V4 = m_vertices[this->opposite_vertex(Ep.id)];
  
//   cout << "------------- V1 ----------------" << endl;
//   cout << V1 << endl;
  // Update vetices of the flipped half-edge pair
  E.from = V4.id;  Ep.from = V3.id;
  E.to = V3.id;    Ep.to = V4.id;
  
  // Update who follows whom  
  E.next = E2.id;
  E2.next = E3.id;
  E3.next = E.id;
  
  Ep.next = E4.id;
  E4.next = E1.id;
  E1.next = Ep.id;
  
  // Update face info for edges
  E3.face = F.id;
  E1.face = Fp.id;
  
  // Update dual info
  E3.dual = E2.dual;
  E1.dual = E4.dual;
  
  // Update face 1
  F.vertices[0] = E.from;   F.edges[0] = E.id;
  F.vertices[1] = E2.from;  F.edges[1] = E2.id;
  F.vertices[2] = E3.from;  F.edges[2] = E3.id;
  this->compute_angles(F.id);
  this->compute_centre(F.id);
  
  // Update face 2
  Fp.vertices[0] = Ep.from;  Fp.edges[0] = Ep.id;
  Fp.vertices[1] = E4.from;  Fp.edges[1] = E4.id;
  Fp.vertices[2] = E1.from;  Fp.edges[2] = E1.id;
  this->compute_angles(Fp.id);
  this->compute_centre(Fp.id);
  
  // Now we need to clean up vertices and their neighbours
  V1.remove_neighbour(V2.id);
  V1.remove_edge(E.id);
  V1.remove_face(Fp.id);
  
  V2.remove_neighbour(V1.id);
  V2.remove_edge(Ep.id);
  V2.remove_face(F.id);
  
  V3.add_neighbour(V4.id);
  V4.add_neighbour(V3.id);
  
  V4.add_edge(E.id);
  V3.add_edge(Ep.id);
  
  V3.add_face(Fp.id);
  V4.add_face(F.id);
  
  // Finally, we need to update edge_map
  map<pair<int,int>, int>::iterator it = m_edge_map.find(make_pair(V1.id,V2.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
     
  it = m_edge_map.find(make_pair(V2.id,V1.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
  
  m_edge_map[make_pair(V3.id,V4.id)]  = Ep.id;
  m_edge_map[make_pair(V4.id,V3.id)]  = E.id;
  
  // Make sure that the vertex stars are all properly ordered
  
//   cout << "########## after ##########" << endl;
//   cout << "------ F -------" << endl;
//   cout << F << endl;
//   cout << "------ Fp -------" << endl;
//   cout << Fp << endl;
//   
//   cout << "------------- V1 ----------------" << endl;
//   cout << V1 << endl;
  
  
  this->order_star(V1.id);
  this->order_star(V2.id);
  this->order_star(V3.id);
  this->order_star(V4.id);
  
  // Update area and permiters for the dual
  this->dual_area(V1.id);   this->dual_perimeter(V1.id);
  this->dual_area(V2.id);   this->dual_perimeter(V2.id);
  this->dual_area(V3.id);   this->dual_perimeter(V3.id);
  this->dual_area(V4.id);   this->dual_perimeter(V4.id);
  
//   cout << "########## after ##########" << endl;
//   cout << "------ F -------" << endl;
//   cout << F << endl;
//   cout << "------ Fp -------" << endl;
//   cout << Fp << endl;
    
}

/*! Implements the equiangulation of the mesh. This is a procedure where 
 *  all edges that have the sum of their opposing angles larger than pi 
 *  flipped. This procedure is guaranteed to converge and at the end one 
 *  recovers a Delaunday triangulation. 
*/
void Mesh::equiangulate()
{
  //cout << "Entered equiangulate" << endl;
  if (!m_is_triangulation)
    return;   // We cannot equiangulate a non-triangular mesh
  //cout << "Still in equiangulate" << endl;
  bool flips = true;
  while (flips)
  {
    flips = false;
    for (int e = 0; e < m_nedge; e++)
    {
      Edge& E = m_edges[e];
      Edge& Ep = m_edges[E.pair];
      if (!(E.boundary || Ep.boundary))
      {
        Vertex& V1 = m_vertices[this->opposite_vertex(E.id)];
        Vertex& V2 = m_vertices[this->opposite_vertex(Ep.id)];
        Face& F1 = m_faces[E.face];
        Face& F2 = m_faces[Ep.face];
        double angle_1 = F1.get_angle(V1.id);
        double angle_2 = F2.get_angle(V2.id);
        //if (angle_1 + angle_2 > M_PI)
        if (angle_1 + angle_2 < 0.0)
        {
          this->edge_flip(E.id);
          flips = true;
        }
      }
    }
  }
}

/*! For a triangular mesh compute derivatives (gradients) of the
 *  position of the face centre with respect to the position of
 *  each individual triangle vertices. We assume that the face (triangle)
 *  is oriented counterclockwise.
 *  \param f id of the face
*/
void Mesh::fc_jacobian(int f)
{
  Face& face = m_faces[f];
  if (face.n_sides > 3)
    return;
  face.drcdr.clear();
  face.drcdr.push_back(Matrix3d());
  face.drcdr.push_back(Matrix3d());
  face.drcdr.push_back(Matrix3d());
  
  if (face.boundary && face.obtuse)
  {
    // p = i
    face.drcdr[0].M[0][0] = 0.0;  face.drcdr[0].M[0][1] = 0.0;  face.drcdr[0].M[0][2] = 0.0;
    face.drcdr[0].M[1][0] = 0.0;  face.drcdr[0].M[1][1] = 0.0;  face.drcdr[0].M[1][2] = 0.0;
    face.drcdr[0].M[2][0] = 0.0;  face.drcdr[0].M[2][1] = 0.0;  face.drcdr[0].M[2][2] = 0.0;

    // p = j
    face.drcdr[1].M[0][0] = 0.5;  face.drcdr[1].M[0][1] = 0.0;  face.drcdr[1].M[0][2] = 0.0;
    face.drcdr[1].M[1][0] = 0.0;  face.drcdr[1].M[1][1] = 0.5;  face.drcdr[1].M[1][2] = 0.0;
    face.drcdr[1].M[2][0] = 0.0;  face.drcdr[1].M[2][1] = 0.0;  face.drcdr[1].M[2][2] = 0.5;
    
    // p = k
    face.drcdr[2].M[0][0] = 0.5;  face.drcdr[2].M[0][1] = 0.0;  face.drcdr[2].M[0][2] = 0.0;
    face.drcdr[2].M[1][0] = 0.0;  face.drcdr[2].M[1][1] = 0.5;  face.drcdr[2].M[1][2] = 0.0;
    face.drcdr[2].M[2][0] = 0.0;  face.drcdr[2].M[2][1] = 0.0;  face.drcdr[2].M[2][2] = 0.5;
    
    return;
  }

  Vector3d& ri = m_vertices[face.vertices[0]].r;
  Vector3d& rj = m_vertices[face.vertices[1]].r;
  Vector3d& rk = m_vertices[face.vertices[2]].r;

  Vector3d rjk = rj - rk;
  Vector3d rki = rk - ri;
  Vector3d rij = ri - rj;
  
  double rjk_2 = rjk.len2(), rki_2 = rki.len2(), rij_2 = rij.len2();
  double L_2 = rjk_2 + rki_2 + rij_2;
  double lambda_1 = rjk_2*(L_2 - 2*rjk_2);
  double lambda_2 = rki_2*(L_2 - 2*rki_2);
  double lambda_3 = rij_2*(L_2 - 2*rij_2);
  
  double Lambda = lambda_1 + lambda_2 + lambda_3;
    
  Vector3d dl1_dri =  2*rjk_2*(-rki + rij);
  Vector3d dl2_dri = -2*(rjk_2 + rij_2 - 2*rki_2)*rki + 2*rki_2*rij;
  Vector3d dl3_dri =  2*(rjk_2 + rki_2 - 2*rij_2)*rij - 2*rij_2*rki;
  
  Vector3d dl1_drj =  2*(rki_2 + rij_2 - 2*rjk_2)*rjk - 2*rjk_2*rij;
  Vector3d dl2_drj =  2*rki_2*(rjk - rij);
  Vector3d dl3_drj = -2*(rjk_2 + rki_2 - 2*rij_2)*rij + 2*rij_2*rjk;
  
  Vector3d dl1_drk = -2*(rki_2 + rij_2 - 2*rjk_2)*rjk + 2*rjk_2*rki;
  Vector3d dl2_drk =  2*(rjk_2 + rij_2 - 2*rki_2)*rki - 2*rki_2*rjk;
  Vector3d dl3_drk =  2*rij_2*(-rjk + rki);
  
  Vector3d dLam_dri = dl1_dri + dl2_dri + dl3_dri;
  Vector3d dLam_drj = dl1_drj + dl2_drj + dl3_drj;
  Vector3d dLam_drk = dl1_drk + dl2_drk + dl3_drk;
  
  double Lambda_2 = Lambda*Lambda;
  double inv_Lambda_2 = 1.0/Lambda_2;
  Vector3d dl1_div_Lam_dri = inv_Lambda_2*(Lambda*dl1_dri - lambda_1*dLam_dri);
  Vector3d dl2_div_Lam_dri = inv_Lambda_2*(Lambda*dl2_dri - lambda_2*dLam_dri);
  Vector3d dl3_div_Lam_dri = inv_Lambda_2*(Lambda*dl3_dri - lambda_3*dLam_dri);
  
  Vector3d dl1_div_Lam_drj = inv_Lambda_2*(Lambda*dl1_drj - lambda_1*dLam_drj);
  Vector3d dl2_div_Lam_drj = inv_Lambda_2*(Lambda*dl2_drj - lambda_2*dLam_drj);
  Vector3d dl3_div_Lam_drj = inv_Lambda_2*(Lambda*dl3_drj - lambda_3*dLam_drj);
  
  Vector3d dl1_div_Lam_drk = inv_Lambda_2*(Lambda*dl1_drk - lambda_1*dLam_drk);
  Vector3d dl2_div_Lam_drk = inv_Lambda_2*(Lambda*dl2_drk - lambda_2*dLam_drk);
  Vector3d dl3_div_Lam_drk = inv_Lambda_2*(Lambda*dl3_drk - lambda_3*dLam_drk);
  
  double l1_div_Lam = lambda_1/Lambda, l2_div_Lam = lambda_2/Lambda, l3_div_Lam = lambda_3/Lambda;
  
  // p = i
  face.drcdr[0].M[0][0] = dl1_div_Lam_dri.x*ri.x + l1_div_Lam + dl2_div_Lam_dri.x*rj.x + dl3_div_Lam_dri.x*rk.x;
  face.drcdr[0].M[0][1] = dl1_div_Lam_dri.y*ri.x              + dl2_div_Lam_dri.y*rj.x + dl3_div_Lam_dri.y*rk.x;
  face.drcdr[0].M[0][2] = dl1_div_Lam_dri.z*ri.x              + dl2_div_Lam_dri.z*rj.x + dl3_div_Lam_dri.z*rk.x;

  face.drcdr[0].M[1][0] = dl1_div_Lam_dri.x*ri.y  + dl2_div_Lam_dri.x*rj.y              + dl3_div_Lam_dri.x*rk.y;
  face.drcdr[0].M[1][1] = dl1_div_Lam_dri.y*ri.y  + dl2_div_Lam_dri.y*rj.y + l1_div_Lam + dl3_div_Lam_dri.y*rk.y;
  face.drcdr[0].M[1][2] = dl1_div_Lam_dri.z*ri.y  + dl2_div_Lam_dri.z*rj.y              + dl3_div_Lam_dri.z*rk.y;
  
  face.drcdr[0].M[2][0] = dl1_div_Lam_dri.x*ri.z + dl2_div_Lam_dri.x*rj.z + dl3_div_Lam_dri.x*rk.z;
  face.drcdr[0].M[2][1] = dl1_div_Lam_dri.y*ri.z + dl2_div_Lam_dri.y*rj.z + dl3_div_Lam_dri.y*rk.z;
  face.drcdr[0].M[2][2] = dl1_div_Lam_dri.z*ri.z + dl2_div_Lam_dri.z*rj.z + dl3_div_Lam_dri.z*rk.z + l1_div_Lam;
  
  // p = j
  face.drcdr[1].M[0][0] = dl1_div_Lam_drj.x*ri.x + l2_div_Lam + dl2_div_Lam_drj.x*rj.x + dl3_div_Lam_drj.x*rk.x;
  face.drcdr[1].M[0][1] = dl1_div_Lam_drj.y*ri.x              + dl2_div_Lam_drj.y*rj.x + dl3_div_Lam_drj.y*rk.x;
  face.drcdr[1].M[0][2] = dl1_div_Lam_drj.z*ri.x              + dl2_div_Lam_drj.z*rj.x + dl3_div_Lam_drj.z*rk.x;

  face.drcdr[1].M[1][0] = dl1_div_Lam_drj.x*ri.y + dl2_div_Lam_drj.x*rj.y              + dl3_div_Lam_drj.x*rk.y;
  face.drcdr[1].M[1][1] = dl1_div_Lam_drj.y*ri.y + dl2_div_Lam_drj.y*rj.y + l2_div_Lam + dl3_div_Lam_drj.y*rk.y;
  face.drcdr[1].M[1][2] = dl1_div_Lam_drj.z*ri.y + dl2_div_Lam_drj.z*rj.y              + dl3_div_Lam_drj.z*rk.y;
  
  face.drcdr[1].M[2][0] = dl1_div_Lam_drj.x*ri.z + dl2_div_Lam_drj.x*rj.z + dl3_div_Lam_drj.x*rk.z;
  face.drcdr[1].M[2][1] = dl1_div_Lam_drj.y*ri.z + dl2_div_Lam_drj.y*rj.z + dl3_div_Lam_drj.y*rk.z;
  face.drcdr[1].M[2][2] = dl1_div_Lam_drj.z*ri.z + dl2_div_Lam_drj.z*rj.z + dl3_div_Lam_drj.z*rk.z + l2_div_Lam;
  
  // p = k
  face.drcdr[2].M[0][0] = dl1_div_Lam_drk.x*ri.x + l3_div_Lam + dl2_div_Lam_drk.x*rj.x + dl3_div_Lam_drk.x*rk.x;
  face.drcdr[2].M[0][1] = dl1_div_Lam_drk.y*ri.x              + dl2_div_Lam_drk.y*rj.x + dl3_div_Lam_drk.y*rk.x;
  face.drcdr[2].M[0][2] = dl1_div_Lam_drk.z*ri.x              + dl2_div_Lam_drk.z*rj.x + dl3_div_Lam_drk.z*rk.x;

  face.drcdr[2].M[1][0] = dl1_div_Lam_drk.x*ri.y + dl2_div_Lam_drk.x*rj.y              + dl3_div_Lam_drk.x*rk.y;
  face.drcdr[2].M[1][1] = dl1_div_Lam_drk.y*ri.y + dl2_div_Lam_drk.y*rj.y + l3_div_Lam + dl3_div_Lam_drk.y*rk.y;
  face.drcdr[2].M[1][2] = dl1_div_Lam_drk.z*ri.y + dl2_div_Lam_drk.z*rj.y              + dl3_div_Lam_drk.z*rk.y;
  
  face.drcdr[2].M[2][0] = dl1_div_Lam_drk.x*ri.z + dl2_div_Lam_drk.x*rj.z + dl3_div_Lam_drk.x*rk.z;
  face.drcdr[2].M[2][1] = dl1_div_Lam_drk.y*ri.z + dl2_div_Lam_drk.y*rj.z + dl3_div_Lam_drk.y*rk.z;
  face.drcdr[2].M[2][2] = dl1_div_Lam_drk.z*ri.z + dl2_div_Lam_drk.z*rj.z + dl3_div_Lam_drk.z*rk.z + l3_div_Lam;
}

/*! This is an auxiliary memebr fucntion that loops over all faces
 *  and updates its properties such as the boundary and obtuse flags.
**/
void Mesh::update_face_properties()
{
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    face.boundary = false;
    face.obtuse = false;
    if (!face.is_hole)
    {
      for (int i = 0; i < face.n_sides; i++)
        if (m_edges[m_edges[face.edges[i]].pair].boundary) // Face is at the boundary if the pair of one of its half-edges is a boundary edge
        { 
          face.boundary = true;
          break;
        }
      for (int i = 0; i < face.n_sides; i++)
        if (!m_vertices[face.vertices[i]].boundary && face.get_angle(face.vertices[i]) < 0.0)
        {
          face.obtuse = true;
          break;
        }
      if (face.boundary && face.obtuse)
      {
        for (int e = 0; e < face.n_sides; e++) 
        {
          Edge& E = m_edges[face.edges[e]];
          if (m_vertices[E.from].boundary && m_vertices[E.to].boundary)
          {
            Vector3d& rj = m_vertices[E.from].r;
            Vector3d& rk = m_vertices[E.to].r;
            face.rc = 0.5*(rj+rk);
            m_dual[E.dual] = face.rc;
            this->fc_jacobian(f);
            break;
          }
        }
      }
    }
  }
}

/*! Loop over all boundary faces. If the face is obtuse,
 *  remove the edge boundary edge. This leaves the face 
 *  information invalid. All faces need to be rebuilt.
*/
bool Mesh::remove_obtuse_boundary()
{
  bool removed = false;
  
  double average_area = 0.0;
  for (int f = 0; f < m_nface; f++)
    average_area += this->face_area(f); 
  
  average_area /= m_nface;
  
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    //cout << "average area : " << average_area << endl;
    //cout << face << endl;
    //if (face.boundary && face.obtuse && face.area < 0.1*average_area)
    //if (face.boundary && face.obtuse)
    if (face.boundary && face.obtuse && face.area < 0.1*average_area)
    {
      //cout << "TO REMOVE " << endl << face << endl;
      for (int e = 0; e < face.n_sides; e++)
      {
        if (m_edges[m_edges[face.edges[e]].pair].boundary) 
        {
          this->remove_edge_pair(face.edges[e]);
          removed = true;
          break;
        }
      }
      if (removed)
        for (int v = 0; v < face.n_sides; v++)
          m_vertices[face.vertices[v]].boundary = true;   // upon removal all vertices are boundary
    }
   }
   
   if (removed)
   {
     m_nface = 0;
     m_faces.clear();
     for (int e = 0; e < m_nedge; e++)
     {
        m_edges[e].visited = false;
        m_edges[e].boundary = false;
     }
   }
   
   return removed;
}

/*! Find the angle deficit at a boundary vertex. This is 
 *  the angle between two boundary edges that that originate from the vertex.
 *  \param v vertex id
*/
double Mesh::angle_deficit(int v)
{
  Vertex& V = m_vertices[v];
  if (!V.boundary)
    return 1.0;
  
  double angle = 0.0;
  for (int f = 0; f < V.n_faces; f++)
  {
    Face& face = m_faces[V.faces[f]];
    if (!face.is_hole)
      angle += std::acos(face.get_angle(v));
  }
  return 1.0 - (2.0*M_PI-angle)/(2.0*M_PI);
}


// Private members

/*! Computes circumcenter of a face (assumes that face is a triangle).
 *  Coordinates are stored in the Face object.
 *  Formula for circumcenter is given as:
 *  Let \f$ \mathbf{A} \f$, \f$ \mathbf{B} \f$, and \f$ \mathbf{C} \f$ be 3-dimensional points, which form the vertices of a triangle.
 *  If we define,
 *  \f$ \mathbf{a} = \mathbf{A}-\mathbf{C}, \f$ and 
 *  \f$ \mathbf{b} = \mathbf{B}-\mathbf{C}. \f$, the circumcentre is given as
 *  \f$ p_0 = \frac{(\left\|\mathbf{a}\right\|^2\mathbf{b}-\left\|\mathbf{b}\right\|^2\mathbf{a})
                      \times (\mathbf{a} \times \mathbf{b})}
                  {2 \left\|\mathbf{a}\times\mathbf{b}\right\|^2} + \mathbf{C} \f$.
 *  \param id of the face 
*/
void Mesh::compute_circumcentre(int f)
{
  Face& face = m_faces[f];
  if (face.n_sides > 3)
    return;
  
  Vector3d& ri = m_vertices[face.vertices[0]].r;
  Vector3d& rj = m_vertices[face.vertices[1]].r;
  Vector3d& rk = m_vertices[face.vertices[2]].r;

  Vector3d rjk = rj - rk;
  Vector3d rki = rk - ri;
  Vector3d rij = ri - rj;
  
  double rjk_2 = rjk.len2(), rki_2 = rki.len2(), rij_2 = rij.len2();
  double L_2 = rjk_2 + rki_2 + rij_2;
  double lambda_1 = rjk_2*(L_2 - 2*rjk_2);
  double lambda_2 = rki_2*(L_2 - 2*rki_2);
  double lambda_3 = rij_2*(L_2 - 2*rij_2);
  
  double Lambda = lambda_1 + lambda_2 + lambda_3;

  face.rc = (lambda_1/Lambda)*ri + (lambda_2/Lambda)*rj + (lambda_3/Lambda)*rk; 
}

/*! Computes geometric centre of a face. Coordinates are stored in 
 *  the Face object.
 *  \param id of the face 
*/
void Mesh::compute_geometric_centre(int f)
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

/*! This is an auxiliary function used to remove a pair of edges.
 *  It updates the edge and vertex information, but leaves faces 
 *  inconsistent. It is to be used in conjunction with the function
 *  the weeds out obtuse boundary faces. After its application, all faces
 *  have to be rebuilt.
 *  \param e index of one of the edges
*/
void Mesh::remove_edge_pair(int e)
{
  Edge& E = m_edges[e];
  Edge& Ep = m_edges[E.pair];
  
  Vertex& V1 = m_vertices[E.from];
  Vertex& V2 = m_vertices[Ep.from];
  
  V1.remove_neighbour(V2.id);
  V2.remove_neighbour(V1.id);
  
  V1.remove_edge(E.id);
  V2.remove_edge(Ep.id);
  
  map<pair<int,int>, int>::iterator it = m_edge_map.find(make_pair(V1.id,V2.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
     
  it = m_edge_map.find(make_pair(V2.id,V1.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
     
  int e1 = (E.id < Ep.id) ? E.id : Ep.id;
  int e2 = (e1 == E.id) ? Ep.id : E.id;
  
  // Remove first edge
  vector<Edge>::iterator it_e = find(m_edges.begin(), m_edges.end(), e1);
  if (it_e != m_edges.end()) m_edges.erase(it_e);
  
  // Remove second edge
  it_e = find(m_edges.begin(), m_edges.end(), e2);
  if (it_e != m_edges.end()) m_edges.erase(it_e);
  
  // Relabel edges
  m_nedge = m_edges.size();
  for (int ee = 0; ee < m_nedge; ee++)
  {
    m_edges[ee].id = ee;
    if (m_edges[ee].pair > e1 && m_edges[ee].pair < e2) m_edges[ee].pair -= 1;
    else if (m_edges[ee].pair > e2) m_edges[ee].pair -= 2;
  }
  
  // Relabel vertex edge info
  for (int i = 0; i < m_size; i++)
  {
     Vertex& V = m_vertices[i];
     for (int ee = 0; ee < V.n_edges; ee++)
        if (V.edges[ee] > e1 && V.edges[ee] < e2) V.edges[ee] -= 1;
        else if (V.edges[ee] > e2) V.edges[ee] -= 2;
  }
  
  // Relabel edge_map info
  for (it = m_edge_map.begin(); it != m_edge_map.end(); it++)
    if ((*it).second > e1 && (*it).second < e2) (*it).second -= 1;
    else if ((*it).second > e2) (*it).second -= 2;
  
  // Relabel face edge info
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    for (int ee = 0; ee < face.n_sides; ee++)
      if (face.edges[ee] > e1 && face.edges[ee] < e2) face.edges[ee] -= 1;
      else if (face.edges[ee] > e2) face.edges[ee] -= 2;
  }
}

/*! Compute ares of a face.
 *  \param f face id
*/
double Mesh::face_area(int f)
{
  Face& face = m_faces[f];
  
  Vector3d& r0 = m_vertices[face.vertices[0]].r;
  
  face.area = 0.0;
  for (int i = 1; i < face.n_sides-1; i++)
  {
    Vector3d& r1 = m_vertices[face.vertices[i]].r;
    Vector3d& r2 = m_vertices[face.vertices[i+1]].r;
    face.area += cross(r1-r0,r2-r0).len(); 
  }
  
  face.area *= 0.5;
  
  return face.area;
}

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
 * \file mesh.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 18-Nov-2015
 * \brief Implementation of Mesh class memebr functions
 */ 

#include "mesh.hpp"

using std::unique;
using std::sort;
using std::rotate;

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
 *  \param vi index of the 1st vertex
 *  \param vj index of the 2nd vertex
*/
void Mesh::add_edge(int vi, int vj)
{
  m_edges.push_back(Edge(m_nedge,vi,vj));
  m_vertices[vi].add_edge(m_nedge);
  m_vertices[vi].add_neighbour(vj);
  m_edge_map[make_pair(vi,vj)] = m_nedge;
  m_nedge++;
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
    if (face.vertices.size() > 0)
    {
      if (face.vertices.size() > 3) face.is_hole = true;
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
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (!face.is_hole)
    {
      this->compute_angles(f);
      this->compute_centre(f);
    }
  }    
}

/*! Update position of the dual vertices as well as the cell centre Jacobian */
void Mesh::update_dual_mesh()
{
  for (int f = 0; f < m_nface; f++)
  {
    Face& face = m_faces[f];
    if (!face.is_hole)
    {
      this->compute_angles(f);
      this->compute_centre(f);     
    }
    this->fc_jacobian(f);
  }   
}


/*! Once the mesh is read in, we need to set things like
 *  boundary flags.
 *  \param flag if true order vertex star
*/ 
void Mesh::postprocess(bool flag)
{
  m_size = m_vertices.size();
  m_nedge = m_edges.size();
  m_nface = m_faces.size();
  m_boundary.clear();
  m_boundary_edges.clear();
  for (int i = 0; i < m_size; i++)
    m_vertices[i].boundary = false;
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
        m_boundary_edges.push_back(E.id);
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
  if (flag)
  {
    for (int v = 0; v < m_size; v++)
      this->order_star(v);
  }
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
  
  bool geom = false;
  if (face.n_sides > 3)
    geom = true;
  
  if (geom) this->compute_geometric_centre(f);
  else this->compute_circumcentre(f);  
  
   
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

/*! Order faces, edges and neighbours in the vertex star. 
 *  \param v vertex index
*/
void Mesh::order_star(int v)
{
  Vertex& V = m_vertices[v];
  V.neigh.clear();
  V.faces.clear();
  vector<int> edges;
  
  if (V.edges.size() == 0)
    V.attached = false;
  else
  {
    if (!V.boundary)      // For internal vertices just pick the first edge
      edges.push_back(V.edges[0]);
    else                  // For a boundary vertex, first edge in the star is always the one with its pair being a boundary edge
    {
      for (unsigned int e = 0; e < V.edges.size(); e++)
      {
        Edge& E = m_edges[V.edges[e]];
        if (m_edges[E.pair].boundary)
        {
          edges.push_back(E.id);
          break;
        }
      }
    }
    int i = 0;
    while(edges.size() < V.edges.size())
    {
      Edge& Ei = m_edges[edges[i++]];
      Vertex& Vi = m_vertices[Ei.to];
      Vector3d ri = Vi.r - V.r;
      double min_angle = 2.0*M_PI;
      int next_edge;
      for (unsigned int e = 0; e < V.edges.size(); e++)
      {
        Edge& Ej = m_edges[V.edges[e]];
        if (find(edges.begin(), edges.end(), Ej.id) == edges.end())  // Check if the edge has not already been included (note that most likely we don't even need this test, since angles are ordered)
        {
          Vertex& Vj = m_vertices[Ej.to];
          Vector3d rj = Vj.r - V.r;
          double ang = angle(ri,rj,V.N);
          ang = (ang >= 0.0) ? ang : 2.0*M_PI+ang;
          if (ang < min_angle)
          {
            min_angle = ang;
            next_edge = Ej.id;
          }
        }
      }
      edges.push_back(next_edge);
    }
    copy(edges.begin(),edges.end(),V.edges.begin());
    // Here we handle neighbours and faces 
    for (int e = 0; e < V.n_edges; e++)
    {
      Edge& E = m_edges[V.edges[e]];
      V.neigh.push_back(E.to);
      V.faces.push_back(E.face);
    }
    
    // Vertex star is not ordered
    V.ordered = true;
    // Make sure that the star of the boundary is in proper order
    if (V.boundary)
      this->order_boundary_star(V.id);
    this->order_dual(v);
  }
  
}

/*! Order dual vertices in the star. This is important as it may happen that 
 *  triangles within the star become obtuse which may lead to wrong ordering of 
 *  dual centres and areas and perimeters not being properly ordered.
 *  We need to check this after every time step.
 *  \param v vertex id
*/
void Mesh::order_dual(int v)
{
  Vertex& V = m_vertices[v];
  V.dual.clear();
  V.dual_neighbour_map.clear();
  if (V.n_faces == 0)
    V.attached = false;
  if (!V.attached) return;

  if (!V.ordered)
  {
    cout << V << endl;
    throw runtime_error("Vertex star has to be ordered before we can order its dual.");
  }
  
  Vector3d fc(0.0,0.0,0.0);
  for (unsigned int f = 0; f < V.faces.size(); f++)
  {
    Face& face = m_faces[V.faces[f]];
    if (face.is_hole)
      fc += V.r;
    else
      fc += face.rc;
  }
  fc.scale(1.0/static_cast<double>(V.faces.size()));
  
  // Now we order duals. This is important for proper computation of areas and force
  if (!V.boundary)
  {
    V.dual.push_back(V.faces[0]);
    V.dual_neighbour_map.push_back(V.neigh[0]);
  }
  else
  {
    Vector3d r0 =  m_vertices[V.neigh[0]].r - V.r;
    double min_angle = M_PI;
    int fface = 0;
    for (unsigned int f = 0; f < V.faces.size(); f++)
    {
       Face& face = m_faces[V.faces[f]];
       if (!face.is_hole)
       {
         Vector3d r1 = face.rc - V.r;
         double ang = angle(r0,r1,V.N);
         if (ang > -0.5*M_PI && ang < min_angle)
         {
           min_angle = ang;
           fface = f;
         }
       }
     }
     V.dual.push_back(V.faces[fface]);
     V.dual_neighbour_map.push_back(V.neigh[fface]);
  }
  int i = 0;
  while(V.dual.size() < V.faces.size())
  {
    Face& Fi = m_faces[V.dual[i]];
    Vector3d ri;
    if (V.boundary)
      ri = Fi.rc - V.r;
    else
      ri = Fi.rc - fc;
    double min_angle = 2.0*M_PI;
    int next_dual;
    bool can_add = false;
    int face_index;
    for (unsigned int f = 0; f < V.faces.size(); f++)
    {
      Face& Fj = m_faces[V.faces[f]];
      if (!Fj.is_hole)
      {
        if (find(V.dual.begin(), V.dual.end(), Fj.id) == V.dual.end())  
        {
          Vector3d rj;
          if (V.boundary)
            rj = Fj.rc - V.r;
          else
            rj = Fj.rc - fc;
          double ang = angle(ri,rj,V.N);
          ang = (ang >= 0.0) ? ang : 2.0*M_PI+ang;
          if (ang < min_angle)
          {
            min_angle = ang;
            next_dual = Fj.id;
            face_index = f;
            can_add = true;
          }
        }
      }
    }
    if (can_add)
    {
      V.dual.push_back(next_dual);
      V.dual_neighbour_map.push_back(V.neigh[face_index]);
      i++;
    }
    if (!can_add && (static_cast<int>(V.dual.size()) != V.n_faces-1))
    {
      cerr << endl;
      cerr << "Unable to order dual vertices of vertex " << V.id << ". There is likely a problem with the mesh itself. ";
      cerr << "Such problems typically arise if the parameters lead to situations such as a very thin (one-cell-wide) neck of cells or ";
      cerr << "a part of the system trying to detach from the bulk, which is currently not supported. " << endl;
      //this->debug_dump("test.off");
      throw runtime_error("Unable to order vertex dual.");
    }
    // for boundary vertices add the hole to the end
    if (V.boundary && (static_cast<int>(V.dual.size()) == V.n_faces-1))
      V.dual.push_back(V.faces[V.n_faces-1]);
  }
  //Not necessary, faces and neighbours line up closely enough
  //rotate(V.dual_neighbour_map.begin(), V.dual_neighbour_map.begin()+1, V.dual_neighbour_map.end());

  // Make sure that dual has the same number of elements as there are faces in the star 
  assert(V.dual.size() == V.faces.size());

  // And update area
  this->dual_area(V.id);
}

/*! Compute dual area by using expression 
 *  \f$ A_i = \frac{1}{2}\sum_{\mu}\left(\vec r_\mu\times\vec r_{\mu+1}\right)\cdot\vec n_i \f$
 *  were \f$ \vec r_\mu \f$ is the coodinate of the cetre of face \f$ \mu \f$ and \f$ \vec n_i \f$ 
 *  is the normal vector to the vertex.
 *  \note We assume that faces are ordered, otherwise the result will be 
 *  wrong.
 *  \param v vertex index
*/
double Mesh::dual_area(int v)
{

  Vertex& V = m_vertices[v];
  if (!V.attached) return 0.0;

  if (!V.ordered)
  {
    cout << V << endl;
    throw runtime_error("Vertex star has to be ordered before dual area can be computed.");
  }

  V.area = 0.0;
  if (V.dual.size() < 3) return V.area;
  if (!V.boundary)
  {
    for (int f = 0; f < V.n_faces; f++)
    {
      int fn = (f == V.n_faces - 1) ? 0 : f + 1;
      Face& F = m_faces[V.dual[f]];
      Face& Fn = m_faces[V.dual[fn]];
      V.area += dot(cross(F.rc,Fn.rc),V.N);
    }
  }
  else
  {
    V.area = dot(cross(V.r,m_faces[V.dual[0]].rc),V.N);
    for (int f = 0; f < V.n_faces-2; f++)
    {
      Face& F = m_faces[V.dual[f]];
      Face& Fn = m_faces[V.dual[f+1]];
      V.area += dot(cross(F.rc,Fn.rc),V.N);
    }
    V.area += dot(cross(m_faces[V.dual[V.n_faces-2]].rc,V.r),V.N);
  }
  
  V.area *= 0.5;
  
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
  
  V.perim = 0.0;
  if (!V.boundary)
  {
    for (int f = 0; f < V.n_faces; f++)
    {
      int fn = (f == V.n_faces - 1) ? 0 : f + 1;
      Face& F = m_faces[V.dual[f]];
      Face& Fn = m_faces[V.dual[fn]];
      V.perim += (F.rc-Fn.rc).len();
    }
  }
  else
  {
    V.perim = (V.r-m_faces[V.dual[0]].rc).len();
    for (int f = 0; f < V.n_faces-2; f++)
    {
      Face& F = m_faces[V.dual[f]];
      Face& Fn = m_faces[V.dual[f+1]];
      V.perim += (F.rc-Fn.rc).len();
    }
    V.perim += (m_faces[V.dual[V.n_faces-2]].rc-V.r).len();
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
    throw runtime_error("Vertex opposite to an edge: Mesh is not consistent. Paramters are likley wrong and there are overlapping vertices causing face construction to fail.");
  }
  return -1;  // If edge is boundary, return -1.
}

/*! This is one of the simplest mesh moves that changes it topology. 
 *  Namely we flip an edge (pair of half-edges) shred by two trangles. Needless to say,
 *  this move is only defined for triangulations.
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
  
  this->order_star(V1.id);
  this->order_star(V2.id);
  this->order_star(V3.id);
  this->order_star(V4.id);
  
  // Update area and permiters for the dual
  this->dual_area(V1.id);   this->dual_perimeter(V1.id);
  this->dual_area(V2.id);   this->dual_perimeter(V2.id);
  this->dual_area(V3.id);   this->dual_perimeter(V3.id);
  this->dual_area(V4.id);   this->dual_perimeter(V4.id);
  
  // Update face angles and centers
  this->compute_angles(F.id);
  this->compute_centre(F.id);
  this->compute_angles(Fp.id);
  this->compute_centre(Fp.id);

  if ((V1.n_edges <= 2) || (V2.n_edges <= 2) || (V3.n_edges <= 2) || (V4.n_edges <= 2))
    this->m_has_dangling = true;
    
}

/*! Implements the equiangulation of the mesh. This is a procedure where 
 *  all edges that have the sum of their opposing angles larger than pi 
 *  flipped. This procedure is guaranteed to converge and at the end one 
 *  recovers a Delaunday triangulation. 
*/
bool Mesh::equiangulate()
{
  //cout << "Entered equiangulate" << endl;
  if (!m_is_triangulation)
    return true;   // We cannot equiangulate a non-triangular mesh
  //cout << "Still in equiangulate" << endl;
  this->m_has_dangling = false;
  bool flips = true;
  bool no_flips = true;
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
          no_flips = false;
        }
      }
    }
  }
  return no_flips;
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

/*! This is an auxiliary member function that loops over all faces
 *  and updates its properties such as the boundary and obtuse flags.
**/
void Mesh::update_face_properties()
{
  m_obtuse_boundary.clear();
  for (unsigned int e = 0; e < m_boundary_edges.size(); e++)
  {
    Edge& E = m_edges[m_boundary_edges[e]];
    Edge& Ep = m_edges[E.pair];
    if (m_vertices[E.from].n_faces < 3) m_vertices[E.from].attached = false;
    else m_vertices[E.from].attached = true;
    Face& face = m_faces[Ep.face];
    face.boundary = true;
    if (face.get_angle(this->opposite_vertex(Ep.id)) < 0)
    {
      face.obtuse = true;
      if (!E.attempted_removal)
        m_obtuse_boundary.push_back(E.id);
    }
  }
  
}

/*! Loop over all boundary faces. If the face is obtuse,
 *  remove the edge boundary edge. This leaves the face 
 *  information invalid. All faces need to be rebuilt.
*/
bool Mesh::remove_obtuse_boundary()
{
  bool no_removals = true;
  for (int e = 0; e < m_nedge; e++)
    m_edges[e].attempted_removal = false;
  this->update_face_properties();
  while (m_obtuse_boundary.size() > 0)
  {
    bool removed = this->remove_edge_pair(*(m_obtuse_boundary.begin()));
    this->update_face_properties();
    if (removed)
      no_removals = false;
  }
  return no_removals;
}

/*! Loop over all boundary edges. If the angle opposite to the boundary 
 * edge is obtuse, mirror the vertex belonging to the obtuse angle across the 
 * boundary and add two edges and one triangle. In the next step, equiangulation 
 * move will pick and flip this edge. This prevents sudden jumps in the force. 
 * This function returns list of newly added vertices, which will be used to 
 * add actual particles to the system. 
*/
vector<Vector3d> Mesh::fix_obtuse_boundary()
{
  vector<Vector3d> new_verts;
  for (vector<int>::iterator it_b = m_obtuse_boundary.begin(); it_b != m_obtuse_boundary.end(); it_b++)
  {
    Edge& E = m_edges[*it_b];
    Edge& Ep  = m_edges[E.pair];
    m_faces[Ep.face].boundary = false;
    Vertex& Vf = m_vertices[E.from];
    Vertex& Vt = m_vertices[E.to];
    Vector3d v = this->mirror_vertex(this->opposite_vertex(Ep.id));
    new_verts.push_back(v);
    // add vertex
    this->add_vertex(m_size,v.x,v.y,v.z);
    Vertex& V = *(m_vertices.end() - 1);
    V.boundary = true;
    V.neigh.push_back(Vt.id);
    V.neigh.push_back(Vf.id);
    Vt.neigh.insert(Vt.neigh.begin(),V.id);
    Vf.neigh.push_back(V.id);
    // add two edges 
    this->add_edge(Vt.id,V.id);
    this->add_edge(V.id,Vt.id);
    Edge& E1 = *(m_edges.end()-2);
    Edge& E2 = *(m_edges.end()-1);
    Vt.edges.insert(Vt.edges.begin(),E1.id);
    V.edges.push_back(E2.id);
    E1.boundary = false;
    E2.boundary = true;
    E1.pair = E2.id; 
    E2.pair = E1.id; 
    // add two more edges
    this->add_edge(V.id,Vf.id);
    this->add_edge(Vf.id,V.id);
    Edge& E3 = *(m_edges.end()-2);
    Edge& E4 = *(m_edges.end()-1);
    Vf.edges.push_back(E4.id);
    V.edges.push_back(E3.id);
    E3.boundary = false;
    E4.boundary = true;
    E3.pair = E4.id; 
    E4.pair = E3.id;
    E.boundary = false;
    // update list of boundary edges
    m_boundary_edges.erase(find(m_boundary_edges.begin(),m_boundary_edges.end(),E.id));
    m_boundary_edges.push_back(E4.id);
    m_boundary_edges.push_back(E2.id);
    // handle new face 
    Face face = Face(m_nface++);
    face.boundary = true;
    face.add_vertex(Vt.id);
    face.add_vertex(V.id);
    face.add_vertex(Vf.id);
    face.add_edge(E.id);
    face.add_edge(E1.id);
    face.add_edge(E3.id); 
    E1.face = face.id; 
    E3.face = face.id; 
    E2.face = E.face; 
    E4.face = E.face;
    V.faces.push_back(face.id);
    V.faces.push_back(E.face);
    V.n_faces = 2;
    Vt.faces.insert(Vt.faces.begin(),face.id);
    vector<int>::iterator it_v = find(Vf.faces.begin(),Vf.faces.end(),E.face);
    Vf.faces.insert(it_v,face.id);
    V.dual.push_back(face.id);
    V.dual.push_back(E.face);
    Vt.dual.insert(Vt.dual.begin(),face.id);
    Vf.dual.push_back(face.id);
    // Insert two new edges and new vertex into the hole face making sure ordering is not messed up
    vector<int>::iterator it_e = find(m_faces[E.face].edges.begin(),m_faces[E.face].edges.end(),E.id);
    m_faces[E.face].edges.insert(it_e,E4.id);
    it_e = find(m_faces[E.face].edges.begin(),m_faces[E.face].edges.end(),E.id);
    m_faces[E.face].edges.insert(it_e,E2.id);
    it_e = find(m_faces[E.face].edges.begin(),m_faces[E.face].edges.end(),E.id);
    m_faces[E.face].edges.erase(it_e);
    it_v = find(m_faces[E.face].vertices.begin(),m_faces[E.face].vertices.end(),Vt.id);
    m_faces[E.face].vertices.insert(it_v,V.id);
    E.face = face.id;
    face.n_sides = 3;
    m_faces.push_back(face);
    this->compute_angles(face.id);
    this->compute_centre(face.id);  
    cout << "Adding vertex : " << V << endl;
  }
  this->update_face_properties();
  this->update_dual_mesh();
  return new_verts;
}

/*! Remove edge triangles.
 *  A triangle is considered an edge triangle if it has at least 
 *  to of its edges with having boundary pairs. That is, if on of 
 *  its vertices have only two neighbours.
*/
bool Mesh::remove_edge_triangles()
{
  bool no_removals = true;
  bool done = false;
  while (!done)
  {
    done = true;
    for (int f = 0; f < m_nface; f++)
      if (this->remove_edge_face(f))
      {
        done = false;
        no_removals = false;
        break;
      }
  }
  return no_removals;
}

/*! Compute radius of a circumscribed circle 
 *  \param f face id
*/
double Mesh::circum_radius(int f)
{
  Face& face = m_faces[f];
  
  if (face.n_sides > 3)
    face.radius = 0.0;
  else
  {
    face.radius = (m_vertices[face.vertices[0]].r-face.rc).len();
  }
  return face.radius;
}

/*! This member function is used to produce data 
 *  for plotting polygons into a VTK file.
 *  \param boundary if true, include boundary vertices as well
*/
PlotArea& Mesh::plot_area(bool boundary)
{
  m_plot_area.points.clear();
  m_plot_area.sides.clear();
  m_plot_area.area.clear();
  m_plot_area.perim.clear();
  m_plot_area.circum_radius.clear();
  m_plot_area.boundary_faces.clear();
  m_plot_area.type.clear();
  map<int,int> bnd_vert;
  map<int,int> face_idx;
  vector<int> added_faces;
  map<int,vector<int> > bnd_neigh;  // collect points that belong to the boundary of the dual mesh (for AJM simulations output)
  vector<int> bnd_faces;
  int idx = 0;
  for (int v = 0; v < m_size; v++)
  {
    Vertex& V = m_vertices[v];
    if (V.attached)
      if (V.boundary && boundary)
      {
        bnd_vert[V.id] = idx++;
        m_plot_area.points.push_back(V.r);
        m_plot_area.circum_radius.push_back(0.0);
        //m_plot_area.type.push_back(V.type);
      }
  }
  
  // Here we collect all faces (cells) that need to be ploted
  idx = m_plot_area.points.size();
  for (int v = 0; v < m_size; v++)
  {
    Vertex& V = m_vertices[v];
    if (V.attached)
    {
      for (int f = 0; f < V.n_faces; f++)
      {
        Face& face = m_faces[V.dual[f]];
        if (!face.is_hole)
        {
          if (find(added_faces.begin(),added_faces.end(),face.id) == added_faces.end())
          {
            this->compute_circumcentre(face.id);
            m_plot_area.points.push_back(face.rc);
            m_plot_area.circum_radius.push_back(this->circum_radius(face.id));
            added_faces.push_back(face.id);
            face_idx[face.id] = idx++;
          }
          if (V.boundary)
            if (find(bnd_faces.begin(), bnd_faces.end(), face_idx[face.id]) == bnd_faces.end()) bnd_faces.push_back(face_idx[face.id]);
        }
     }
    }
  }

  Vector3d rc(0.0,0.0,0.0);
  for (int i = 0; i < bnd_faces.size(); i++)
    rc += m_plot_area.points[bnd_faces[i]];

  rc.scale(1.0/bnd_faces.size());

  vector<int> sides;
  for (int v = 0; v < m_size; v++)
  {
    Vertex& V = m_vertices[v];
    if (V.attached)
    {
      if (!V.boundary)
      {
        sides.clear();
        for (int f = 0; f < V.n_faces; f++)
          sides.push_back(face_idx[V.dual[f]]);
        m_plot_area.sides.push_back(sides);
        m_plot_area.area.push_back(this->dual_area(V.id));
        m_plot_area.perim.push_back(this->dual_perimeter(V.id));
        m_plot_area.type.push_back(V.type);
      }
      else
      {
        for (int f = 0; f < V.n_faces-2; f++)
        {
          if (find(bnd_neigh[face_idx[V.dual[f]]].begin(), bnd_neigh[face_idx[V.dual[f]]].end(), face_idx[V.dual[f+1]]) == bnd_neigh[face_idx[V.dual[f]]].end()) bnd_neigh[face_idx[V.dual[f]]].push_back(face_idx[V.dual[f+1]]); 
          if (find(bnd_neigh[face_idx[V.dual[f+1]]].begin(), bnd_neigh[face_idx[V.dual[f+1]]].end(), face_idx[V.dual[f]]) == bnd_neigh[face_idx[V.dual[f+1]]].end()) bnd_neigh[face_idx[V.dual[f+1]]].push_back(face_idx[V.dual[f]]); 
        }
      }
      if (V.boundary && boundary)
      {
        sides.clear();
        sides.push_back(bnd_vert[V.id]);
        for (int f = 0; f < V.n_faces-1; f++)
          sides.push_back(face_idx[V.dual[f]]);
        m_plot_area.sides.push_back(sides);
        m_plot_area.area.push_back(this->dual_area(V.id));
        m_plot_area.perim.push_back(this->dual_perimeter(V.id));
        m_plot_area.type.push_back(V.type);
      }
    }
  }

  // Order boundary face centres in the clockwise fashion
  vector<vert_angle> angles;
  for (int i = 0; i < bnd_faces.size(); i++)
  {
    double dx = m_plot_area.points[bnd_faces[i]].x - rc.x;
    double dy = m_plot_area.points[bnd_faces[i]].y - rc.y;
    angles.push_back(make_pair(bnd_faces[i],atan2(dy,dx)));
  }
  
  sort(angles.begin(), angles.end(), comp);
  
  for (int i = 0; i < angles.size(); i++)
    m_plot_area.boundary_faces.push_back(angles[i].first);
  
  reverse(m_plot_area.boundary_faces.begin(), m_plot_area.boundary_faces.end());
  
  /*
  map<int,vector<int> >::iterator it = bnd_neigh.begin();
  int start = (*it).first;
  Vector3d r1 = m_plot_area.points[start];
  int next = ((*it).second)[0];
  Vector3d r2 = m_plot_area.points[next];
  Vector3d r01 = r1 - rc;
  Vector3d r02 = r2 - rc;
  Vector3d rcross = cross(r01,r02);
  if (rcross.z > 0)
    next = ((*it).second)[1];

  // This is used for output in AJM format
  m_plot_area.boundary_faces.push_back(start);
  m_plot_area.boundary_faces.push_back(next);

  while (m_plot_area.boundary_faces.size() < bnd_neigh.size())
  {
    next = m_plot_area.boundary_faces[m_plot_area.boundary_faces.size()-1];
    int nnext = bnd_neigh[next][0];
    if (find(m_plot_area.boundary_faces.begin(), m_plot_area.boundary_faces.end(), nnext) == m_plot_area.boundary_faces.end()) m_plot_area.boundary_faces.push_back(nnext);
    else m_plot_area.boundary_faces.push_back(bnd_neigh[next][1]);
  }
  */
  return m_plot_area;
}

// Private members

/*! Computes circumcentre of a face (assumes that face is a triangle).
 *  Coordinates are stored in the Face object.
 *  Formula for circumcentre is given as:
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
bool Mesh::remove_edge_pair(int e)
{
  Edge& E = m_edges[e];
  Edge& Ep = m_edges[E.pair];
  vector<int> affected_vertices;   // Vertices that are affected by the removal and need reordering
 
  E.attempted_removal = true;
  Ep.attempted_removal = true;
  // We can only remove boundary edge pairs
  if (!E.boundary)
    return false;
   
  // Get the face to be removed
  Face& face = m_faces[Ep.face];
  // And its pair
  Face& face_pair = m_faces[E.face];
  
  assert(!face.is_hole);
  assert(face_pair.is_hole);
  assert(face.n_sides == 3);
  
  // Check is the face is regular (at least one of its vertices are internal)
  bool non_regular = true;
  for (int v = 0; v < face.n_sides; v++)
    non_regular = non_regular && m_vertices[face.vertices[v]].boundary;
  
  if (non_regular)
    return false;
  
  Vertex& V1 = m_vertices[E.from];
  Vertex& V2 = m_vertices[Ep.from];
  
  V1.remove_neighbour(V2.id);
  V2.remove_neighbour(V1.id);
  
  V1.remove_edge(E.id);
  V2.remove_edge(Ep.id);
  
  V1.remove_face(face.id);
  V2.remove_face(face.id);
  
  map<pair<int,int>, int>::iterator it = m_edge_map.find(make_pair(V1.id,V2.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
     
  it = m_edge_map.find(make_pair(V2.id,V1.id));
  if (it != m_edge_map.end())
     m_edge_map.erase(it);
  
  // Remove edge from the list of boundary edges
  m_boundary_edges.erase(find(m_boundary_edges.begin(),m_boundary_edges.end(),E.id));
  
  // Update edges and vertices for the face to be removed
  for (int v = 0; v < face.n_sides; v++)
  {
    Vertex& VV = m_vertices[face.vertices[v]];
    if (VV.id != V1.id && VV.id != V2.id)
    {
      VV.remove_face(face.id);
      VV.add_face(face_pair.id);
      face_pair.add_vertex(VV.id);
      VV.boundary = true;
    }
    affected_vertices.push_back(VV.id);
  }
  
  for (int e = 0; e < face.n_sides; e++)
  {
    Edge& EE = m_edges[face.edges[e]];
    if (EE.id != Ep.id)
    {
      EE.face = face_pair.id;
      EE.boundary = true;
      face_pair.add_edge(EE.id);
      m_boundary_edges.push_back(EE.id);
    }
  }
  
  // Identify edge ids to be removed     
  int e1 = (E.id < Ep.id) ? E.id : Ep.id;
  int e2 = (e1 == E.id) ? Ep.id : E.id;
  
  // Identify face id to be removed
  int f = face.id;
  
  // Remove first edge
  vector<Edge>::iterator it_e = find(m_edges.begin(), m_edges.end(), e1);
  assert(it_e != m_edges.end());
  m_edges.erase(it_e);
  
  // Remove second edge
  it_e = find(m_edges.begin(), m_edges.end(), e2);
  assert(it_e != m_edges.end());
  m_edges.erase(it_e);
  
  // Update total number of edges
  assert(m_nedge == static_cast<int>(m_edges.size()) + 2);
  m_nedge = m_edges.size();
  
  // Remove face
  vector<Face>::iterator it_f = find(m_faces.begin(), m_faces.end(), f); 
  assert(it_f != m_faces.end());
  m_faces.erase(it_f);
  
  // Update total number of faces
  assert(m_nface == static_cast<int>(m_faces.size()) + 1);
  m_nface = m_faces.size();
  
  // Relabel edges
  for (int ee = 0; ee < m_nedge; ee++)
  {
    if (m_edges[ee].id > e1 && m_edges[ee].id < e2) m_edges[ee].id--;
    else if (m_edges[ee].id > e2) m_edges[ee].id -= 2;
    if (m_edges[ee].pair > e1 && m_edges[ee].pair < e2) m_edges[ee].pair -= 1;
    else if (m_edges[ee].pair > e2) m_edges[ee].pair -= 2;
    if (m_edges[ee].next > e1 && m_edges[ee].next < e2) m_edges[ee].next--;
    else if (m_edges[ee].next > e2) m_edges[ee].next -= 2;
    if (m_edges[ee].face > f) m_edges[ee].face--;
  }
  
  // Relabel boundary edges
  for (unsigned int ee = 0; ee < m_boundary_edges.size(); ee++)
  {
    if (m_boundary_edges[ee] > e1 && m_boundary_edges[ee] < e2) m_boundary_edges[ee] -= 1;
    else if (m_boundary_edges[ee] > e2) m_boundary_edges[ee] -= 2;
  }
  
  
  // Relabel vertex edge and face info
  for (int i = 0; i < m_size; i++)
  {
    Vertex& V = m_vertices[i];
    assert(V.n_edges == V.n_faces);
    for (int ee = 0; ee < V.n_edges; ee++)
    {
      if (V.edges[ee] > e1 && V.edges[ee] < e2) 
        V.edges[ee] -= 1;
      else if (V.edges[ee] > e2)
        V.edges[ee] -= 2;
    }
    for (int ff = 0; ff < V.n_faces; ff++)
    {
      if (V.faces[ff] > f) 
        V.faces[ff]--;
      if (V.dual[ff] > f) 
        V.dual[ff]--;
    }
  }
  
  // Relabel edge_map info
  for (it = m_edge_map.begin(); it != m_edge_map.end(); it++)
  {
    if ((*it).second > e1 && (*it).second < e2) (*it).second -= 1;
    else if ((*it).second > e2) (*it).second -= 2;
  }
  
  // Relabel face edge info
  for (int ff = 0; ff < m_nface; ff++)
  {
    Face& face = m_faces[ff];
    face.id = ff;
    for (int ee = 0; ee < face.n_sides; ee++)
      if (face.edges[ee] > e1 && face.edges[ee] < e2) face.edges[ee] -= 1;
      else if (face.edges[ee] > e2) face.edges[ee] -= 2;
  }
  
  
  // Order affected vertices
  for (unsigned int v = 0; v < affected_vertices.size(); v++)
    this->order_star(affected_vertices[v]);
    
  return true;
  
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

/*! This is an auxiliary private member function that 
 *  ensures that the boundary star is ordered such that
 *  the "hole" face appears last in the list of all faces 
 * belonging to the vertex.
 * \param v vertex index
*/
void Mesh::order_boundary_star(int v)
{
  Vertex& V = m_vertices[v];
  if (!V.boundary)
    return;
  
  unsigned int pos = 0;
  for (unsigned int f = 0; f < V.faces.size(); f++)
  {
    Face& face = m_faces[V.faces[f]];
    if (face.is_hole)
    {
      pos = (f == V.faces.size()-1) ? 0 : f + 1;
      break;
    }
  }
  rotate(V.edges.begin(), V.edges.begin()+pos, V.edges.end());
  //rotate(V.dual.begin(), V.dual.begin()+pos, V.dual.end());
  rotate(V.neigh.begin(), V.neigh.begin()+pos, V.neigh.end());
  rotate(V.faces.begin(), V.faces.begin()+pos, V.faces.end());
  
}

/*! Remove a single edge and update lables of all affected 
 *  vertices, edges and faces.
 *  \note This function removes only one edge (NOT its pair!).
 *  Therefore, if called alone it will leave the mesh in an 
 *  inconsistent state. Caller must assure that edges are always removed in 
 *  pairs.
 *  \param e index of the edge to be removed
*/
void Mesh::remove_edge(int e)
{
  // Remove edge from the list of all edges
  vector<Edge>::iterator it_e = find(m_edges.begin(), m_edges.end(), e);
  assert(it_e != m_edges.end());
  m_edges.erase(it_e);
  
  // Update total number of edges
  assert(m_nedge == static_cast<int>(m_edges.size()) + 1);
  m_nedge = m_edges.size();
  
  // Relabel edges
  for (int ee = 0; ee < m_nedge; ee++)
  {
    if (m_edges[ee].id > e) m_edges[ee].id--;
    if (m_edges[ee].pair > e) m_edges[ee].pair--;
    if (m_edges[ee].next > e) m_edges[ee].next--;
  }
  
  // Relabel boundary edges
  for (unsigned int ee = 0; ee < m_boundary_edges.size(); ee++)
    if (m_boundary_edges[ee] > e) m_boundary_edges[ee]--;
    
 // Relabel vertex edge and face info
  for (int i = 0; i < m_size; i++)
  {
    Vertex& V = m_vertices[i];
    for (int ee = 0; ee < V.n_edges; ee++)
      if (V.edges[ee] > e)  
        V.edges[ee]--;
  }
  
  // Relabel edge_map info
  for (map<pair<int,int>, int>::iterator it = m_edge_map.begin(); it != m_edge_map.end(); it++)
    if ((*it).second > e) (*it).second--;
    
  // Relabel face edge info
  for (int ff = 0; ff < m_nface; ff++)
  {
    Face& face = m_faces[ff];
    for (int ee = 0; ee < face.n_sides; ee++)
      if (face.edges[ee] > e) face.edges[ee]--;
  }   
}

/*! Remove a face from the mesh and update corresponding vertex, edge
 *  and face information.
 *  \param f id of the face to be removed
*/
void Mesh::remove_face(int f)
{
  // Remove face
  vector<Face>::iterator it_f = find(m_faces.begin(), m_faces.end(), f); 
  assert(it_f != m_faces.end());
  m_faces.erase(it_f);
  
  // Update total number of faces
  assert(m_nface == static_cast<int>(m_faces.size()) + 1);
  m_nface = m_faces.size();
 
  // Relabel face labels
  for (int ff = 0; ff < m_nface; ff++)
    m_faces[ff].id = ff;
 
  // Relabel edge face info
  for (int ee = 0; ee < m_nedge; ee++)
    if (m_edges[ee].face > f) m_edges[ee].face--;
  
  // Relabel vertex edge and face info
  for (int i = 0; i < m_size; i++)
  {
    Vertex& V = m_vertices[i];
    for (int ff = 0; ff < V.n_faces; ff++)
    {
      if (V.faces[ff] > f)
        V.faces[ff]--;
      if (V.dual[ff] > f)
        V.dual[ff]--;
    }
  }
} 

/*! Remove edge face.
 *  This function is used in conjunction with remove_edge_triangles function 
 *  to remove all triangles that live at the boundary and and have one 
 *  of its vertices having only two neighbours. This situation leads to 
 *  ill-defined angle deficits.
 *  \param f id of the face to remove
 *  \return true is the face was removed
*/
bool Mesh::remove_edge_face(int f)
{
  Face& face = m_faces[f];
  
  if (face.n_sides > 3)
    return false;
    
  bool can_remove = false;
  int vid;
  for (int v = 0; v < face.n_sides; v++)
  {
    Vertex& V = m_vertices[face.vertices[v]];
    if (V.n_edges < 3) 
    {
      can_remove = true;
      vid = face.vertices[v];
      break;
    }
  }
  
  if (!can_remove)
    return false;    // This triangle is OK and we don't need to touch it.
  
  vector<int> affected_vertices;   // Vertices that are affected by the removal and need reordering
   
  // We are now set to do actual removal
  Vertex& V = m_vertices[vid];
  assert(V.n_faces == 2);
  int face_pair = (V.faces[0] == f) ? V.faces[1] : V.faces[0];
  assert(m_faces[face_pair].is_hole);

  // Identify edge opposite to this vertex 
  int opposite_edge;
  for (int ee = 0; ee < face.n_sides; ee++)
    if (this->opposite_vertex(face.edges[ee]) == V.id)
    {
      opposite_edge = face.edges[ee];
      break;
    }
 
  // Opposite edge becomes boundary edge
  assert(!m_edges[opposite_edge].boundary);
  m_edges[opposite_edge].boundary = true;
  m_edges[opposite_edge].face = face_pair;
  m_boundary_edges.push_back(opposite_edge);
  
  // Vertex is no longer attached
  V.attached = false;
  V.boundary = false;
  // Update neighbour information
  for (unsigned int v = 0; v < V.neigh.size(); v++)
  {
    assert(m_vertices[V.neigh[v]].boundary);
    m_vertices[V.neigh[v]].remove_neighbour(V.id);
    m_vertices[V.neigh[v]].remove_face(f);
    map<pair<int,int>, int>::iterator it = m_edge_map.find(make_pair(V.id,m_vertices[V.neigh[v]].id));
    if (it != m_edge_map.end())
      m_edge_map.erase(it);
    it = m_edge_map.find(make_pair(m_vertices[V.neigh[v]].id,V.id));
    if (it != m_edge_map.end())
      m_edge_map.erase(it);
    affected_vertices.push_back(V.neigh[v]);
  }
  V.neigh.clear();
  V.faces.clear();
  V.dual.clear();
  
  // Remove edges
  while (V.edges.size() > 0)
  {
    Edge& E = m_edges[V.edges[0]];
    Edge& Ep = m_edges[E.pair];
    V.remove_edge(E.id);
    m_vertices[Ep.from].remove_edge(Ep.id);
    if (E.boundary)
      m_boundary_edges.erase(find(m_boundary_edges.begin(),m_boundary_edges.end(),E.id));
    if (Ep.boundary)
      m_boundary_edges.erase(find(m_boundary_edges.begin(),m_boundary_edges.end(),Ep.id));
    
    // Order edge ids for easy removal
    int e1 = (E.id < Ep.id) ? E.id : Ep.id;
    int e2 = (e1 == E.id) ? Ep.id : E.id;
    this->remove_edge(e1);
    this->remove_edge(e2-1);  // note that -1 in here is due to edge ids being shifted in the previous removal
  }
  // Vertex V is now detached. We set umber of its faces and edges to zero
  V.n_edges = 0;
  V.n_faces = 0;
  
  // Finally remove the face
  this->remove_face(f);
  
  
  // Order affected vertices
  for (unsigned int v = 0; v < affected_vertices.size(); v++)
    this->order_star(affected_vertices[v]);
    
  return true;
  
}

/*! Compute coordinates of a vertex that is a mirror image (with respect to an edge)
 *  of a vertex opposite to a boundary edge edge.
 *  
 *  \param edge index of the edge 
 */
 Vector3d Mesh::mirror_vertex(int edge)
 {
   Edge& E = m_edges[edge];
   Vertex& V0 = m_vertices[this->opposite_vertex(edge)];
   Vertex& V1 = m_vertices[E.from];
   Vertex& V2 = m_vertices[E.to];
   Vector3d r21 = V2.r - V1.r;
   Vector3d r01 = V0.r - V1.r; 
   Vector3d rp = (dot(r21,r01)/r21.len2())*r21;
   Vector3d rn = r01 - rp;
   return Vector3d(V1.r + (r01 - 2*rn));
 }


/*! Dump entire mesh into an OFF file for debugging purposes. 
 *  \param name file name for output 
 */
void Mesh::debug_dump(const string& name)
{
  ofstream out(name.c_str());
  out << "OFF" << endl;
  out << "# SAMoS debug OFF file." << endl;
  out << m_vertices.size() << " " << m_faces.size()-1 << " 0" << endl;
  for (unsigned int v = 0; v < m_vertices.size(); v++)
  {
    Vertex& V = m_vertices[v];
    out << V.r.x << " " << V.r.y << " " << V.r.z << endl;
  }
  for (unsigned int f = 0; f < m_faces.size(); f++)
  {
    Face& F = m_faces[f];
    if (F.n_sides < 7)  // some number larger than 3 to catch bugs but omit outter face
    {
      out << F.n_sides << " ";
      for (int v = 0; v < F.n_sides; v++)
        out << F.vertices[v] << " ";
      out << endl;
    }
  }
  out.close();
}

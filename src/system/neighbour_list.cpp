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
 * \file neighbour_list.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of the NeighbourList class
*/

#include "neighbour_list.hpp"

/* Auxiliary function which checks the side particle is on with respect to a line
*/
static int check_side(Particle& pi, Particle& pj, Particle& pk)
{
  double val = (pj.x-pi.x)*(pk.y-pi.y) - (pj.y-pi.y)*(pk.x-pi.x);
  return (val < 0.0) ? -1 : 1;
}

/* Auxiliary function which checks if the projection onto the segment is within the segment
*/
static int check_projection(Particle& pi, Particle& pj, Particle& pk)
{
  double dx = pj.x - pi.x, dy = pj.y - pi.y, dz = pj.z - pi.z;
  double p_dot_AB = (pk.x-pi.x)*dx+(pk.y-pi.y)*dy+(pk.z-pi.z)*dz;
  double AB_2 = dx*dx + dy*dy + dz*dz;
  double t = p_dot_AB/AB_2;
  return (t >= 0.0 && t <= 1.0);
}

/* Get neighbour opposite to an edge
*/
static int get_opposite(const Delaunay::Face_handle f, const int i, const int j)
{
  int I = f->vertex(0)->info(), J = f->vertex(1)->info(), K = f->vertex(2)->info();
  if (I == i)
  {
    if (J == j) return K;
    else return J;
  }
  else if (J == i)
  {
    if (I == j) return K;
    else return I;
  }
  else
  {
    if (J == j) return I;
    else return J;
  }
}

/*! Build neighbour list.
*/
void NeighbourList::build()
{
 m_list.clear();
 
 if (m_remove_detached)
   this->remove_detached();
  
 for (int i = 0; i < m_system->size(); i++)
 {
   m_list.push_back(vector<int>());
   if (m_build_contacts)
     m_contact_list.push_back(vector<int>());
 }
 
 if (!m_disable_nlist)
 {
  if (m_use_cell_list) 
    this->build_cell();
  else
  {
    m_old_state.clear();
    for (int i = 0; i < m_system->size(); i++)
      this->build_nsq(i);
  }
 }
 else
 {
   m_old_state.clear();
   for (int i = 0; i < m_system->size(); i++)
    {
      Particle& pi = m_system->get_particle(i);
      m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
    }
 }
  
 this->build_mesh();
}

/*! Build faces of the mesh from the particle locations */
void NeighbourList::build_mesh()
{
  if (m_build_contacts)
  {
    m_contact_list.clear();
    for (int i = 0; i < m_system->size(); i++)
      m_contact_list.push_back(vector<int>());
#ifdef HAS_CGAL
    if (m_triangulation)
    {
      this->build_triangulation();
    }
    else
#endif
      this->build_contacts();
  }
  if (m_build_faces)
    this->build_faces(true);
}

// Private methods below

/* Do actual building. */

//! Builds neighbour list using N^2 (all pairs algorithm).
//! \param i id of the particle for which to build the list
void NeighbourList::build_nsq(int i)
{
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  
  Particle& pi = m_system->get_particle(i);
  pi.coordination = 0;
  for (int j = 0; j < pi.get_id(); j++)
  {
    Particle& pj = m_system->get_particle(j);
    bool exclude = false;
    double dx = pi.x - pj.x;
    double dy = pi.y - pj.y;
    double dz = pi.z - pj.z;
    m_system->apply_periodic(dx,dy,dz);
    d2 = dx*dx + dy*dy + dz*dz;
    if (m_system->has_exclusions())
      if (m_system->in_exclusion(pi.get_id(), pj.get_id()))
        exclude = true;
    if (d2 < cut2 && (!exclude))
      m_list[j].push_back(i);
    if (d2 < cut2)
    {
      double r = pi.get_radius() + pj.get_radius();
      if (d2 < r*r)
      {
        pi.coordination++;
        pj.coordination++;
      }
    }
    
  }
  m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
  
}

//! Build neighbour list using cell list
void NeighbourList::build_cell()
{
  int N = m_system->size();
  BoxPtr box = m_system->get_box();
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  
  m_old_state.clear();
  
  m_cell_list->populate();
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    pi.coordination = 0;
  }
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    int cell_idx = m_cell_list->get_cell_idx(pi);
    vector<int>& neigh_cells = m_cell_list->get_cell(cell_idx).get_neighbours();  // per design includes this cell as well
    for (vector<int>::iterator it = neigh_cells.begin(); it != neigh_cells.end(); it++)
    {
      Cell& c = m_cell_list->get_cell(*it);
      vector<int>& p_idx_vec = c.get_particles(); 
      for (unsigned int j = 0; j < p_idx_vec.size(); j++)
      {
        Particle& pj = m_system->get_particle(p_idx_vec[j]);
        bool exclude = false;
        if (pj.get_id() > pi.get_id())
        {
          double dx = pi.x - pj.x;
          double dy = pi.y - pj.y;
          double dz = pi.z - pj.z;
          m_system->apply_periodic(dx,dy,dz);
          d2 = dx*dx + dy*dy + dz*dz;
          if (m_system->has_exclusions())
            if (m_system->in_exclusion(pi.get_id(), pj.get_id()))
              exclude = true;
          if (d2 < cut2 && (!exclude))
            m_list[i].push_back(pj.get_id());
          if (d2 < cut2)
          {
            double r = pi.get_radius() + pj.get_radius();
            if (d2 < r*r)
            {
             pi.coordination++;
             pj.coordination++;
            }
          } 
        }
      }
    }
    m_old_state.push_back(PartPos(pi.x,pi.y,pi.z));
  }
}

//! Check is neighbour list of the given particle needs update
//! \param p particle to check 
//! \return true if the list needs update
bool NeighbourList::need_update(Particle& p)
{
  int id = p.get_id();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  
  double dx = m_old_state[id].x - p.x;
  double dy = m_old_state[id].y - p.y;
  double dz = m_old_state[id].z - p.z;
  
  if (periodic)
  {
    if (dx > box->xhi) dx -= box->Lx;
    else if (dx < box->xlo) dx += box->Lx;
    if (dy > box->yhi) dy -= box->Ly;
    else if (dy < box->ylo) dy += box->Ly;
    if (dz > box->zhi) dz -= box->Lz;
    else if (dz < box->zlo) dz += box->Lz;
  }
  
  if (dx*dx + dy*dy + dz*dz < 0.25*m_pad*m_pad)
    return false;
  else
    return true;
}

//! Builds contact list based on particle distance
void NeighbourList::build_contacts()
{
  if (m_disable_nlist)
    throw runtime_error("Nlist build has to be enabled to build contacts.");
  int N = m_system->size();
  double dist = m_contact_dist;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    double ri = pi.get_radius();
    vector<int>& neigh = this->get_neighbours(i);
    for (unsigned int j = 0; j < neigh.size(); j++)
    {
      Particle& pj = m_system->get_particle(neigh[j]);
      double rj = pj.get_radius();
      if (m_contact_dist == 0.0)  dist = ri+rj;
      double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
      m_system->apply_periodic(dx,dy,dz);
      double r_sq = dx*dx + dy*dy + dz*dz;
      if (r_sq <= dist*dist)
      {
          m_contact_list[i].push_back(pj.get_id());
          m_contact_list[pj.get_id()].push_back(i);
      }
    }
  }
  
  this->remove_dangling();
}

//! Check if two edges intersect
//! \param i index of fist particle
//! \param j index of second particle
bool NeighbourList::contact_intersects(int i, int j)
{
  Particle& v1_i = m_system->get_particle(i);
  Vector3d p(v1_i.x, v1_i.y, v1_i.z);
  Particle& v1_j = m_system->get_particle(j);
  Vector3d s(v1_j.x-v1_i.x, v1_j.y-v1_i.y, v1_j.z-v1_i.z);
  
  
  //! Check all contacts of neighbours of particle i
  vector<int>& neigh_i = this->get_neighbours(i);
  for (unsigned int n = 0; n < neigh_i.size(); n++)
  {
    if (!(i == neigh_i[n] || j == neigh_i[n]))
    {
      cout << i << " , " << j << " --> " << neigh_i[n] << endl;
      Particle& v2_i = m_system->get_particle(neigh_i[n]);
      Vector3d q(v2_i.x, v2_i.y, v2_i.z);
      Vector3d q_m_p = q - p;
      for (unsigned int k = 0; k < m_contact_list[neigh_i[n]].size(); k++)
      {
        if (!(i == m_contact_list[neigh_i[n]][k] || j == m_contact_list[neigh_i[n]][k]))
        {
          Particle& v2_j = m_system->get_particle(m_contact_list[neigh_i[n]][k]);
          Vector3d r(v2_j.x-v2_i.x, v2_j.y-v2_i.y, v2_j.z-v2_i.z);
          double r_cross_s = cross(r,s).len();
          cout << r_cross_s << endl;
          if (r_cross_s != 0)
          {
            double t = cross(q_m_p,s).len()/r_cross_s;
            double u = cross(q_m_p,r).len()/r_cross_s;
            cout << t << " " << u << endl;
            if ((t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0))
              return true;
          }
        }
      }
    }
  }
  
  //! Check all contacts of neighbours of particle i
  vector<int>& neigh_j = this->get_neighbours(j);
  for (unsigned int n = 0; n < neigh_j.size(); n++)
  {
    if (!(i == neigh_j[n] || j == neigh_j[n]))
    {
      cout << i << " , " << j << " --> " << neigh_j[n] << endl;
      Particle& v2_i = m_system->get_particle(neigh_j[n]);
      Vector3d q(v2_i.x, v2_i.y, v2_i.z);
      Vector3d q_m_p = q - p;
      for (unsigned int k = 0; k < m_contact_list[neigh_j[n]].size(); k++)
      {
        if (!(i == m_contact_list[neigh_j[n]][k] || j == m_contact_list[neigh_j[n]][k]))
        {
          Particle& v2_j = m_system->get_particle(m_contact_list[neigh_j[n]][k]);
          Vector3d r(v2_j.x-v2_i.x, v2_j.y-v2_i.y, v2_j.z-v2_i.z);
          double r_cross_s = cross(r,s).len();
          cout << r_cross_s << endl;
          if (r_cross_s != 0)
          {
            double t = cross(q_m_p,s).len()/r_cross_s;
            double u = cross(q_m_p,r).len()/r_cross_s;
            cout << t << " " << u << endl;
            if ((t >= 0.0 && t <= 1.0) && (u >= 0.0 && u <= 1.0))
              return true;
          }
        }
      }
    }
  }

  return false;
}

/*! Build faces using contact network 
 *  Assumes that contacts have been built. 
 *  \param flag if true do the postprocessing of the mesh
**/
void NeighbourList::build_faces(bool flag)
{
  Mesh& mesh = m_system->get_mesh();
  mesh.reset();
  int N = m_system->size();
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    mesh.add_vertex(pi);
    for (unsigned int j = 0; j < m_contact_list[i].size(); j++)
      mesh.add_edge(i,m_contact_list[i][j]);
  }
  
  try
  {
  mesh.set_circumcenter(m_circumcenter);
  mesh.set_max_face_perim(m_max_perim);
  mesh.generate_faces();
  mesh.generate_dual_mesh();
  mesh.postprocess(flag);
  // Here we remove obtuse and edge trianges
  // Strictly speaking, we should do this iteratively, until 
  // there are no more triangles or edges to remove. This would be
  // simple to implement, but we postpone it unless it causes problems in
  // actual simulations. 
  //mesh.remove_obtuse_boundary();
  //mesh.remove_edge_triangles();
  m_system->update_mesh();
  }   catch (...)
      {
        mesh.debug_dump("mesh_dump.off");
        throw runtime_error("Still an error.");
      }
}

#ifdef HAS_CGAL
/*! Use CGAL to build 2D triangulation. This can be done only in 
 *  the plane, so this function will check if all z-coordinates are zero
 *  before proceeding. Otherwise, it will raise an exception.
*/
bool NeighbourList::build_triangulation()
{
  vector< pair<Point,unsigned> > points;
  //vector<int> on_convex_hull;
  //vector<int>& boundary = m_system->get_boundary();
  int N = m_system->size();
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.z != 0.0) 
    {
      m_msg->msg(Messenger::ERROR,"Delaunay triangulation is only supported in plane. All z components must be set to zero.");
      throw runtime_error("Unable to build Delaunay triangulation for non-planar systems.");
    }
    points.push_back( make_pair( Point(pi.x,pi.y), pi.get_id() ) );
  }

  
  Delaunay triangulation;  
  triangulation.insert(points.begin(),points.end());
  
  /*
  ofstream out("dela.off");
  out << "OFF" << endl;
  out << "# SAMoS debug OFF file." << endl;
  int del_size = 0;
  for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); fit++)
    del_size++;
  out << m_system->size() << " " << del_size << " 0" << endl;
  for (int i = 0; i < m_system->size(); i++)
  {
    Particle& p = m_system->get_particle(i);
    out << p.x << " " << p.y << " " << p.z << endl;
  }
  for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); fit++)
  {
    Delaunay::Face_handle face = fit;
    int i = face->vertex(0)->info();
    int j = face->vertex(1)->info();
    int k = face->vertex(2)->info();
    out << 3 << " " << i << " " << j << " " << k << endl;
  }
  out.close();
  int idx = 1; 
  for(Delaunay::Finite_edges_iterator eit = triangulation.finite_edges_begin(); eit != triangulation.finite_edges_end(); eit++)
  {
    Delaunay::Face_handle f1 = eit->first;
    Delaunay::Face_handle f2 = f1->neighbor(eit->second);
    int i = f1->vertex(f1->cw(eit->second))->info();
    int j = f1->vertex(f1->ccw(eit->second))->info();
    if (triangulation.is_infinite(f1) || triangulation.is_infinite(f2))
      cout << idx++ << " " << i << " " << j << " " << triangulation.is_infinite(f1) << " " << triangulation.is_infinite(f2) << endl;
  }
  */
  for(Delaunay::Finite_edges_iterator eit = triangulation.finite_edges_begin(); eit != triangulation.finite_edges_end(); eit++)
  {
    Delaunay::Face_handle f1 = eit->first;
    Delaunay::Face_handle f2 = f1->neighbor(eit->second);
    int i = f1->vertex(f1->cw(eit->second))->info();
    int j = f1->vertex(f1->ccw(eit->second))->info();
    int k, l;
    int to_flip;
    bool can_add = false;
    bool simple_add = true;
    // if f1 or f2 are boudnary handle this first
    if (triangulation.is_infinite(f1) || triangulation.is_infinite(f2))
    {
      if (triangulation.is_infinite(f2))
        k = f1->vertex(eit->second)->info();                //  k is the index of the index opposite to the edge (i,j)
      else
        k = get_opposite(f2, i, j);
      Particle& pi = m_system->get_particle(i);
      Particle& pj = m_system->get_particle(j);
      Particle& pk = m_system->get_particle(k);
      can_add = !(pi.boundary && pj.boundary && pk.boundary);        // We can add this edge only if at least of the vertices is not boundary
      if (can_add && dot(pk, pi, pj) < 0)    // In this case, we know that the Internal vertex in pk, so we compute angle at it
      {
        simple_add = false;       
        to_flip = k;
      }
    }
    else
    {
      k = f1->vertex(eit->second)->info();                //  k is the index of the index opposite to the edge (i,j) for face f1
      l = get_opposite(f2, i, j);                         //  l is the index of the index opposite to the edge (i,j) for face f2 
      Particle& pi = m_system->get_particle(i);
      Particle& pj = m_system->get_particle(j);
      Particle& pk = m_system->get_particle(k);
      Particle& pl = m_system->get_particle(l);
      bool all_f1_boundary = pi.boundary && pj.boundary && pk.boundary;
      bool all_f2_boundary = pi.boundary && pj.boundary && pl.boundary;
      can_add = !(all_f1_boundary && all_f2_boundary);
      if (can_add)
      {
        if ((all_f1_boundary && !all_f2_boundary))
          if (dot(pl, pi, pj) < 0) 
          {
            simple_add = false; 
            to_flip = l;
          }
        if ((!all_f1_boundary && all_f2_boundary))
          if (dot(pk, pi, pj) < 0) 
          {
            simple_add = false; 
            to_flip = k;
          }
      }
    }
    if (can_add)
    {
      if (simple_add)
      {
        if (find(m_contact_list[i].begin(),m_contact_list[i].end(),j) == m_contact_list[i].end()) m_contact_list[i].push_back(j);
        if (find(m_contact_list[j].begin(),m_contact_list[j].end(),i) == m_contact_list[j].end()) m_contact_list[j].push_back(i);
      }
      else
      {
        int i1 = i, i2 = j, i3 = to_flip, i4;
        Particle& p1 = m_system->get_particle(i1);
        Particle& p2 = m_system->get_particle(i2);
        Particle& p3 = m_system->get_particle(i3);

        double x, y, z;                  // contains coordinates of mirrored particles 
        mirror(p3, p1, p2, x, y, z);     // compute poistion of mirrored particle 
        Particle p(m_system->size(),p3.get_type(), p3.get_radius());    // generate new particle with the "last" id and inhereted type and radus from p3
        i4 = p.get_id();                                            
        // set parameters for the new particle
        p.x = x; p.y = y; p.z = z;
        p.Nx = p3.Nx;  p.Ny = p3.Ny;  p.Nz = p3.Nz;
        p.nx = p3.nx;  p.ny = p3.ny;  p.nz = p3.nz;
        p.coordination = 0;
        p.groups.push_back("all");
        // Note: Make sure that new particle is added to all necessary groups. 
        p.boundary = true;
        if (p1.boundary_neigh.size() == 2)  p1.boundary_neigh[(p1.boundary_neigh[0] == i2) ? 0 : 1] = i4;
        if (p2.boundary_neigh.size() == 2)  p2.boundary_neigh[(p2.boundary_neigh[0] == i1) ? 0 : 1] = i4;
        p.boundary_neigh.push_back(i1);
        p.boundary_neigh.push_back(i2);
        p3.boundary = false; 
        p3.boundary_neigh.clear();
        // add it to the system
        m_system->add_particle(p);
        // Extend the contact and neighbour list
        m_contact_list.push_back(vector<int>());
        m_list.push_back(vector<int>());

        // update neighbour list
        this->build_nsq(p.get_id());
        
        if (find(m_contact_list[i3].begin(),m_contact_list[i3].end(),i1) == m_contact_list[i3].end()) m_contact_list[i3].push_back(i1);
        if (find(m_contact_list[i1].begin(),m_contact_list[i1].end(),i3) == m_contact_list[i1].end()) m_contact_list[i1].push_back(i3);
        
        if (find(m_contact_list[i3].begin(),m_contact_list[i3].end(),i2) == m_contact_list[i3].end()) m_contact_list[i3].push_back(i2);
        if (find(m_contact_list[i2].begin(),m_contact_list[i2].end(),i3) == m_contact_list[i2].end()) m_contact_list[i2].push_back(i3);
        
        if (find(m_contact_list[i4].begin(),m_contact_list[i4].end(),i1) == m_contact_list[i4].end()) m_contact_list[i4].push_back(i1);
        if (find(m_contact_list[i1].begin(),m_contact_list[i1].end(),i4) == m_contact_list[i1].end()) m_contact_list[i1].push_back(i4);
        
        if (find(m_contact_list[i4].begin(),m_contact_list[i4].end(),i2) == m_contact_list[i4].end()) m_contact_list[i4].push_back(i2);
        if (find(m_contact_list[i2].begin(),m_contact_list[i2].end(),i4) == m_contact_list[i2].end()) m_contact_list[i2].push_back(i4);
        
        if (find(m_contact_list[i4].begin(),m_contact_list[i4].end(),i3) == m_contact_list[i4].end()) m_contact_list[i4].push_back(i3);
        if (find(m_contact_list[i3].begin(),m_contact_list[i3].end(),i4) == m_contact_list[i3].end()) m_contact_list[i3].push_back(i4);
      }
    }

  }
  //this->debug_dump("after_contact_build.mol2");

  this->remove_dangling();
  this->remove_detached();
  
  //this->debug_dump("after_cleanup.mol2");

  return true;
}
#endif


//! Remove all edges that have two or less contacts
void NeighbourList::remove_dangling()
{
  bool done = false;
  int N = m_system->size();
  while(!done)
  {
    done = true;
    for (int i = 0; i < N; i++)
    {
      if (m_contact_list[i].size() == 2) 
      {
        int i1 = m_contact_list[i][0];
        int i2 = m_contact_list[i][1];
        Particle& p1 = m_system->get_particle(i1);
        Particle& p2 = m_system->get_particle(i2);
        if (p1.boundary && p2.boundary)
        {
          if (find(m_contact_list[i1].begin(), m_contact_list[i1].end(),i2) == m_contact_list[i1].end()) m_contact_list[i1].push_back(i2);
          if (find(m_contact_list[i2].begin(), m_contact_list[i2].end(),i1) == m_contact_list[i2].end()) m_contact_list[i2].push_back(i1);
        }
      }
      if ((m_contact_list[i].size() > 0) && (m_contact_list[i].size() <= 2))
      {
        for (unsigned int j = 0; j < m_contact_list[i].size(); j++)
        {
          vector<int>& v = m_contact_list[m_contact_list[i][j]];
          vector<int>::iterator it = find(v.begin(),v.end(),i);
          if (it != v.end()) v.erase(it);
        }
        m_contact_list[i].clear();
        done = false;
      }
    }
  }
}

/*! Loop over all vertices and remove particles that belong to cells that are detached from the
 *  rest of the tissue. This is clearly only possible is triangulation has been set, i.e., for
 *  tissue simulations. 
 */
 void NeighbourList::remove_detached()
 {
   Mesh& mesh = m_system->get_mesh();
   if (mesh.size() == 0) return; 
   
   vector<int> to_remove;
   int offset = 0;  // We need to shift vertex ids to match them in the removal
   for (unsigned int i = 0; i < m_contact_list.size(); i++)
     if (m_contact_list[i].size() == 0)
     {
       to_remove.push_back(i-offset);
       offset++;       
     }
   /*
   for (int v = 0; v < mesh.size(); v++)
   {
     Vertex& V = mesh.get_vertices()[v];
     if (!V.attached)
     {
       to_remove.push_back(V.id-offset);
       offset++;
     }
   }
   */
   for (unsigned int i = 0; i < to_remove.size(); i++)
   {
     m_system->remove_particle(to_remove[i]);
     if (m_contact_list[to_remove[i]].size() != 0)
       throw runtime_error("Trying to remove connected particle.");
     else
     {
       m_contact_list.erase(m_contact_list.begin() + to_remove[i]);
       for (unsigned int j = 0; j < m_contact_list.size(); j++)
       {
         for (unsigned int k = 0; k < m_contact_list[j].size(); k++)
           if (m_contact_list[j][k] > to_remove[i]) m_contact_list[j][k]--;
       }
     }
   }
 } 

/*! Aufiliary function which checks if all negbours are on the same side a line connecting two 
 *  particles. This is used in constructing meshes to make sure that some edges are not 
 *  removed by accident
**/
bool NeighbourList::same_side_line(Particle& pi, Particle& pj, vector<int>& neigh)
{
  int side = 0;
  bool same = true;
  for (unsigned int n = 0; n < neigh.size(); n++)
  {
    Particle& pk = m_system->get_particle(neigh[n]);
    if (pi.get_id() != pk.get_id() && pj.get_id() != pk.get_id())
      if (check_projection(pi,pj,pk))
      {
        if (side == 0)
          side = check_side(pi,pj,pk);
        else
        {
          if (side != check_side(pi,pj,pk))
          {
            //cout << pi.get_id() << " " <<  pj.get_id()  << " " << pk.get_id() << " " <<  side << "  " << check_side(pi,pj,pk) << endl;
            same = false;
            break;
          }
        }
      }
  }
  return same;
}

/*! Auxiliary function for computing dot product between two vectors defined by three points.
 *  This is used to to determine if a particle needs to be mirrored when building the intial triangulation.
 *  \param p1 particle 1
 *  \param p2 particle 2
 *  \param p3 particle 3
 *  Dot product is computed between vector (p2-p1) and vector (p3-p1)
*/
double dot(const Particle& p1, const Particle& p2, const Particle& p3)
{
  double x21 = p2.x - p1.x, y21 = p2.y - p1.y, z21 = p2.z - p1.z;
  double x31 = p3.x - p1.x, y31 = p3.y - p1.y, z31 = p3.z - p1.z;
  return (x21*x31 + y21*y31 + z21*z31); 
}

/*! Auxiliary function for finding a mirror image of a point relative to the edge defined by two other points
 *  This is used to to determine the position of the mirrored particle that sits opposite to a boundary edge 
 *  and has its associated angle being obtuse. 
 *  \param p1 particle 1 (mirror this particle)
 *  \param p2 particle 2 (with repsect to a line connecting this)
 *  \param p3 particle 3 (and this particle)
 *  \param x x coordiante of the mirrored particle
 *  \param y y coordiante of the mirrored particle
 *  \param z z coordiante of the mirrored particle
*/
void mirror(const Particle& p1, const Particle& p2, const Particle& p3, double& x, double& y, double& z)
{
  double x32 = p3.x - p2.x, y32 = p3.y - p2.y, z32 = p3.z - p2.z;
  double len = sqrt(x32*x32 + y32*y32 + z32*z32);
  double Nx = x32/len, Ny = y32/len, Nz = z32/len;
  double n_xx = Nx*Nx, n_xy = Nx*Ny, n_xz = Nx*Nz;
  double n_yx = Ny*Nx, n_yy = Ny*Ny, n_yz = Ny*Nz;
  double n_zx = Nz*Nx, n_zy = Nz*Ny, n_zz = Nz*Nz;
  
  double P_m_Q_x = p2.x - p1.x, P_m_Q_y = p2.y - p1.y, P_m_Q_z = p2.z - p1.z;
  x = p1.x + 2.0*((1.0 - n_xx)*P_m_Q_x         - n_xy*P_m_Q_y -         n_xz*P_m_Q_z);
  y = p1.y + 2.0*(      - n_yx*P_m_Q_x + (1.0 - n_yy)*P_m_Q_y -         n_yz*P_m_Q_z); 
  z = p1.z + 2.0*(      - n_zx*P_m_Q_x         - n_zy*P_m_Q_y + (1.0 - n_zz)*P_m_Q_z); 
}

/*! Dump state of the system into a MOL2 file. Internal paricles will have type 1, boundary type 2
 *  bonds insed will be type 1, and bonds on the boundary will be type 2. 
 *  \param name name of the MOL2 file
 */
 void NeighbourList::debug_dump(const string& name)
 {
   ofstream out(name.c_str());
  int N = m_system->size();
  int Nbonds = 0;
  for (int i = 0; i < N; i++)
    Nbonds += m_contact_list[i].size();
  Nbonds /= 2;
  out << "@<TRIPOS>MOLECULE" << endl;
  out << "Generated by SAMoS code" << endl;
  out << N << "  " << Nbonds << endl;
  out << "NO_CHARGES" << endl;
  out << "@<TRIPOS>ATOM" << endl;
  for (int i = 0; i < N; i++)
  {
    Particle& p = m_system->get_particle(i);
    if (!p.boundary)
      out << format("%d\t%d\t%10.6f\t%10.6f\t%10.6f\t%d") % (p.get_id()+1) % 1 % p.x % p.y % p.z % p.get_type() << endl;
    else
      out << format("%d\t%d\t%10.6f\t%10.6f\t%10.6f\t%d") % (p.get_id()+1) % 2 % p.x % p.y % p.z % p.get_type() << endl;
  }
  out << "@<TRIPOS>BOND" << endl;
  int bidx = 1;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    for (unsigned int j = 0; j < m_contact_list[i].size(); j++)
    {
      Particle& pj = m_system->get_particle(m_contact_list[i][j]);
      if (pj.get_id() != m_contact_list[i][j])
        throw runtime_error("Index/neighbour mismatch.");
      if (pi.get_id() < pj.get_id())
      {
        if (find(pi.boundary_neigh.begin(),pi.boundary_neigh.end(),pj.get_id()) == pi.boundary_neigh.end())
          out << format("%d\t%d\t%d\t%d") % (bidx++) % (pi.get_id() + 1) % (pj.get_id() + 1) % 1 << endl;
        else
          out << format("%d\t%d\t%d\t%d") % (bidx++) % (pi.get_id() + 1) % (pj.get_id() + 1) % 2 << endl;
      }
    }
  }
   out.close();
 }


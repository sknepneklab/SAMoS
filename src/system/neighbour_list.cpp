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

/*! Build neighbour list using N^2 algorithm 
*/
void NeighbourList::build()
{
 //cout << "Building NL. " << endl;
 Mesh& mesh = m_system->get_mesh();
 m_list.clear();
 
 if (m_build_contacts)
 {
   mesh.reset();
   m_contact_list.clear();
 }
 
 if (m_build_faces)
   mesh.reset();
  
 for (int i = 0; i < m_system->size(); i++)
 {
   Particle& p = m_system->get_particle(i);
   p.coordination = 0;
   m_list.push_back(vector<int>());
   if (m_build_contacts)
     m_contact_list.push_back(vector<int>());
 }
  
 if (m_use_cell_list) 
    this->build_cell();
  else
    this->build_nsq();
  
  if (m_build_contacts)
  {
#ifdef HAS_CGAL
    if (m_triangulation)
      this->build_triangulation();
    else
#endif
      this->build_contacts();
  }
  if (m_build_faces)
    this->build_faces();
}

/*! Build faces of the mesh from the particle locations */
bool NeighbourList::build_faces()
{
  return this->build_mesh();
}

// Private methods below

/* Do actual building. */

//! Builds neighbour list using N^2 (all pairs algorithm).
void NeighbourList::build_nsq()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  Mesh& mesh = m_system->get_mesh();
  
  m_old_state.clear();
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_build_contacts)
      mesh.add_vertex(pi);
    for (int j = i + 1; j < N; j++)
    {
      Particle& pj = m_system->get_particle(j);
      bool exclude = false;
      double dx = pi.x - pj.x;
      double dy = pi.y - pj.y;
      double dz = pi.z - pj.z;
      if (periodic)
      {
        if (dx > box->xhi) dx -= box->Lx;
        else if (dx < box->xlo) dx += box->Lx;
        if (dy > box->yhi) dy -= box->Ly;
        else if (dy < box->ylo) dy += box->Ly;
        if (dz > box->zhi) dz -= box->Lz;
        else if (dz < box->zlo) dz += box->Lz;
      }
      d2 = dx*dx + dy*dy + dz*dz;
      if (m_system->has_exclusions())
        if (m_system->in_exclusion(pi.get_id(), pj.get_id()))
          exclude = true;
      if (d2 < cut2 && (!exclude))
        m_list[i].push_back(j);
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
  //m_msg->msg(Messenger::INFO, "Rebuilt neighbour list (N^2 algorithm).");
}

//! Build neighbour list using cell list
void NeighbourList::build_cell()
{
  int N = m_system->size();
  bool periodic = m_system->get_periodic();
  BoxPtr box = m_system->get_box();
  double cut = m_cut+m_pad;
  double cut2 = cut*cut;
  double d2;
  Mesh& mesh = m_system->get_mesh();
  
  m_old_state.clear();
  
  m_cell_list->populate();
  
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (m_build_contacts)
      mesh.add_vertex(pi);
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
          if (periodic)
          {
            if (dx > box->xhi) dx -= box->Lx;
            else if (dx < box->xlo) dx += box->Lx;
            if (dy > box->yhi) dy -= box->Ly;
            else if (dy < box->ylo) dy += box->Ly;
            if (dz > box->zhi) dz -= box->Lz;
            else if (dz < box->zlo) dz += box->Lz;
          }
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
  //m_msg->msg(Messenger::INFO, "Rebuilt neighbour list (using cell list).");
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
  //Mesh& mesh = m_system->get_mesh();
  int N = m_system->size();
  double dist = m_contact_dist;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    double ri = pi.get_radius();
    //cout << i << " --> ";
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
**/
bool NeighbourList::build_mesh()
{
    
  Mesh& mesh = m_system->get_mesh();
  int N = m_system->size();
  
  for (int i = 0; i < N; i++)
    for (unsigned int j = 0; j < m_contact_list[i].size(); j++)
      mesh.add_edge(i,m_contact_list[i][j]);
  
  mesh.set_max_face_perim(m_max_perim);
  mesh.generate_faces();
  mesh.generate_dual_mesh();
  mesh.postprocess();
  m_system->update_mesh();
  
  return true;
   
}

#ifdef HAS_CGAL
/*! Use CGAL to build 2D triangulation. This can be done only in 
 *  the plane, so this function will check if all z-coordinates are zero
 *  before proceeding. Otherwise, it will raise an exception.
*/
bool NeighbourList::build_triangulation()
{
  vector< pair<Point,unsigned> > points;
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

  for(Delaunay::Finite_faces_iterator fit = triangulation.finite_faces_begin(); fit != triangulation.finite_faces_end(); fit++)
  {
    Delaunay::Face_handle face = fit;
    unsigned int i = face->vertex(0)->info();
    unsigned int j = face->vertex(1)->info();
    unsigned int k = face->vertex(2)->info();
    //cout << i << " " << j << " " << k << endl;
    Particle& pi = m_system->get_particle(i);
    Particle& pj = m_system->get_particle(j);
    Particle& pk = m_system->get_particle(k);
    double dx = pi.x - pj.x, dy = pi.y - pj.y, dz = pi.z - pj.z;
    m_system->apply_periodic(dx,dy,dz);
    if (std::sqrt(dx*dx + dy*dy + dz*dz) < m_max_edge_len)
    {
      if (find(m_contact_list[i].begin(),m_contact_list[i].end(),j) == m_contact_list[i].end()) m_contact_list[i].push_back(j);
      if (find(m_contact_list[j].begin(),m_contact_list[j].end(),i) == m_contact_list[j].end()) m_contact_list[j].push_back(i);
    }
    dx = pj.x - pk.x, dy = pj.y - pk.y, dz = pj.z - pk.z;
    m_system->apply_periodic(dx,dy,dz);
    if (std::sqrt(dx*dx + dy*dy + dz*dz) < m_max_edge_len)
    {
      if (find(m_contact_list[j].begin(),m_contact_list[j].end(),k) == m_contact_list[j].end()) m_contact_list[j].push_back(k);
      if (find(m_contact_list[k].begin(),m_contact_list[k].end(),j) == m_contact_list[k].end()) m_contact_list[k].push_back(j);
    }
    dx = pk.x - pi.x, dy = pk.y - pi.y, dz = pk.z - pi.z;
    m_system->apply_periodic(dx,dy,dz);
    if (std::sqrt(dx*dx + dy*dy + dz*dz) < m_max_edge_len)
    {
      if (find(m_contact_list[k].begin(),m_contact_list[k].end(),i) == m_contact_list[k].end()) m_contact_list[k].push_back(i);
      if (find(m_contact_list[i].begin(),m_contact_list[i].end(),k) == m_contact_list[i].end()) m_contact_list[i].push_back(k);
    }
  }
 
  this->remove_dangling();
  
  return true;
  
}

//! Remove all edges that have one contanct
void NeighbourList::remove_dangling()
{
  bool done = false;
  int N = m_system->size();
  while(!done)
  {
    done = true;
    for (int i = 0; i < N; i++)
    {
      if (m_contact_list[i].size() == 1)
      {
        vector<int>& v = m_contact_list[m_contact_list[i][0]];
        vector<int>::iterator it = find(v.begin(),v.end(),i);
        if (it != v.end()) v.erase(it);
        m_contact_list[i].clear();
        done = false;
      }
    }
  }
}

#endif

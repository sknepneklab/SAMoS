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
 * \file neighbour_list.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Definition of the NeighbourList class
*/

#include "neighbour_list.hpp"

/* Get neighbour opposite to an edge
 * \note This is extremely inelegant. However, I cannot figure out
 * a simple way to tall CGAL directly to get me the other vertex 
 * opposite to an edge.
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

// Following two functions are used to determine which trainges are inside the boundary. 
// These are adopted from: http://doc.cgal.org/latest/Triangulation_2/index.html#title29
// Example: Triangulating a Polygonal Domain
static void mark_domains(Delaunay& ct, Delaunay::Face_handle start, int index, list<Delaunay::Edge>& border )
{
  if(start->info().nesting_level != -1)
  {
    return;
  }
  list<Delaunay::Face_handle> queue;
  queue.push_back(start);
  while(! queue.empty())
  {
    Delaunay::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1)
    {
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++)
      {
        Delaunay::Edge e(fh,i);
        Delaunay::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1)
        {
          if(ct.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}
//explore set of facets connected with non constrained edges,
//and attribute to each such set a nesting level.
//We start from facets incident to the infinite vertex, with a nesting
//level of 0. Then we recursively consider the non-explored facets incident 
//to constrained edges bounding the former set and increase the nesting level by 1.
//Facets in the domain are those with an odd nesting level.
void mark_domains(Delaunay& cdt)
{
  for(Delaunay::All_faces_iterator it = cdt.all_faces_begin(); it != cdt.all_faces_end(); ++it)
  {
    it->info().nesting_level = -1;
  }
  list<Delaunay::Edge> border;
  mark_domains(cdt, cdt.infinite_face(), 0, border);
  while(!border.empty())
  {
    Delaunay::Edge e = border.front();
    border.pop_front();
    Delaunay::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1)
    {
      mark_domains(cdt, n, e.first->info().nesting_level+1, border);
    }
  }
}
// End for functions imported from CGAL documentation. 

/*! Build neighbour list.
*/
void NeighbourList::build()
{
 m_list.clear();
 
 if (m_remove_detached)
   this->remove_detached();
  
 for (int i = 0; i < m_system->size(); i++)
   m_list.push_back(vector<int>());
 
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
#ifdef HAS_CGAL
  if (m_triangulation)
  {
    if (!m_system->has_boundary_neighbours())
    {
      m_msg->msg(Messenger::INFO,"Boundary neighbours have to be defined in a tissue simulations. Please use command \"read_cell_boundary\" in the config file.");
      throw runtime_error("Boundary neighbours not defined.");
    }
    m_contact_list.clear();
    for (int i = 0; i < m_system->size(); i++)
      m_contact_list.push_back(vector<int>());
    this->build_triangulation();
    this->build_faces(true);
  }
#endif
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

  mesh.set_circumcenter(m_circumcenter);
  mesh.set_max_face_perim(m_max_perim);
  mesh.generate_faces();
  mesh.generate_dual_mesh();
  mesh.postprocess(flag);
  m_system->update_mesh();
  
}

#ifdef HAS_CGAL
/*! Use CGAL to build 2D triangulation. This can be done only in 
 *  the plane, so this function will check if all z-coordinates are zero
 *  before proceeding. Otherwise, it will raise an exception.
*/
bool NeighbourList::build_triangulation()
{
  vector< pair<Point,unsigned> > points;
  vector< Point > points_for_boundary;  // need to double the data structure do the the way CGAL constrained triangulations work
  vector< pair<int,int> > boundary_indices;  // pairs of edges at the boundary 

  int N = m_system->size();
  vector<int> tissue_id;   // Bookkeeping vector for building triangulation
  vector<int> index_map(N,-1);
  int tissue_count = 0;
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.in_tissue && pi.z != 0.0) 
    {
      m_msg->msg(Messenger::ERROR,"Delaunay triangulation is only supported in plane. All z components of tissue particles must be set to zero.");
      throw runtime_error("Unable to build Delaunay triangulation for non-planar systems.");
    }
    if (pi.in_tissue)
    {
      //points.push_back( make_pair( Point(pi.x,pi.y), pi.get_id() ) );
      points.push_back( make_pair( Point(pi.x,pi.y), tissue_count ) );
      points_for_boundary.push_back(Point(pi.x,pi.y));
      tissue_id.push_back(pi.get_id());
      index_map[pi.get_id()] = tissue_count;
      tissue_count++;
    }
  }
  for (int i = 0; i < N; i++)
  {
    Particle& pi = m_system->get_particle(i);
    if (pi.boundary && pi.in_tissue)
    {
      boundary_indices.push_back( make_pair(index_map[pi.get_id()],index_map[pi.boundary_neigh[0]]) );
      boundary_indices.push_back( make_pair(index_map[pi.get_id()],index_map[pi.boundary_neigh[1]]) );
    }
  }
  
  Delaunay triangulation;  
  triangulation.insert_constraints(points_for_boundary.begin(), points_for_boundary.end(), boundary_indices.begin(), boundary_indices.end());
  triangulation.insert(points.begin(),points.end());
  mark_domains(triangulation);

  for(Delaunay::Finite_edges_iterator eit = triangulation.finite_edges_begin(); eit != triangulation.finite_edges_end(); eit++)
  {
    Delaunay::Face_handle f1 = eit->first;
    Delaunay::Face_handle f2 = f1->neighbor(eit->second);
    int i = tissue_id[f1->vertex(f1->cw(eit->second))->info()];
    int j = tissue_id[f1->vertex(f1->ccw(eit->second))->info()];
    int k, l;
    int to_flip;
    bool can_add = false;
    bool simple_add = true;
    // if f1 or f2 are boundary handle this first
    if (triangulation.is_infinite(f1) || triangulation.is_infinite(f2))
    {
      if (triangulation.is_infinite(f2))
        k = tissue_id[f1->vertex(eit->second)->info()];                //  k is the index of the index opposite to the edge (i,j)
      else
        k = tissue_id[get_opposite(f2, f1->vertex(f1->cw(eit->second))->info(), f1->vertex(f1->ccw(eit->second))->info())];
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
      //  k is the index of the index opposite to the edge (i,j) for face f1
      k = tissue_id[f1->vertex(eit->second)->info()];                
      //  l is the index of the index opposite to the edge (i,j) for face f2
      l = tissue_id[get_opposite(f2, f1->vertex(f1->cw(eit->second))->info(), f1->vertex(f1->ccw(eit->second))->info())];   
      Particle& pi = m_system->get_particle(i);   
      Particle& pj = m_system->get_particle(j);  
      Particle& pk = m_system->get_particle(k);  
      Particle& pl = m_system->get_particle(l);
      bool all_f1_boundary = pi.boundary && pj.boundary && pk.boundary;
      bool all_f2_boundary = pi.boundary && pj.boundary && pl.boundary;
      can_add = !(all_f1_boundary && all_f2_boundary) || f1->info().in_domain() || f2->info().in_domain();
      if (can_add)
      {
        if (find(pi.boundary_neigh.begin(),pi.boundary_neigh.end(),pj.get_id()) != pi.boundary_neigh.end()) // Flip only if edge (i,j) is a boundary edge
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
    }
    if (m_static_boundary)  // In case boundary is static, ignore all additions of new vertices
      simple_add = true;
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
        Particle& p1 = m_system->get_particle(i1);  // boundary particle
        Particle& p2 = m_system->get_particle(i2);  // boundary particle
        Particle& p3 = m_system->get_particle(i3);  // internal particle

        double x, y, z;                  // contains coordinates of mirrored particles 
        mirror(p3, p1, p2, x, y, z);     // compute position of mirrored particle 
        Particle p(m_system->size(), p1.get_type(), p1.get_radius());    // generate new particle with the "last" id and inhereted type and radius from p1
        i4 = p.get_id();                                            
        // set parameters for the new particle
        p.x = x; p.y = y; p.z = z;
        p.Nx = p3.Nx;  p.Ny = p3.Ny;  p.Nz = p3.Nz;
        p.nx = p3.nx;  p.ny = p3.ny;  p.nz = p3.nz;
        p.vx = 0.5*(p1.vx+p2.vx);  p.vy = 0.5*(p1.vy+p2.vy);  p.vz = 0.5*(p1.vz+p2.vz);
        p.coordination = 0;
        // copy all the groups of p1 (we need to do this as internal particles might be in the different group than the boundary)
        for(list<string>::iterator it_g = p1.groups.begin(); it_g != p1.groups.end(); it_g++)
          p.add_group(*it_g);
        // also make sure that the new particle belongs to the boundary group
        if (m_system->has_group("boundary"))
          p.add_group("boundary");
        // Note: Make sure that new particle is added to all necessary groups. 
        p.boundary = true;
        p.in_tissue = true;
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
 *  rest of the tissue. This is clearly only possible if triangulation has been set, i.e., for
 *  tissue simulations. 
 */
 void NeighbourList::remove_detached()
 {
   Mesh& mesh = m_system->get_mesh();
   if (mesh.size() == 0) return; 
   

   vector<int> to_remove;
   int offset = 0;  // We need to shift vertex ids to match them in the removal
   for (unsigned int i = 0; i < m_contact_list.size(); i++)
   {
     Particle& pi = m_system->get_particle(i);
     if (pi.in_tissue && m_contact_list[i].size() == 0)
     {
       to_remove.push_back(i-offset);
       offset++;       
     }
   }

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


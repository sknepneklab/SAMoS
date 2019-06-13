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
 * \file group.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 25-Sept-2014
 * \brief Declaration of Group class.
 */ 

#ifndef __GROUP_HPP__
#define __GROUP_HPP__

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <memory>


using std::ostream;
using std::vector;
using std::string;
using std::find;
using std::shared_ptr;

/*! Group class
 *  This class defines groups of particles 
 *  that have same behaviour under integration. For example, 
 *  they can all be integrated with the NVE integrator
 */
class Group
{
public:
  
  Group() { }
  
  //! Construct a Particle object
  //! \param id group id
  //! \param name group name
  Group(int id, const string name) : m_id(id), m_name(name) 
  {
    m_size = 0;
  }
  
  //! Get group id
  int get_id() const { return m_id; } //!< \return group id (m_id)
  
  //! Get group name
  string get_name() const { return m_name; } //!< \return group name
  
  //! Get size of the group (number of particles)
  int get_size() const { return m_particles.size(); } //!< \return group size
  
  //! Add particle to a group
  //! \param id particle id to add
  void add_particle(int id) 
  { 
    // protect against same particle being added two or more times to a group
    if (find(m_particles.begin(),m_particles.end(),id) == m_particles.end())
    {
      m_particles.push_back(id);
      m_size++;
    }
  }
  
  //! Remove particle from group
  //! \param id particle id to add
  void remove_particle(int id) 
  { 
    vector<int>::iterator it_p = find(m_particles.begin(), m_particles.end(), id);
    if (it_p != m_particles.end())
    {
      m_particles.erase(it_p);
      m_size--;
    }
  }
  
  
  //! Shift indices
  //! When a particle is removed from the system, we need to 
  //! shift down all indices that are larger than its original value
  //! \param id if of the removed particle
  void shift(int id)
  {
    vector<int>::iterator it_e = find(m_particles.begin(), m_particles.end(), id);
    if (it_e != m_particles.end())
    {
      m_particles.erase(it_e);
      m_size--;
    }
    for (unsigned int i = 0; i < m_particles.size(); i++)
      if (m_particles[i] > id)
        m_particles[i]--;
  }
  
  //! Get particles in the group
  vector<int>& get_particles() { return m_particles; } //!< \return reference to the vector containing indices of all particles in this group
    
private:  // Make these attributes immutable 
  
  int m_id;                        //!< Unique id
  string m_name;                   //!< Name of the group
  int m_size;                      //!< Number of particles in the group 
  vector<int> m_particles;         //!< Contains indices of all particles in the group
  
};

typedef shared_ptr<Group> GroupPtr;

#endif

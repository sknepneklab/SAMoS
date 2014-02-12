/* *************************************************************
 *  
 *   Active Particles on Curved Spaces (APCS)
 *   
 *   Author: Rastko Sknepnek
 *  
 *   Division of Physics
 *   School of Engineering, Physics and Mathematics
 *   University of Dundee
 *   
 *   (c) 2013
 *   
 *   This program cannot be used, copied, or modified without
 *   explicit permission of the author.
 * 
 * ************************************************************* */

/*!
 * \file neighbour_list.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 17-Oct-2013
 * \brief Declaration of NeighbourList class.
 */ 

#ifndef __NEIGHBOUR_LIST_HPP__
#define __NEIGHBOUR_LIST_HPP__

#include <vector>
#include <stdexcept>

#include "messenger.hpp"
#include "system.hpp"
#include "cell_list.hpp"

using std::vector;


/*! Auxiliary structure to keep coordinates of the first state 
 *  after the build so we can check if it is necessary to rebuild.
*/
struct PartPos
{
  //! Constructor
  //! \param x x coordinate of the stored position 
  //! \param y y coordinate of the stored position 
  //! \param z z coordinate of the stored position 
  PartPos(double x, double y, double z) : x(x), y(y), z(z) { }
  //@{
  double x, y, z;  //!< Stores actual position 
  //@}
};


/*! This class handles neighbour lists for fast potential and force 
 *  calculations. It is implemented as a vector of vectors (for performance). 
 *  Since the system is on a curved surface
 *  it would be hard to make a general cell list. Therefore, we relay on
 *  the slower but generic N^2 list generation 
*/
class NeighbourList
{
public:
  
  //! Construct NeighbourList object
  //! \param sys Reference to the System object
  //! \param msg Constant reference to the Messenger object
  //! \param cutoff Cutoff distance (should be set to potential cutoff distance + padding distance)
  //! \param pad Padding distance
  NeighbourList(SystemPtr sys, MessengerPtr msg, double cutoff, double pad) : m_system(sys), m_msg(msg), m_cut(cutoff), m_pad(pad)
  {
    for (int i = 0; i < m_system->size(); i++)
      m_list.push_back(vector<int>());
    // Check if box is large enough for cell list
    if (m_system->get_box()->Lx > 2.0*(cutoff+pad) && m_system->get_box()->Ly > 2.0*(cutoff+pad) && m_system->get_box()->Lz > 2.0*(cutoff+pad))
    {
      m_use_cell_list = true;
      m_cell_list = boost::shared_ptr<CellList>(new CellList(m_system,m_msg,cutoff+pad));
      m_msg->msg(Messenger::INFO,"Using cell lists for neighbour list builds.");
    }
    else
    {
      m_use_cell_list = false;
      m_msg->msg(Messenger::INFO,"Box dimensions are too small to be able to use cell lists. Neighbour list will be built using N^2 algorithm.");
    }
    this->build();
  }
  
  //! Destructor
  ~NeighbourList()
  {
    for(unsigned int i = 0; i < m_list.size(); i++)
      m_list[i].clear();
    m_list.clear();
    m_old_state.clear();
  }
  
  //! Check is neighbour list of the given particle needs update
  //! \param p particle to check 
  //! \return true if the list needs update
  bool need_update(Particle& p)
  {
    int id = p.get_id();
    double dx = m_old_state[id].x - p.x;
    double dy = m_old_state[id].y - p.y;
    double dz = m_old_state[id].z - p.z;
    if (dx*dx + dy*dy + dz*dz < 0.25*m_pad*m_pad)
      return false;
    else
      return true;
  }
  
  //! Get neighbour list for a give particle
  //! \param id Particle id
  //! \return Reference to the particle's neighbour list
  vector<int>& get_neighbours(int id) { return m_list[id]; }
  
  //! Get neighbour list cutoff distance
  double get_cutoff() { return m_cut;  }  //!< \return neighbour list cutoff distance
  
  //! Build neighbour list
  void build();
  
private:
  
  SystemPtr m_system;              //!< Pointer to the System object
  MessengerPtr m_msg;              //!< Handles messages sent to output
  CellListPtr m_cell_list;         //!< Pointer to the cell list
  vector<vector<int> >  m_list;    //!< Holds the list for each particle 
  vector<PartPos> m_old_state;     //!< Coordinates of particles right after the build
  double m_cut;                    //!< List build cutoff distance 
  double m_pad;                    //!< Padding distance (m_cut should be set to potential cutoff + m_pad)
  bool m_use_cell_list;            //!< If true, use cell list to speed up neighbour list builds
  
  // Actual neighbour list builds
  void build_nsq();    //!< Build with N^2 algorithm
  void build_cell();   //!< Build using cells list
  
};

typedef shared_ptr<NeighbourList> NeighbourListPtr;

#endif
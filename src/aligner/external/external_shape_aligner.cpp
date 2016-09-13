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
 *   (c) 2015, 2016
 * 
 *   Author: Silke Henkes
 * 
 *   Department of Physics 
 *   Institute for Complex Systems and Mathematical Biology
 *   University of Aberdeen  
 * 
 *   (c) 2014, 2015, 2016
 *  
 *   This program cannot be used, copied, or modified without
 *   explicit written permission of the authors.
 * 
 * ************************************************************* */

/*!
 * \file external_shape_aligner.cpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 16-Oct-2014
 * \brief Declaration of ExternalShapeAlign class
 */ 


#include "external_shape_aligner.hpp"

void ExternalShapeAlign::compute()
{
  int N = m_system->size();
  double J = m_J;
  Mesh& mesh = m_system->get_mesh();

  if (mesh.size() > 0) 
  {
    vector<Vertex>& vertices = mesh.get_vertices();
    vector<Face>& faces = mesh.get_faces();
    double ax = 1.0, ay;  // components of the shape eigenvector
    for  (int i = 0; i < N; i++)
    {
      Particle& pi = m_system->get_particle(i);
      if (!pi.boundary)
      {
        Vertex& V = vertices[pi.get_id()];
        // Compute components of the shape tensor defined by outter product of cell edges
        double A = 0.0; // xx component
        double B = 0.0; // xy component 
        double C = 0.0; // yy component
        for (unsigned int d = 0; d < V.dual.size(); d++)
        {
          Face& f1 = faces[V.dual[d]];
          int d_p_1 = (d == V.dual.size()-1) ? 0 : d + 1;
          Face& f2 = faces[V.dual[d_p_1]];
          Vector3d dr = f2.rc - f1.rc;
          A += dr.x*dr.x;  B += dr.x*dr.y;   C += dr.y*dr.y;
        }
        A /= V.dual.size();  B /= V.dual.size();  C /= V.dual.size();
        // Maximum eigenvalue is 
        double l1 = 0.5*(A+C + sqrt((A-C)*(A-C)+4*B*B));
        double l2 = 0.5*(A+C - sqrt((A-C)*(A-C)+4*B*B));
        double lambda = (l1 > l2) ? l1 : l2;
        // compute corresponding eigen vector
        ay = -(A-lambda)/B;
        // normalise it
        double len_a = sqrt(ax*ax + ay*ay);
        ax /= len_a;  ay /= len_a;
        double tau_z = pi.nx*ay - pi.ny*ax;
        if (m_has_params)
          J = m_type_params[pi.get_type()].J;
        pi.tau_x += 0.0;
        pi.tau_y += 0.0;
        pi.tau_z += J*tau_z;  // per construction torque is only in the z direction
      }  
    }
  }
}

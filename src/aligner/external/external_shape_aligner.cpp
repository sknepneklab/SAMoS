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
        // compute centre of the face 
        Vertex& V = vertices[pi.get_id()];
        Vector3d rcm(0.0,0.0,0.0);
        for (int f = 0; f < V.n_faces; f++)
        {
          Face& face = faces[V.faces[f]];
          rcm += face.rc;
        }
        rcm.scale(1.0/V.n_faces);
        // Now compute components of the gyration tensor
        double A = 0.0; // xx component
        double B = 0.0; // xy component 
        double C = 0.0; // yy component
        for (int f = 0; f < V.n_faces; f++)
        {
          Face& face = faces[V.faces[f]];
          Vector3d r = face.rc - rcm;
          A += r.x*r.x;  B += r.x*r.y;  C += r.y*r.y;
        } 
        A /= V.n_faces;  B /= V.n_faces;  C /= V.n_faces;
        // Maximum eigenvalue is 
        double lambda = 0.5*(A+C + sqrt((A-C)*(A-C)+4*B*B));
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

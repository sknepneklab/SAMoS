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
 * \file samos.hpp
 * \author Rastko Sknepnek, sknepnek@gmail.com
 * \date 03-Mar-2015
 * \brief All include files for the main code.
*/

#ifndef __SAMOS_HPP__
#define __SAMOS_HPP__

#include "defaults.hpp"
#include "messenger.hpp"
#include "dump.hpp"
#include "logger.hpp"
#include "parse_command.hpp"
#include "parse_constraint.hpp"
#include "parse_external.hpp"
#include "parse_input.hpp"
#include "parse_run.hpp"
#include "parse_potential.hpp"
#include "parse_aux.hpp"
#include "parse_parameters.hpp"
#include "parse_rng_seed.hpp"
#include "parse_box.hpp"
#include "parse_integrator.hpp"
#include "parse_log_dump.hpp"
#include "parse_align.hpp"
#include "parse_external_align.hpp"
#include "parse_group.hpp"
#include "parse_disable.hpp"
#include "parse_population.hpp"
#include "parse_bond.hpp"
#include "parse_angle.hpp"
#include "parse_timestep.hpp"
#include "constraint.hpp"
#include "constraint_sphere.hpp"
#include "constraint_plane.hpp"
#include "constraint_plane_walls.hpp"
#include "constraint_cylinder.hpp"
#include "constraint_peanut.hpp"
#include "constraint_torus.hpp"
#include "constraint_ellipsoid.hpp"
#include "constraint_gyroid.hpp"
#include "constraint_actomyo.hpp"
#include "constraint_hourglass.hpp"
#include "constraint_gaussian_bump.hpp"
#include "constraint_none.hpp"
#include "constraint_tetrahedron.hpp"
#include "constraint_slab.hpp"
#include "rng.hpp"
#include "particle.hpp"
#include "vector3d.hpp"
#include "box.hpp"
#include "system.hpp"
#include "neighbour_list.hpp"
#include "external_potential.hpp"
#include "external_gravity_potential.hpp"
#include "external_harmonic_potential.hpp"
#include "external_self_propulsion.hpp"
#include "pair_potential.hpp"
#include "pair_coulomb_potential.hpp"
#include "pair_soft_potential.hpp"
#include "pair_lj_potential.hpp"
#include "pair_gaussian_potential.hpp"
#include "pair_morse_potential.hpp"
#include "pair_active_potential.hpp"
#include "pair_rod_potential.hpp"
#include "pair_ljrod_potential.hpp"
#include "pair_soft_attractive_potential.hpp"
#include "pair_vertex_particle_potential.hpp"
#include "pair_line_tension_potential.hpp"
#include "pair_boundary_bending_potential.hpp"
#include "pair_boundary_attraction_potential.hpp"
#include "potential.hpp"
#include "integrator.hpp"
#include "integrator_brownian.hpp"
#include "integrator_brownian_rod.hpp"
#include "integrator_vicsek.hpp"
#include "integrator_nve.hpp"
#include "integrator_nematic.hpp"
#include "integrator_actomyo.hpp"
#include "integrator_brownian_pos.hpp"
#include "integrator_brownian_align.hpp"
#include "integrator_langevin.hpp"
#include "aligner.hpp"
#include "pair_aligner.hpp"
#include "pair_polar_aligner.hpp"
#include "pair_nematic_aligner.hpp"
#include "pair_vicsek_aligner.hpp"
#include "external_aligner.hpp"
#include "external_ajpolar_aligner.hpp"
#include "external_field_aligner.hpp"
#include "external_ajnematic_aligner.hpp"
#include "population.hpp"
#include "population_random.hpp"
#include "population_density.hpp"
#include "population_growth.hpp"
#include "population_elongation.hpp"
#include "population_cell.hpp"
#include "bond_potential.hpp" 
#include "bond_harmonic_potential.hpp"
#include "bond_fene_potential.hpp"
#include "bond_active_force.hpp"
#include "angle_potential.hpp"
#include "angle_harmonic_potential.hpp"
#include "angle_cosine_potential.hpp"
#include "value.hpp"
#include "vertex.hpp"
#include "edge.hpp"
#include "face.hpp"
#include "mesh.hpp"

#endif
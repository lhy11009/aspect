/*
 Copyright (C) 2015 - 2022 by the authors of the ASPECT code.

 This file is part of ASPECT.

 ASPECT is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 ASPECT is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ASPECT; see the file LICENSE.  If not see
 <http://www.gnu.org/licenses/>.
 */

#ifndef _aspect_particle_generator_world_builder_feature_surface_h
#define _aspect_particle_generator_world_builder_feature_surface_h

#include <aspect/particle/generator/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      /**
       * Generate a distribution of particles that is determined by the
       * coordinates given in an ascii data file.
       *
       * @ingroup ParticleGenerators
       */
      template <int dim>
      class WorldBuilderFeatureSurface : public Interface<dim>
      {
        public:
          /**
           * Reads in a file and generate a set of particles at the prescribed
           * positions.
           *
           * @param [in,out] particle_handler The particle handler into which
           * the generated particles should be inserted.
           */
          void
          generate_particles(Particles::ParticleHandler<dim> &particle_handler) override;

          /**
           * The method of generating particles by looking at decreasing radius
           * */
          void
          generate_particles_by_radius(Particles::ParticleHandler<dim> &particle_handler, const double);

          /**
           * The method of generating particles by looking at a constant distance
           * */
          void
          generate_particles_by_distance(Particles::ParticleHandler<dim> &particle_handler, const double);

          // avoid -Woverloaded-virtual
          // TODO: remove this using directive once the following deprecated
          // function in the interface class has been removed:
          // generate_particles(std::multimap<Particles::internal::LevelInd, Particle<dim>> &particles)
          using Generator::Interface<dim>::generate_particles;

          /**
           * Declare the parameters this class takes through input files.
           */
          static
          void
          declare_parameters (ParameterHandler &prm);

          /**
           * Read the parameters this class declares from the parameter file.
           */
          void
          parse_parameters (ParameterHandler &prm) override;

        private:
          /**
           * Number of particles to create
           */
          types::particle_index n_particles;

          /**
           * Number of particles created
           */
          types::particle_index n_particles_created;

          /**
           * Number of particles created to maintain constant distance
           */
          types::particle_index n_insert_particles;

          /**
          * The name of the geometry model
          */
          std::string geometry_model_name;

          /**
           * The method to generate particles
          */
          int generate_method;

          /**
           * Number of sets of particles
          */
          int n_set_of_particles;

          /**
           * Number of particles created in each set
           */
          std::vector<int> n_particles_in_each_set;

          /**
          * position of the feature start
          */
          double feature_start;

          /**
          * The distance of the particles
          * to the feature surface
          */
          std::vector<double> feature_surface_distances;

          /**
          * The minimum distance of the particles
          * to the feature surface
          */
          double minimum_feature_surface_distance;

          /**
          * The maximum distance of the particles
          * to the feature surface
          */
          double maximum_feature_surface_distance;

          /**
          * The interval of the distance of the particles
          * to the feature surface
          */
          double interval_feature_surface_distance;

          /**
          * The Maximum radius to generate particles
          */
          double maximum_radius;

          /**
          * The Minimum radius to generate particles
          */
          double minimum_radius;

          /**
          * The coordinate of the query point
          * box geometry: x coordinate
          * chunk geometry: theta coordinate
          */
          double query_coordinate;

          /**
          * The length of search
          * box geometry: x coordinate
          * chunk geometry: theta coordinate
          */
          double search_length;

          /**
          * The maximum step of search
          */
          int search_max_step;

          /**
          * The start point of the search
          * box geometry: x coordinate
          * chunk geometry: theta coordinate
          */
          double search_start;

          /**
          * The tolerance of the searching
          */
          double search_tolerance;

          double search_distance;

          double search_angle;

          bool positive_dipping;
      };

    }
  }
}

#endif

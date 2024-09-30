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

#include <world_builder/world.h>
#include <aspect/utilities.h>
#include <aspect/particle/generator/world_builder_feature_surface.h>


namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      template <>
      void
      WorldBuilderFeatureSurface<2>::generate_particles(Particles::ParticleHandler<2> &particle_handler)
      {
        //todo_wb
        double x = 0.0;
        double y = 0.0;
        double x_last = 0.0;
        double y_last = 0.0;
        types::particle_index n_particles_on_plate = n_particles - n_particles_on_slab;
        types::particle_index particle_index = 0;
        /*
        First part: generate particles on the plate
        */
        // compute the interval of points on the plate for later steps
        double particle_interval_plate = 0.0;
        if (geometry_model_name == "chunk")
          {
            particle_interval_plate = maximum_radius * feature_start / n_particles_on_plate;
          }
        else
          {
            // other geometries are not implemented yet
            Assert (false, ExcNotImplemented());
          }

        for (unsigned i = 0; i < n_particles_on_plate; ++i)
          {
            const double depth = feature_surface_distance;
            const double radius = maximum_radius - depth;
            const double theta = feature_start * i / (n_particles_on_plate-1);
            if (geometry_model_name == "chunk")
              {
                x = radius*cos(theta);
                y = radius*sin(theta);
              }
            else
              {
                // other geometries are not implemented yet
                Assert (false, ExcNotImplemented());
              }
            const Point<2> particle_position(x, y);
            this->insert_particle_at_position(particle_position, particle_index, particle_handler);
            ++particle_index;
            if (i == n_particles_on_plate - 1)
              {
                x_last = x;
                y_last = y;
              }
          }
        /*
        Second part: generate particles on the slab.
        */
        // Initiate variables
        // The depth interval between adjacent particles
        const double interval = (maximum_radius - minimum_radius) / (n_particles_on_slab - 1.0); // 10 km
        unsigned n_insert_particles = 0;
        for (unsigned i = 0; i < n_particles_on_slab; ++i)
          {
            const double depth = i * interval;
            const double radius = maximum_radius - depth;
            if (geometry_model_name == "chunk")
              {
                x = radius*cos(search_start);
                y = radius*sin(search_start);
              }
            else
              {
                // other geometries are not implemented yet
                Assert (false, ExcNotImplemented());
              }
            // get the distance of the initial point to the surface of the slab
            double diff, min_diff, initial_distance, distance;
            const std::array<double, 3> point = {{x, y, 0.0}};
            auto plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
            initial_distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
            min_diff = (initial_distance - feature_surface_distance) / feature_surface_distance ;
            // earch for a point that approximates the desinated distance to the feature surface
            int search_step = 0;
            for (int j = -search_max_step; j < search_max_step; j++)
              {
                if (geometry_model_name == "chunk")
                  {
                    x = radius*cos(search_start + j*search_length);
                    y = radius*sin(search_start + j*search_length);
                  }
                else
                  {
                    // other geometries are not implemented yet
                    Assert (false, ExcNotImplemented());
                  }
                const std::array<double, 3> point = {{x, y, 0.0}};
                plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
                distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
                diff = (distance - feature_surface_distance)/feature_surface_distance;
                if (abs(diff) < abs(min_diff))
                  {
                    min_diff = diff;
                    search_step = j;
                  }
                // std::cerr << "i = " << i << ", j = "<< j << ", distance = " << distance << std::endl;  //debug
              }
            // first, check this point is valid (i.e. a non-infinite minimum difference)
            // then, pair it with an adjacent point (also check this point has a valid distance)
            // to create a range and iterate for an approximation of the desinated point.
            if (isinf(min_diff))
              {
                continue;
              }
            double search_result_0, search_result_1, search_result_temp;
            if (min_diff < 0.0)
              {
                search_result_0 = search_start + (search_step-1)*search_length;
                search_result_1 = search_start + search_step*search_length;
                if (geometry_model_name == "chunk")
                  {
                    x = radius*cos(search_result_0);
                    y = radius*sin(search_result_0);
                  }
                else
                  {
                    // other geometries are not implemented yet
                    Assert (false, ExcNotImplemented());
                  }
                const std::array<double, 3> point1 = {{x, y, 0.0}};
                auto plane_distances = this->get_world_builder().distance_to_plane(point1, depth, "Slab");
                distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
              }
            else
              {
                search_result_0 = search_start + search_step*search_length;
                search_result_1 = search_start + (search_step+1)*search_length;
                if (geometry_model_name == "chunk")
                  {
                    x = radius*cos(search_result_1);
                    y = radius*sin(search_result_1);
                  }
                else
                  {
                    // other geometries are not implemented yet
                    Assert (false, ExcNotImplemented());
                  }
                const std::array<double, 3> point1 = {{x, y, 0.0}};
                auto plane_distances = this->get_world_builder().distance_to_plane(point1, depth, "Slab");
                distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
              }
            //std::cerr << "distance: " << distance << std::endl;  // debug
            if (isinf(distance))
              {
                continue;
              }
            // Search for a point with a desinated distance to the surface
            diff = min_diff;
            bool found = false;
            while (true)
              {
                search_result_temp = (search_result_0 + search_result_1) / 2.0;
                if (geometry_model_name == "chunk")
                  {
                    x = radius*cos(search_result_temp);
                    y = radius*sin(search_result_temp);
                  }
                else
                  {
                    Assert (false, ExcNotImplemented());
                  }
                const std::array<double, 3> point_temp = {{x, y, 0.0}};
                plane_distances = this->get_world_builder().distance_to_plane(point_temp, depth, "Slab");
                distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
                diff = (distance - feature_surface_distance)/feature_surface_distance;
                //std::cerr << "search_result_temp: " << search_result_temp << ", diff: " << diff << std::endl;  // debug
                if (isinf(diff))
                  {
                    // Not lucky, it's possible to mess up with the curvature of the surface that
                    // the medium of two valid points may not be valid.
                    break;
                  }
                else if (abs(diff) < search_tolerance)
                  {
                    // Good, a point approximating the desinated distance to the surface is food
                    found = true;
                    break;
                  }
                else if (diff < 0.0)
                  {
                    // search in the direction of increasing coordinate
                    search_result_1 = search_result_temp;
                  }
                else
                  {
                    // search in the direction of decreasing coordinate
                    search_result_0 = search_result_temp;
                  }
              }
            // add a new particle at the point found in the last step
            if (found)
              {
                if (geometry_model_name == "chunk")
                  {
                    x = radius*cos(search_result_temp);
                    y = radius*sin(search_result_temp);
                  }
                else
                  {
                    // other geometries are not implemented yet
                    Assert (false, ExcNotImplemented());
                  }
                // make sure the distance between adjacent points are smaller than the distance between points
                // on the plate by inserting additional points.
                const double n_insert = ceil(pow((x - x_last) * (x - x_last) +\
                                                 (y - y_last) * (y - y_last), 0.5) / particle_interval_plate);
                for (unsigned j = 1; j < n_insert; ++j)
                  {
                    const double x_i = (x_last * (n_insert - j) + x * j) / n_insert  ;
                    const double y_i = (y_last * (n_insert - j) + y * j) / n_insert  ;
                    const Point<2> particle_position(x_i, y_i);
                    this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                    ++particle_index;
                    ++n_insert_particles;
                  }
                const Point<2> particle_position(x, y);
                this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                ++particle_index;
                x_last = x;
                y_last = y;
              }
          }
        particle_handler.update_cached_numbers();
        std::cout << "WorldBuilderFeatureSurface::generate_particles: " << n_particles_on_plate <<
                  " are generated on the plate; Expect " << n_particles_on_slab << " on the slab, "
                  << (particle_index - n_particles_on_plate) << " are successfully generated on the slab, including "
                  << n_insert_particles << " inserted separtately to maintain a regular distances between adjacent points."
                  << std::endl << std::endl;
      }

      template <>
      void
      WorldBuilderFeatureSurface<3>::generate_particles(Particles::ParticleHandler<3> &)
      {
        Assert (false, ExcNotImplemented());
      }

      template <int dim>
      void
      WorldBuilderFeatureSurface<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Geometry model");
        {
          prm.declare_entry ("Model name", "unspecified",  Patterns::Anything (), "Name of the geometry model");
        }
        prm.leave_subsection();

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            prm.declare_entry ("Number of particles", "1000",
                               Patterns::Double (0.),
                               "Total number of particles to create (not per processor or per element). "
                               "The number is parsed as a floating point number (so that one can "
                               "specify, for example, '1e4' particles) but it is interpreted as "
                               "an integer, of course.");

            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("World builder feature surface");
              {
                prm.declare_entry ("Number of particles on the slab", "1000",
                                   Patterns::Double (0.),
                                   "Total number of particles to generate(not per processor or per element) on the subducting slab surface"
                                   "The number is parsed as a floating point number (so that one can "
                                   "specify, for example, '1e4' particles) but it is interpreted as "
                                   "an integer, of course.");
                prm.declare_entry ("Feature surface distance", "0.0",
                                   Patterns::Double (0.),
                                   "The distance of particles to the feature surface. Units: \\si{\\meter}.");
                prm.declare_entry ("Maximum radius", "1.0",
                                   Patterns::Double (0.),
                                   "Radius at the shallow limit of particles. Units: \\si{\\meter}.");
                prm.declare_entry ("Minimum radius", "1.0",
                                   Patterns::Double (0.),
                                   "Radius at the lower limit of particles. Units: \\si{\\meter}.");
                prm.declare_entry ("Query point coordinate", "0.0",
                                   Patterns::Double (0.),
                                   "The coordinate of the query point. If the geometry is \"box\", then"
                                   " The x coordinate of the point should be entered."
                                   " If the geometry is \"chunk\", then the theta coordinate"
                                   " of the point should be entered. Units: \\si{\\meter} or {\\radian}.");
                prm.declare_entry ("Feature start", "1.0",
                                   Patterns::Double (0.),
                                   "The position where the feature starts,"
                                   " this is given as the first coordinate of the point (x or theta)."
                                   " Units: \\si{\\meter} or {\\radian}.");
                prm.declare_entry ("Search length", "1.0",
                                   Patterns::Double (0.),
                                   "The length of each step of searching. When a step of searching is performed,"
                                   " The query point is altered by this length along the first coordinate (x or theta)."
                                   " Units: \\si{\\meter} or {\\radian}.");
                prm.declare_entry ("Search max step", "10",
                                   Patterns::Integer(),
                                   "The total steps of searching. The code searches around the Query point by this number of steps");
                prm.declare_entry ("Search start", "1.0",
                                   Patterns::Double (0.),
                                   "The position of the query point to start the search,"
                                   " this is given as the first coordinate of the point (x or theta)."
                                   " Units: \\si{\\meter} or {\\radian}.");
                prm.declare_entry ("Search tolerance", "0.01",
                                   Patterns::Double (0.),
                                   "The tolerance of the searching."
                                   " this is the tolerance on the relative difference between the query points' distance"
                                   " to the surface and the desinated distance");
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }


      template <int dim>
      void
      WorldBuilderFeatureSurface<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Geometry model");
        {
          geometry_model_name = prm.get("Model name");
          AssertThrow((geometry_model_name == "box" || geometry_model_name == "chunk"),
                      ExcMessage("You need to select a Geometry model of"
                                 "either 'box' or 'chunk' to use the particle generator of"
                                 "world builder feature surface"));
        }
        prm.leave_subsection();

        prm.enter_subsection("Postprocess");
        {
          prm.enter_subsection("Particles");
          {
            n_particles    = static_cast<types::particle_index>(prm.get_double ("Number of particles"));
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("World builder feature surface");
              {
                n_particles_on_slab    = static_cast<types::particle_index>(prm.get_double ("Number of particles on the slab"));
                feature_surface_distance = prm.get_double("Feature surface distance");
                maximum_radius = prm.get_double("Maximum radius");
                minimum_radius = prm.get_double("Minimum radius");
                query_coordinate = prm.get_double("Query point coordinate");
                feature_start = prm.get_double("Feature start");
                search_start = prm.get_double("Search start");
                search_length = prm.get_double("Search length");
                search_max_step = prm.get_integer("Search max step");
                search_tolerance = prm.get_double("Search tolerance");
                AssertThrow(maximum_radius > minimum_radius,
                            ExcMessage("The value for the 'Maximum radius' has to be bigger than"
                                       " the value for the 'Minimum radius"));
              }
              prm.leave_subsection();
            }
            prm.leave_subsection();
          }
          prm.leave_subsection();
        }
        prm.leave_subsection();
      }
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Generator
    {
      ASPECT_REGISTER_PARTICLE_GENERATOR(WorldBuilderFeatureSurface,
                                         "world builder feature surface",
                                         "Generates a distribution of particles from coordinates "
                                         "specified along a feature surface defined in WorldBuilder."
                                         "All of the values that define this generator are read "
                                         "from a section ``Postprocess/Particles/Generator/World builder feature surface'' in the "
                                         "input file")
    }
  }
}

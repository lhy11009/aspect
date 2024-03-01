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
      WorldBuilderFeatureSurface<2>::generate_particles_by_radius(Particles::ParticleHandler<2> &particle_handler, const double feature_surface_distance)
      {
        double x = 0.0;
        double y = 0.0;
        types::particle_index particle_index = n_particles_created;
        // generate particles on the slab.
        // Initiate variables
        // The depth interval between adjacent particles
        const double interval = (maximum_radius - minimum_radius) / (n_particles - 1.0); // 10 km
        for (unsigned i = 0; i < n_particles; ++i)
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
            double diff = 0.0;
            double distance = 0.0;
            const std::array<double, 3> point = {{x, y, 0.0}};
            auto plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
            double initial_distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
            double min_diff = (initial_distance - feature_surface_distance) / feature_surface_distance ;
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
            double search_result_0 = 0.0;
            double search_result_1 = 0.0;
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
            double search_result_temp = 0.0;
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
                const Point<2> particle_position(x, y);
                this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                ++particle_index;
              }
          }
        particle_handler.update_cached_numbers();
        std::cout << "WorldBuilderFeatureSurface::generate_particles: distance = " << feature_surface_distance << "." << std::endl
                  << "assigned " << n_particles << " on the slab, "
                  << particle_index - n_particles_created  << " are successfully generated on the slab, including "
                  << n_insert_particles << " inserted separtately to maintain a regular distances between adjacent points."
                  << std::endl << std::endl;
        n_particles_created = particle_index;
      }

      template <>
      void
      WorldBuilderFeatureSurface<2>::generate_particles_by_distance(Particles::ParticleHandler<2> &particle_handler, const double feature_surface_distance)
      {
        //part I: get the initial point
        double x = 0.0;
        double y = 0.0;
        double x_last = 0.0;
        double y_last = 0.0;
        double theta_last = search_start;
        types::particle_index particle_index = n_particles_created;
        // generate particles on the slab.
        // Initiate variables
        // The depth interval between adjacent particles
        // The strategy is to decrease from the maximum radius
        // and look up with the new radius
        const double interval = 1e3; // 10 km
        int initial_search_steps = 1000;
        for (int i=0; i < initial_search_steps; ++i)
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
            double diff = 0.0;
            double distance = 0.0;
            const std::array<double, 3> point = {{x, y, 0.0}};
            auto plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
            double initial_distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
            double min_diff = (initial_distance - feature_surface_distance) / feature_surface_distance ;
            // search for a point that approximates the desinated distance to the feature surface
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
            double search_result_0 = 0.0;
            double search_result_1 = 0.0;
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
            double search_result_temp = 0.0;
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
                const Point<2> particle_position(x, y);
                this->insert_particle_at_position(particle_position, particle_index, particle_handler);
                ++particle_index;
                x_last = x;
                y_last = y;
                theta_last = search_result_temp;
                break;
              }
          }

        // second part: search for the following particles
        bool reach_end = false;
        std::string ending_scheme = std::string("");
        while (true)
          {
            // first get two search point
            double search_angle_0 = 0.0;
            double search_angle_1 = 0.0;
            if (positive_dipping)
              {
                // dipping towards positive azimuth
                search_angle_0 = theta_last + 1.0 / 2.0 * M_PI;
                search_angle_1 = theta_last + M_PI;
              }
            else
              {
                // dipping towards negative azimuth
                search_angle_0 = theta_last + M_PI;
                search_angle_1 = theta_last + 3.0 / 2.0 * M_PI;
              }
            double depth_last = maximum_radius - sqrt(x_last*x_last + y_last*y_last);
            double x0 = x_last + search_distance * cos(search_angle_0);
            double y0 = y_last + search_distance * sin(search_angle_0);
            const std::array<double, 3> point0 = {{x0, y0, 0.0}};
            double depth0 = maximum_radius - sqrt(x0*x0 + y0*y0);
            auto plane_distances_0 = this->get_world_builder().distance_to_plane(point0, depth0, "Slab");
            double distance_0 = plane_distances_0.get_distance_from_surface(); // the distance of the query point to the feature surface
            double diff_0 = (distance_0 - feature_surface_distance)/feature_surface_distance;
            double x1 = x_last + search_distance * cos(search_angle_1);
            double y1 = y_last + search_distance * sin(search_angle_1);
            const std::array<double, 3> point1 = {{x1, y1, 0.0}};
            double depth1 = maximum_radius - sqrt(x1*x1 + y1*y1);
            auto plane_distances_1 = this->get_world_builder().distance_to_plane(point1, depth1, "Slab");
            double distance_1 = plane_distances_1.get_distance_from_surface(); // the distance of the query point to the feature surface
            double diff_1 = (distance_1 - feature_surface_distance)/feature_surface_distance;
            std::vector<double> diff_history(search_max_step, 0.0);
            // then use the bisect method to look for the right point
            int j = 0;
            for (; j < search_max_step; ++j)
              {
                search_angle = (search_angle_0 + search_angle_1) / 2.0;
                x = x_last + search_distance * cos(search_angle);
                y = y_last + search_distance * sin(search_angle);
                const std::array<double, 3> point = {{x, y, 0.0}};
                double depth = maximum_radius - sqrt(x*x + y*y);
                auto plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
                double distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface
                double diff = (distance - feature_surface_distance)/feature_surface_distance;
                diff_history[j] = diff;
                if (abs(diff) < search_tolerance)
                  {
                    // found
                    break;
                  }
                else if (abs(diff) == std::numeric_limits<double>::infinity())
                  {
                    // reach unreasonable range, end is reached
                    ending_scheme = std::string("infinity value");
                    reach_end = true;
                    break;
                  }

                if (diff_0 * diff < 0)
                  {
                    search_angle_1 = search_angle;
                    diff_1 = diff;
                  }
                else if (diff_1 * diff < 0)
                  {
                    search_angle_0 = search_angle;
                    diff_0 = diff;
                  }
                else
                  {
                    // no result in the range of search, end is reached
                    ending_scheme = std::string("invalid searching range");
                    reach_end = true;
                    break;
                  }
              }
            // check if the end is reached
            if (reach_end)
              break;
            // break if an error happens, comment this to do the assert throw instead
            /*
            if (j >= search_max_step)
              {
                std::string diff_history_string = std::string("");
                for (int jj = 0; jj < j; ++jj)
                  {
                    diff_history_string += (std::to_string(jj) + ": " + std::to_string(diff_history[jj]) + "\n");
                  }
                std::cout << "The max search step is reached when looking for a particle at the WorldBuilder feature surface. "
                          "particle_index: " + std::to_string(particle_index) + ", theta_last: " + std::to_string(theta_last)
                          + ", depth_last: " + std::to_string(depth_last) + "." << std::endl;
                std::cout << diff_history_string << std::endl;
                break;
              }
            */
            // check the maximum searching step is not reached
            AssertThrow(j < search_max_step,
                        ExcMessage ("The max search step is reached when looking for a particle at the WorldBuilder feature surface. "
                                    "particle_index: " + std::to_string(particle_index) + ", theta_last: " + std::to_string(theta_last)
                                    + ", depth_last: " + std::to_string(depth_last) + "."));
            const Point<2> particle_position(x, y);
            this->insert_particle_at_position(particle_position, particle_index, particle_handler);
            ++particle_index;
            x_last = x;
            y_last = y;
          }

        particle_handler.update_cached_numbers();
        std::cout << "WorldBuilderFeatureSurface::generate_particles: distance = " << feature_surface_distance << "." << std::endl
                  << "Generating method = " << ((generate_method == 0)? "along radius": "constant distance") << std::endl
                  << "Assigned " << n_particles << " on the slab, "
                  << particle_index - n_particles_created << " are successfully generated on the slab, including "
                  << n_insert_particles << " inserted separtately to maintain a regular distances between adjacent points." << std::endl
                  << "ending scheme = " << ending_scheme
                  << std::endl << std::endl;
        n_particles_created = particle_index;
      }

      template <>
      void
      WorldBuilderFeatureSurface<2>::generate_particles(Particles::ParticleHandler<2> &particle_handler)
      {
        if (generate_method == 0)
          {
            for (int i_set = 0; i_set < n_set_of_particles ; ++i_set)
              {
                generate_particles_by_radius(particle_handler, feature_surface_distances[i_set]);
              }
          }
        else if (generate_method == 1)
          {
            for (int i_set = 0; i_set < n_set_of_particles ; ++i_set)
              {
                generate_particles_by_distance(particle_handler, feature_surface_distances[i_set]);
              }
          }
        else
          AssertThrow(false, ExcNotImplemented());
      }

      template <>
      void
      WorldBuilderFeatureSurface<3>::generate_particles(Particles::ParticleHandler<3> &)
      {
        Assert (false, ExcNotImplemented());
      }

      template <>
      void
      WorldBuilderFeatureSurface<3>::generate_particles_by_radius(Particles::ParticleHandler<3> &, const double)
      {
        Assert (false, ExcNotImplemented());
      }

      template <>
      void
      WorldBuilderFeatureSurface<3>::generate_particles_by_distance(Particles::ParticleHandler<3> &, const double)
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
                prm.declare_entry ("Generate method", "0",
                                   Patterns::Integer (0.), "The method to generate particles.");
                prm.declare_entry ("Feature surface distances", "0.0",
                                   Patterns::List(Patterns::Double (0.)),
                                   "The distances of particles to the feature surface. Each value is a separate set of particles."
                                   "Units: \\si{\\meter}.");
                prm.declare_entry ("Minimum feature surface distance", "0.0",
                                   Patterns::Double (0.),
                                   "The minimum distance of particles to the feature surface."
                                   "Units: \\si{\\meter}.");
                prm.declare_entry ("Maximum feature surface distance", "0.0",
                                   Patterns::Double (0.),
                                   "The maximum distance of particles to the feature surface."
                                   "Units: \\si{\\meter}.");
                prm.declare_entry ("Interval feature surface distance", "0.0",
                                   Patterns::Double (0.),
                                   "The interval of distance of particles to the feature surface."
                                   "Units: \\si{\\meter}.");
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
                prm.declare_entry ("Search distance", "5e3",
                                   Patterns::Double (0.),
                                   "Distance of searching.");
                prm.declare_entry ("Positive dipping", "true",
                                   Patterns::Bool(),
                                   "Whether the slab dips toward positive azimuth");
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
            n_particles_created = static_cast<types::particle_index>(0.0);
            n_insert_particles = static_cast<types::particle_index>(0.0);
            prm.enter_subsection("Generator");
            {
              prm.enter_subsection("World builder feature surface");
              {
                generate_method = prm.get_integer("Generate method");

                std::vector<double> feature_surface_distances_input = Utilities::string_to_double(Utilities::split_string_list(prm.get("Feature surface distances")));
                minimum_feature_surface_distance = prm.get_double("Minimum feature surface distance");
                maximum_feature_surface_distance = prm.get_double("Maximum feature surface distance");
                interval_feature_surface_distance = prm.get_double("Interval feature surface distance");
                if (feature_surface_distances_input.size() == 1 && feature_surface_distances_input[0] == 0.0)
                  {
                    n_set_of_particles = std::ceil((maximum_feature_surface_distance - minimum_feature_surface_distance) / interval_feature_surface_distance);
                    feature_surface_distances = std::vector<double> (n_set_of_particles, 0.0);
                    for (int i = 0; i < n_set_of_particles; ++i)
                      {
                        feature_surface_distances[i] = minimum_feature_surface_distance + interval_feature_surface_distance * i;
                      }
                  }
                else
                  {
                    feature_surface_distances = feature_surface_distances_input;
                    n_set_of_particles = feature_surface_distances.size();
                  }

                n_particles_in_each_set = std::vector<int> (n_set_of_particles, 0);
                maximum_radius = prm.get_double("Maximum radius");
                minimum_radius = prm.get_double("Minimum radius");
                query_coordinate = prm.get_double("Query point coordinate");
                feature_start = prm.get_double("Feature start");
                search_start = prm.get_double("Search start");
                search_length = prm.get_double("Search length");
                search_max_step = prm.get_integer("Search max step");
                search_tolerance = prm.get_double("Search tolerance");
                search_distance = prm.get_double("Search distance");
                positive_dipping = prm.get_bool("Positive dipping");
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

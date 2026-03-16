/*
  Copyright (C) 2017 - 2023 by the authors of the ASPECT code.

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
#include <aspect/initial_composition/world_builder_feature_surface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    WBFeatureSurface<dim>::initialize()
    {
      AssertThrow (this->introspection().compositional_name_exists("spcrust"),
                   ExcMessage("The 'world builder feature surface' initial composition requires the existence of a compositional field "
                              "named 'spcrust'. This field does not exist."));

      spcrust_index = this->introspection().compositional_index_for_name("spcrust");
    }

    template <int dim>
    double
    WBFeatureSurface<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const

    {
      if (compositional_index == spcrust_index)
        {
          const double depth = this->get_geometry_model().depth(position);

          std::array<double, 3> point = {{position[0], position[1], 0.0}};
          auto plane_distances = this->get_world_builder().distance_to_plane(point, depth, "Slab");
          const double distance = plane_distances.get_distance_from_surface(); // the distance of the query point to the feature surface

          double minimum_distance, maximum_distance;
          minimum_distance = minimum_distance_surface;
          if (function_type==1)
            {
              // linear
              maximum_distance = this->compute_m_dist_at_depth_linear(depth, maximum_distance_surface, maximum_distance_deep, depth_deep);
            }
          else
            {
              // constant
              maximum_distance = maximum_distance_surface;
            }

          if (distance > minimum_distance && distance < maximum_distance)
            return 1.0;
          else
            return 0.0;
        }
      else
        {
          return 0.0;
        }
    }

    template <int dim>
    void
    WBFeatureSurface<dim>::declare_parameters (ParameterHandler &prm)
    {
      //todo_wb
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("World builder feature surface");
        {
          prm.declare_entry("Distribution function", std::string("constant"), Patterns::Anything(),
                            "What type of distribution function to use for the shear zone thickness."
                            "Options are constant, linear.");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    WBFeatureSurface<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("World builder feature surface");
        {
          const std::string function_type_string = prm.get("Distribution function");
          function_type = match_function_string_input_to_unsigned(function_type_string);


          minimum_distance_surface = 0.0;
          maximum_distance_surface = 15e3;  // m
          maximum_distance_deep = 7.5e3;  // m
          depth_deep = 100e3;

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    double
    WBFeatureSurface<dim>::compute_m_dist_at_depth_linear(double depth, double m_dist,
                                                          double m_dist_deep, double depth_m) const
    {
      if (depth <= 0.0)
        {
          return m_dist;
        }
      else if (depth >= depth_m)
        {
          return m_dist_deep;
        }
      else
        {
          return m_dist + (m_dist_deep - m_dist) * (depth / depth_m);
        }
    }

    template <int dim>
    unsigned
    WBFeatureSurface<dim>::match_function_string_input_to_unsigned(const std::string &input) const
    {
      static const std::unordered_map<std::string, unsigned> mapping =
      {
        {"constant", 0},
        {"linear", 1}
      };

      auto it = mapping.find(input);
      if (it != mapping.end())
        {
          return it->second;
        }
      else
        {
          throw std::invalid_argument("Invalid input: " + input);
        }
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(WBFeatureSurface,
                                              "world builder feature surface",
                                              "A class that implements initial conditions by calculating"
                                              "the distance to a feature surface in the world builder."
                                             )
  }
}

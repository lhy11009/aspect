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


#ifndef _aspect_initial_composition_world_builder_feature_surface_h
#define _aspect_initial_composition_world_builder_feature_surface_h

#include <aspect/initial_composition/interface.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>

namespace aspect
{
  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements a method to determine
     * the composition based on a distance to a worldbuilder feature
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class WBFeatureSurface: public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Initialize the plugin.
         */
        void initialize () override;

        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position,
                                    const unsigned int compositional_index) const override;

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
        //todo_wb
        unsigned function_type;

        double minimum_distance_surface;
        double maximum_distance_surface;
        double maximum_distance_deep;
        double depth_deep;

        unsigned spcrust_index;

        /**
         * Compute the distribution of thickness versus depth from a linear function
         */
        double compute_m_dist_at_depth_linear(double depth, double m_dist, double m_dist_deep, double depth_m) const;

        /**
         * Match the string input of distribution function to a unsigned int
         */
        unsigned match_function_string_input_to_unsigned(const std::string &input) const;
    };
  }
}



#endif

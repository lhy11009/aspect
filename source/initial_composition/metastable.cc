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


#include <aspect/initial_composition/metastable.h>
#include <aspect/initial_temperature/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>


namespace aspect
{
  namespace InitialComposition
  {
    template <int dim>
    void
    Metastable<dim>::initialize()
    {
      AssertThrow (this->introspection().compositional_name_exists("metastable"),
                   ExcMessage("The 'metastable table lookup' initial composition requires the existence of a compositional field "
                              "named 'metastable'. This field does not exist."));

      // Make sure we keep track of the initial temperature manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_temperature_manager = this->get_initial_temperature_manager_pointer();

      // Make sure we keep track of the initial composition manager and
      // that it continues to live beyond the time when the simulator
      // class releases its pointer to it.
      initial_composition_manager = this->get_initial_composition_manager_pointer();

      metastable_index = this->introspection().compositional_index_for_name("metastable");

      if (this->introspection().compositional_name_exists("meta_grain_size"))
        {
          has_meta_grain_size = true;
          meta_grain_size_index = this->introspection().compositional_index_for_name("meta_grain_size");
        }
      else
        {
          has_meta_grain_size = false;
          meta_grain_size_index = std::numeric_limits<unsigned int>::max();
        }
    }


    template <int dim>
    double
    Metastable<dim>::
    initial_composition (const Point<dim> &position,
                         const unsigned int compositional_index) const
    {

      const double depth = this->get_geometry_model().depth(position);

      double reference_density;
      if (this->get_adiabatic_conditions().is_initialized())
        {
          reference_density = this->get_adiabatic_conditions().density(position);
        }
      else
        {
          reference_density = 3300.0;
        }

      const double gravity_norm = this->get_gravity_model().gravity_vector(position).norm();

      const double pressure_depth_derivative = gravity_norm*reference_density;

      const double temperature = initial_temperature_manager->initial_temperature(position);

      const double  depth_deviation = depth - transition_depth - transition_slope / pressure_depth_derivative
                                      * (temperature - transition_temperature);
      if (compositional_index == metastable_index)
        {
          // todo_gz
          // we make this variable 1.0 if equilirbrium phase transition happens
          // to make a equilibrium mantle at the start.
          // A detail included here is allowing a transition width as well
          // to make the transition smooth. If this is not designated, assign
          // a 0.0 value despite whatever value is used in the material model.
          double metastable = (depth_deviation > -1.0 * transition_width)? 1.0: 0.0;
          return metastable;
        }
      else if (compositional_index == metastable_index + 5)
        {
          // also handle the is_saturation state, set to 1.0 in the case of equilibrium
          // transition
          double metastable_is = (depth_deviation > -1.0 * transition_width)? 1.0: 0.0;

          return metastable_is;
        }
      else if (has_meta_grain_size && compositional_index == meta_grain_size_index)
        {
          return initial_grain_size;
        }
      return 0.0;
    }

    template <int dim>
    void
    Metastable<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Metastable");
        {
          prm.declare_entry ("Phase transition depth", "410e3", Patterns::Double (),
                             "Depth of the metastable phase transition");

          prm.declare_entry ("Phase transition width", "13e3", Patterns::Double (),
                             "Width of the metastable phase transition. This is only used to "
                             "facilitate an initial equilibrium mantle. The value doesn't need to "
                             "agree with the value in the material model");

          prm.declare_entry ("Phase transition temperature", "1780.0", Patterns::Double (),
                             "Temperature of the metastable phase transition");

          prm.declare_entry ("Phase transition Clapeyron slope", "2e6", Patterns::Double (),
                             "Claperon slope of the metastable phase transition");

          prm.declare_entry ("Initial grain size", "1e-2", Patterns::Double (),
                             "Initial grain size, used to track grain size evolution");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Metastable<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Metastable");
        {
          transition_depth = Utilities::string_to_double(prm.get("Phase transition depth"));
          transition_width = Utilities::string_to_double(prm.get("Phase transition width"));
          transition_temperature = Utilities::string_to_double(prm.get("Phase transition temperature"));
          transition_slope = Utilities::string_to_double(prm.get("Phase transition Clapeyron slope"));
          initial_grain_size = Utilities::string_to_double(prm.get("Initial grain size"));
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Metastable,
                                              "metastable",
                                              "A class that implements initial conditions for the metastable field "
                                              "by converting the initial temperature field through a look up table. "
                                              "Note that this plugin only works if there is a compositional field "
                                              "called `metastable', and an additional look up table that can convert "
                                              "pressure and temperature to metastable. "
                                              "For all compositional fields except metastable this plugin returns 0.0, "
                                              "and they are therefore not changed as long as the default `add' "
                                              "operator is selected for this plugin.")
  }
}

/*
  Copyright (C) 2011 - 2022 by the authors of the ASPECT code.

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

#include <aspect/global.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/gravity_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>

#include "adiabat_table_lookup.h"

namespace aspect
{
  namespace InitialTemperature
  {
    template <int dim>
    AdiabatTableLookup<dim>::AdiabatTableLookup ()
      = default;

    template <int dim>
    void
    AdiabatTableLookup<dim>::
    initialize()
    {
      material_lookup = std::make_unique<Utilities::StructuredDataLookup<2>>(7,1.0);
      material_lookup->load_file(data_directory+material_file_name,
                                 this->get_mpi_communicator());
    }

    template <int dim>
    double
    AdiabatTableLookup<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      const double pressure = this->get_adiabatic_conditions().pressure(position);

      double temperature;
      if (pressure_first)
        temperature = material_lookup->get_data({pressure/1e5, surface_entropy}, 0);
      else
        temperature = material_lookup->get_data({surface_entropy, pressure/1e5}, 0);

      return temperature;
    }

    template <int dim>
    void
    AdiabatTableLookup<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection("Adiabat table lookup");
        {
          prm.declare_entry ("Surface entropy", "0",
                             Patterns::Double(),
                             "The surface entropy for the profile.");
          prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/opxtable/",
                             Patterns::DirectoryName (),
                             "The path to the model data. The path may also include the special "
                             "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                             "in which the ASPECT source files were located when ASPECT was "
                             "compiled. This interpretation allows, for example, to reference "
                             "files located in the `data/' subdirectory of ASPECT.");
          prm.declare_entry ("Material file name", "material_table.txt",
                             Patterns::List (Patterns::Anything()),
                             "The file name of the material data.");
          prm.declare_entry ("Pressure first","false",
                             Patterns::Bool (),
                             "Whether the pressure is the first coordinate in the table.");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }

    template <int dim>
    void
    AdiabatTableLookup<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial temperature model");
      {
        prm.enter_subsection("Adiabat table lookup");
        {
          data_directory              = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
          material_file_name          = prm.get ("Material file name");
          pressure_first              = prm.get_bool ("Pressure first");
          surface_entropy             = prm.get_double ("Surface entropy");
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
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(AdiabatTableLookup,
                                              "adiabat table lookup",
                                              "Specify the initial temperature by looking up an adiabatic profile. "
                                              "from a look-up table. An entropy value is specified for the adiabatic profile.")
  }
}

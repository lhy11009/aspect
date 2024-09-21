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

#ifndef _aspect_initial_temperature_adiabat_table_lookup_h
#define _aspect_initial_temperature_adiabat_table_lookup_h

#include <aspect/global.h>


#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that implements temperature initial conditions based on a
     * look-up for the adiabatic condition in the P,S -> T lookup table
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class AdiabatTableLookup : public Interface<dim>, public SimulatorAccess<dim>
    {
      public:
        /**
         * Constructor.
         */
        AdiabatTableLookup ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program after parse_parameters is run and after
         * the SimulatorAccess (if applicable) is initialized.
         */
        void
        initialize () override;

        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

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
         * Starting entropy for the profile.
         */
        double surface_entropy;


        /**
         * Pointer to the StructuredDataLookup object that holds the material data.
         */
        std::unique_ptr<Utilities::StructuredDataLookup<2>> material_lookup;


        /**
         * Information about the location of data files.
         */
        std::string data_directory;
        std::string material_file_name;

        /*
          * whether pressure is the first indice in the table
        */
        bool pressure_first;

    };
  }
}

#endif

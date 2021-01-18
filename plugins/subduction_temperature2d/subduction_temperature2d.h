/*
  Copyright (C) 2016 - 2019 by the authors of the ASPECT code.

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


#ifndef _aspect_initial_temperature_anelasticity_temperature_h
#define _aspect_initial_temperature_anelasticity_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <cmath>
#include <algorithm>
#include <functional>

namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that describes an initial temperature field for a 2D or 3D shear wave velocity (Vs) model.
     * Vs values are converted to temperature using the anelasticity parameterization of Yamauchi & Takei (2016).
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Subduction2T : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        Subduction2T ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        // Avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return the boundary temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const;

        /**
         * consist model temperature with mantle phase transition
        */
        double
        below_lith_model (const double r, 
                          const double mantle_temperature, 
                          const double lith_depth,
                          const std::vector<double>& phases_below_lith_depth,
                          const std::vector<double>& phases_below_lith_temperature) const;
        
        /**
         * temperature from a slab temperature model 
         */ 
        double slab_model (const double r, const double phi) const;

        double katrina_mantle_temperature(const double r) const;

        /**
         * Calculate temperature from a half-space cooling model
         * We use the internal temperature and the surface temperature so that it is easier to plug in.
         */ 
        double 
        half_space_cooling (const double internal_temperature, 
                            const double surface_temperature,
                            const double depth, 
                            const double thermal_diffusivity,
                            const double age) const;
        /**
         * substract adiabatic temperature if needed
         */ 
        double
        temperature_substract_adiabatic(const double r, const double temperature) const;

        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:

        /**
         * Convert from cartesian to polar coordinates
         */
        Point<dim> cartesian_to_polar(const Point<dim> &cartesian_position) const;

        /**
         * todo
         * method of slab boudnary
         */
        // method_for_slab_boundary;


    };
  }
}

#endif

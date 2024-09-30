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


#ifndef _aspect_initial_composition_metastable_h
#define _aspect_initial_composition_metastable_h

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
     * A class that implements initial conditions for the metastable field
     * Note that this plugin only
     * works if there is a compositional field called 'metastable'.
     * All compositional fields except metastable are not changed by this plugin.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Metastable : public Interface<dim>,
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
        /**
         * Index of the metastable in the compositional fields
         */
        unsigned metastable_index;
        double transition_depth;
        double transition_width;
        double transition_temperature;
        double transition_slope;

        /**
         * A shared pointer to the initial temperature object
         * that ensures that the current object can continue
         * to access the initial temperature object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialTemperature::Manager<dim>> initial_temperature_manager;

        /**
         * A shared pointer to the initial composition object
         * that ensures that the current object can continue
         * to access the initial composition object beyond the
         * first time step.
         */
        std::shared_ptr<const aspect::InitialComposition::Manager<dim>> initial_composition_manager;
    };
  }
}


#endif

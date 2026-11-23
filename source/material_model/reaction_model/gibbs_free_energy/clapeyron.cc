/*
  Copyright (C) 2011 - 2026 by the authors of the ASPECT code.

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

#include <aspect/material_model/reaction_model/gibbs_free_energy/clapeyron.h>

#include <deal.II/base/patterns.h>

#include <cmath>

namespace aspect
{
  namespace MaterialModel
  {
    namespace ReactionModel
    {
      namespace GibbsFreeEnergyModel
      {
        template <int dim>
        double Clapeyron<dim>::compute_delta_gibbs_free_energy(const double temperature,
                                                               const double pressure,
                                                               const GibbsFreeEnergyInputs &,
                                                               const unsigned int reaction_index) const
        {
          AssertIndexRange(reaction_index, transitions.size());
          const auto &transition = transitions[reaction_index];
          const double pressure_offset = (pressure - transition.reference_pressure)
                                         - transition.slope * (temperature - transition.reference_temperature);
          return transition.delta_volume * pressure_offset;
        }



        template <int dim>
        void Clapeyron<dim>::declare_parameters(ParameterHandler &prm)
        {
          const std::string list_description = " Supply one '|'-separated entry per reaction assigned to this model.";
          const Patterns::List signed_values(Patterns::Double(), 1, Patterns::List::max_int_value, "|");
          const Patterns::List nonnegative_values(Patterns::Double(0), 1, Patterns::List::max_int_value, "|");

          prm.enter_subsection("Clapeyron Gibbs free energy");
          {
            prm.declare_entry("Clapeyron slopes", "0", signed_values,
                              "Slope dP/dT of the equilibrium phase boundary. Units: Pa/K." + list_description);
            prm.declare_entry("Reference pressures", "0", signed_values,
                              "Pressure on the equilibrium phase boundary at the corresponding reference temperature. Units: Pa."
                              + list_description);
            prm.declare_entry("Reference temperatures", "0", nonnegative_values,
                              "Temperature on the equilibrium phase boundary at the corresponding reference pressure. Units: K."
                              + list_description);
            prm.declare_entry("Molar volume changes", "0", signed_values,
                              "Signed, nonzero difference V_B-V_A between the molar volumes of phases A and B. Units: m^3/mol."
                              + list_description);
          }
          prm.leave_subsection();
        }



        template <int dim>
        void Clapeyron<dim>::parse_parameters(ParameterHandler &prm, const unsigned int n_reactions)
        {
          prm.enter_subsection("Clapeyron Gibbs free energy");
          {
            const auto read_values = [&prm, n_reactions](const std::string &name)
            {
              const auto values = Utilities::string_to_double(Utilities::split_string_list(prm.get(name), '|'));
              AssertThrow(values.size() == n_reactions,
                          ExcMessage("The 'Clapeyron Gibbs free energy/" + name + "' parameter must have exactly "
                                     + std::to_string(n_reactions) + " '|'-separated entries."));
              for (const double value : values)
                AssertThrow(std::isfinite(value),
                            ExcMessage("The 'Clapeyron Gibbs free energy/" + name + "' parameter must contain finite values."));
              return values;
            };

            const auto slopes = read_values("Clapeyron slopes");
            const auto pressures = read_values("Reference pressures");
            const auto temperatures = read_values("Reference temperatures");
            const auto volumes = read_values("Molar volume changes");

            transitions.resize(n_reactions);
            for (unsigned int r = 0; r < n_reactions; ++r)
              {
                AssertThrow(volumes[r] != 0.0,
                            ExcMessage("The 'Clapeyron Gibbs free energy/Molar volume changes' parameter must contain nonzero values."));
                transitions[r] = {slopes[r], pressures[r], temperatures[r], volumes[r]};
              }
          }
          prm.leave_subsection();
        }

        ASPECT_REGISTER_GIBBS_FREE_ENERGY_MODEL(Clapeyron,
                                                "Clapeyron Gibbs free energy",
                                                "Computes G_B-G_A as the molar volume difference times the signed pressure "
                                                "offset from a linearized equilibrium boundary with a specified Clapeyron slope.")
      }
    }
  }
}

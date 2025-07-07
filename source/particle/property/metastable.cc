/*
  Copyright (C) 2022 by the authors of the ASPECT code.

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

#include <aspect/particle/property/metastable.h>
// #include <aspect/material_model/metastable.h>
#include <aspect/initial_composition/interface.h>

namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      template <int dim>
      Metastable<dim>::Metastable ()
        :
        material_inputs(1,0),
        material_outputs(1,0)
      {}



      template <int dim>
      void
      Metastable<dim>::initialize ()
      {
        material_inputs  = MaterialModel::MaterialModelInputs<dim>(1, this->n_compositional_fields());
        material_outputs = MaterialModel::MaterialModelOutputs<dim>(1, this->n_compositional_fields());

        AssertThrow(this->introspection().compositional_name_exists("metastable"),
                    ExcMessage("This particle property only makes sense if "
                               "there is a compositional field named 'metastable'."));

        metastable_index = this->introspection().compositional_index_for_name("metastable");



        if (this->introspection().compositional_name_exists("meta_x0"))
          {
            AssertThrow(this->introspection().compositional_name_exists("meta_x1") &&
                        this->introspection().compositional_name_exists("meta_x2") &&
                        this->introspection().compositional_name_exists("meta_x3") &&
                        this->introspection().compositional_name_exists("meta_is") &&
                        this->introspection().compositional_name_exists("meta_rate"),
                        ExcMessage("This particle property only makes sense if "
                                   "the compositional fields named 'meta_x0', 'meta_x1',"
                                   "'meta_x2', 'meta_x3' and `meta_is` are included together."));

            AssertThrow(this->introspection().compositional_index_for_name("meta_x0") == metastable_index + 1 &&
                        this->introspection().compositional_index_for_name("meta_x1") == metastable_index + 2 &&
                        this->introspection().compositional_index_for_name("meta_x2") == metastable_index + 3 &&
                        this->introspection().compositional_index_for_name("meta_x3") == metastable_index + 4 &&
                        this->introspection().compositional_index_for_name("meta_is") == metastable_index + 5 &&
                        this->introspection().compositional_index_for_name("meta_rate") == metastable_index + 6,
                        ExcMessage("This particle property only makes sense if "
                                   "the compositional fields named 'metastable', 'meta_x0', 'meta_x1',"
                                   "'meta_x2', 'meta_x3' and 'meta_is' are consequtive."));

            solve_fourth_order_kinetics = true;

          }
        else
          solve_fourth_order_kinetics = false;

        // saved for later: include different kinetics
        // const unsigned int metastable_gradient_index_diff = solve_fourth_order_kinetics? 6: 1;
        // assert(this->introspection().compositional_index_for_name("meta_rate") == metastable_index + metastable_gradient_index_diff);
      }

      template <int dim>
      void
      Metastable<dim>::initialize_one_particle_property(const Point<dim> &position,
                                                        std::vector<double> &data) const
      {
        // Set the initial composition to the initial metastable.
        data.push_back(this->get_initial_composition_manager().initial_composition(position,metastable_index));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_x0")));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_x1")));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_x2")));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_x3")));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_is")));
        data.push_back(this->get_initial_composition_manager().initial_composition(position,this->introspection().compositional_index_for_name("meta_rate")));
        if (this->introspection().compositional_name_exists("meta_grain_size"))
          {
            const unsigned int meta_grain_size_index = this->introspection().compositional_index_for_name("meta_grain_size");
            data.push_back(this->get_initial_composition_manager().initial_composition(position,meta_grain_size_index));
          }
      }



      template <int dim>
      void
      Metastable<dim>::update_particle_property(const unsigned int data_position,
                                                const Vector<double> &solution,
                                                const std::vector<Tensor<1,dim>> &gradients,
                                                typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        // Prepare material inputs
        material_inputs.compute_metastable_reaction = true;

        material_inputs.position[0] = particle->get_location();

        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                      &(this->get_dof_handler()));

        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];

        for (unsigned int i = 0; i <= 5; ++i)
          material_inputs.composition[0][metastable_index+i] = particle->get_properties()[data_position+i];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        material_inputs.requested_properties = MaterialModel::MaterialProperties::reaction_terms;

        // Decide whether this particle is beyond the equilibrium state
        const double eq_P = this->get_material_model().computeEqP(material_inputs.temperature[0]);
        const bool beyond_eq = (material_inputs.pressure[0] > eq_P);

        // Evaluate the material model and get the kinetic change in metastable status
        if (beyond_eq)
          this->get_material_model().evaluate(material_inputs,
                                              material_outputs);

        // Assign values to the particle properties
        // Note for particles below the equilibrium state, all entries are assigned 0.0
        const bool debug_reaction_output = false;
        if (solve_fourth_order_kinetics)
          {
            double metastable_old = 0.0;

            // debug outputs
            if (debug_reaction_output && material_inputs.pressure[0] > 1e10)
              {
                std::cout << "T = " << material_inputs.temperature[0] << std::endl;
                std::cout << "P = " << material_inputs.pressure[0] << std::endl;
              }

            for (unsigned int i = 0; i <= 6; ++i)
              {
                if (beyond_eq)
                  {
                    // todo_gz
                    // record the old value
                    if (i==0)
                      metastable_old = particle->get_properties()[data_position];

                    // update value in case it's beyond the quilibrium
                    // take a threshold on the 3rd component.
                    // because the expresion is (1- exp(-meta_x3))
                    particle->get_properties()[data_position+i] += material_outputs.reaction_terms[0][metastable_index+i];
                  }
                else
                  particle->get_properties()[data_position+i] = 0.0;
                // debug information
                if (debug_reaction_output && i==0 && material_inputs.temperature[0] > 917 &&  material_inputs.temperature[0] < 918.0 && this->get_timestep_number() > 0)
                  {
                    std::cout << "Find the point" << std::endl;
                    std::cout << "T = " << material_inputs.temperature[0] << std::endl;
                    std::cout << "P = " << material_inputs.pressure[0] << std::endl;
                  }
                if (debug_reaction_output && material_inputs.pressure[0] > 1e10)
                  {
                    std::cout << "i = " << i << ", reaction_terms = " << material_outputs.reaction_terms[0][metastable_index+i] << std::endl;
                    if (i == 6)
                      std::cout<<"\n" << std::endl;
                  }
              }
            if (this->introspection().compositional_name_exists("meta_grain_size"))
              {
                const unsigned int meta_grain_size_index = this->introspection().compositional_index_for_name("meta_grain_size");

                const double metastable_threshold_for_grain_growth = 0.7;
                const double grain_density_threshold_for_grain_growth = 1000;
                const double metastable_threshold_singular_min = 0.01;
                const double metastable_threshold_singular_max = 0.99;

                const double metastable = particle->get_properties()[data_position];
                const double grain_density = particle->get_properties()[data_position+1];

                // todo
                // the first option check at least 70 % of grains are converted and grain density is large
                // the second (combined) condition makes sure that kinetics doesn't finish in one single step
                // the thirs condition makes sure that kinetics haven't finished in the last step
                if (metastable > metastable_threshold_for_grain_growth && grain_density > grain_density_threshold_for_grain_growth &&
                    !(metastable_old < metastable_threshold_singular_min && metastable > metastable_threshold_singular_max)&&
                    !(metastable_old > metastable_threshold_singular_max))
                  {
                    const double initial_grain_size_from_kinetics = std::pow(6.0 * (metastable/grain_density) / M_PI, 1/3.0);
                    particle->get_properties()[data_position+meta_grain_size_index-metastable_index] = compute_grain_size(material_inputs.pressure[0], material_inputs.temperature[0], initial_grain_size_from_kinetics);

                    // todo_gz
                    // std::cout << "Assigning grain size:" << std::endl;
                    // std::cout << "grain_size = " << particle->get_properties()[data_position+meta_grain_size_index-metastable_index] << std::endl;
                  }
              }
          }

        if (beyond_eq)
          particle->get_properties()[data_position+6] = material_outputs.reaction_terms[0][metastable_index+6];
        else
          particle->get_properties()[data_position+6] = 0.0;

        // Reset negative values to 0.0 value. Because we use the cell to compute the reaction terms, it's possible the updated value of a single particle is negative
        for (unsigned int i = 0; i <= 5; ++i)
          {
            if (particle->get_properties()[data_position+i] < 0.0)
              particle->get_properties()[data_position+i] = 0.0;
          }
      }

      template <int dim>
      double Metastable<dim>::compute_grain_size(const double pressure, const double temperature, const double initial_grain_size) const
      {
        // estimated timescale for subduction and grain growth
        const double grain_growth_timescale = 1e6 * year_in_seconds;

        // growth kinetics from Wd, from Dannberg 2017 paper supplementary
        const double grain_growth_activation_energy = 6.62e5; // Eg (j/mol)
        const double grain_growth_activation_volume = 0;
        const double m = 3;  // pg
        const double grain_growth_rate_constant = 3.02e-4; // k0 (m^pg / s)

        double A = std::exp(-(grain_growth_activation_energy + pressure * grain_growth_activation_volume) / (constants::gas_constant * temperature));
        double inside = std::pow(initial_grain_size, m) + grain_growth_rate_constant * A * grain_growth_timescale;

        double grain_size = std::pow(inside, 1.0 / m);

        return grain_size;
      }


      template <int dim>
      InitializationModeForLateParticles
      Metastable<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
      }

      template <int dim>
      AdvectionField
      Metastable<dim>::advection_field_for_boundary_initialization (const unsigned int property_component) const
      {
        if (this->introspection().compositional_name_exists("meta_grain_size"))
          {
            Assert(property_component <= 7,
                   ExcMessage("The metastable particle property "
                              "only has eight components"));
          }
        else
          {
            Assert(property_component <= 6,
                   ExcMessage("The metastable particle property "
                              "only has seven components"));
          }
        return AdvectionField::composition(metastable_index+property_component);
      }


      template <int dim>
      UpdateTimeFlags
      Metastable<dim>::need_update() const
      {
        return update_time_step;
      }



      template <int dim>
      UpdateFlags
      Metastable<dim>::get_needed_update_flags () const
      {
        return update_values | update_gradients;
      }



      template <int dim>
      std::vector<std::pair<std::string, unsigned int>>
      Metastable<dim>::get_property_information() const
      {
        if (this->introspection().compositional_name_exists("meta_grain_size"))
          return {{"kinetic metastable",1}, {"kinetic meta_x0", 1}, {"kinetic meta_x1", 1}, {"kinetic meta_x2", 1}, {"kinetic meta_x3", 1}, {"kinetic meta_is", 1}, {"kinetic meta_rate", 1}, {"kinetic meta_grain_size", 1}};
        else
          return {{"kinetic metastable",1}, {"kinetic meta_x0", 1}, {"kinetic meta_x1", 1}, {"kinetic meta_x2", 1}, {"kinetic meta_x3", 1}, {"kinetic meta_is", 1}, {"kinetic meta_rate", 1}};
      }
    }
  }
}

// explicit instantiations
namespace aspect
{
  namespace Particle
  {
    namespace Property
    {
      ASPECT_REGISTER_PARTICLE_PROPERTY(Metastable,
                                        "metastable",
                                        "A plugin in which the particle property is "
                                        "defined as the evolving metastable of a particle. "
                                        "See the metastable material model "
                                        "documentation for more detailed information.")

    }
  }
}

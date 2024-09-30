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
                        this->introspection().compositional_name_exists("meta_is"),
                        ExcMessage("This particle property only makes sense if "
                                   "the compositional fields named 'meta_x0', 'meta_x1',"
                                   "'meta_x2', 'meta_x3' and `meta_is` are included together."));

            AssertThrow(this->introspection().compositional_index_for_name("meta_x0") == metastable_index + 1 &&
                        this->introspection().compositional_index_for_name("meta_x1") == metastable_index + 2 &&
                        this->introspection().compositional_index_for_name("meta_x2") == metastable_index + 3 &&
                        this->introspection().compositional_index_for_name("meta_x3") == metastable_index + 4 &&
                        this->introspection().compositional_index_for_name("meta_is") == metastable_index + 5,
                        ExcMessage("This particle property only makes sense if "
                                   "the compositional fields named 'metastable', 'meta_x0', 'meta_x1',"
                                   "'meta_x2', 'meta_x3' and 'meta_is' are consequtive."));

            solve_fourth_order_kinetics = true;

          }
        else
          solve_fourth_order_kinetics = false;

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
      }



      template <int dim>
      void
      Metastable<dim>::update_particle_property(const unsigned int data_position,
                                                const Vector<double> &solution,
                                                const std::vector<Tensor<1,dim>> &gradients,
                                                typename ParticleHandler<dim>::particle_iterator &particle) const
      {
        // todo_metastable
        material_inputs.position[0] = particle->get_location();

        material_inputs.current_cell = typename DoFHandler<dim>::active_cell_iterator(*particle->get_surrounding_cell(),
                                                                                      &(this->get_dof_handler()));

        material_inputs.temperature[0] = solution[this->introspection().component_indices.temperature];

        material_inputs.pressure[0] = solution[this->introspection().component_indices.pressure];

        for (unsigned int d = 0; d < dim; ++d)
          material_inputs.velocity[0][d] = solution[this->introspection().component_indices.velocities[d]];

        for (unsigned int n = 0; n < this->n_compositional_fields(); ++n)
          material_inputs.composition[0][n] = solution[this->introspection().component_indices.compositional_fields[n]];

        material_inputs.composition[0][metastable_index] = particle->get_properties()[data_position];

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = gradients[d];
        material_inputs.strain_rate[0] = symmetrize (grad_u);

        material_inputs.requested_properties = MaterialModel::MaterialProperties::reaction_terms;

        this->get_material_model().evaluate(material_inputs,
                                            material_outputs);


        const bool debug_reaction_output = true;

        if (solve_fourth_order_kinetics)
          {
            for (unsigned int i = 0; i <= 5; ++i)
            {
              particle->get_properties()[data_position+i] += material_outputs.reaction_terms[0][metastable_index+i];
              // debug information
              if (debug_reaction_output && material_inputs.pressure[0] > 1e10){
                std::cout << "i = " << i << ", reaction_terms = " << material_outputs.reaction_terms[0][metastable_index+i] << std::endl;
              }
              if (material_inputs.temperature[0] > 917 &&  material_inputs.temperature[0] < 918.0 && this->get_timestep_number() > 0){
                std::cout << "Find the point" << std::endl;
              }
            }
          }
      }



      template <int dim>
      InitializationModeForLateParticles
      Metastable<dim>::late_initialization_mode () const
      {
        return interpolate_respect_boundary;
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
        return {{"metastable",1}, {"meta_x0", 1}, {"meta_x1", 1}, {"meta_x2", 1}, {"meta_x3", 1}, {"meta_is", 1}};
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


#include "slab2d_statistics.h"
#include <aspect/boundary_temperature/interface.h>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>


namespace aspect
{
  namespace Postprocess
  {
    template <int dim>
    std::pair<std::string,std::string>
    Slab2dStatistics<dim>::execute (TableHandler &statistics)
    {
      // create a quadrature formula based on the temperature element alone.
      const QGauss<dim> quadrature_formula (this->get_fe().base_element(this->introspection().base_elements.temperature).degree+1);
      const unsigned int n_q_points = quadrature_formula.size();
      const QGauss<dim> quadrature_formula_composition (this->get_fe().base_element(this->introspection().base_elements.compositional_fields).degree+1);
      const unsigned int n_q_points_composition = quadrature_formula_composition.size();
      AssertThrow(n_q_points == n_q_points_composition,
                  ExcMessage("In plugin slab_statistics, I take for granted that the quadrature formula for temperature and composition are the same")
                 );

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

//      std::vector<double> temperature_values(n_q_points);
      std::vector<double> compositional_values(n_q_points);

      double local_slab_volume_integral = 0;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
              {
                if (this->introspection().name_for_compositional_index(c) == "spcrust" ||
                    this->introspection().name_for_compositional_index(c) == "spharz"){
                  fe_values[this->introspection().extractors.compositional_fields[c]].get_function_values (this->get_solution(),
                     compositional_values);
                  for (unsigned int q=0; q<n_q_points; ++q){
                    if (this->get_geometry_model().depth(fe_values.quadrature_point(q)) > 40.0e3)
                      local_slab_volume_integral += compositional_values[q]*fe_values.JxW(q);
                  }
                }
              }
          }
      
      // compute the sum over all processors
      const double global_slab_volume_integral
        = Utilities::MPI::sum (local_slab_volume_integral, this->get_mpi_communicator());

      statistics.add_value ("Slab volume",
                            global_slab_volume_integral);

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Slab volume"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }

      }

      std::ostringstream output;
      output.precision(4);
      output << global_slab_volume_integral;

      return std::pair<std::string, std::string> ("Slab Statistics volume:",
                                                  output.str());
    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace Postprocess
  {
    ASPECT_REGISTER_POSTPROCESSOR(Slab2dStatistics,
                                  "slab2d statistics",
                                  "A postprocessor that computes some statistics about "
                                  "the 2d slab field.")
  }
}

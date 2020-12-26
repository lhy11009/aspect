
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

      FEValues<dim> fe_values (this->get_mapping(),
                               this->get_fe(),
                               quadrature_formula,
                               update_values   |
                               update_quadrature_points |
                               update_JxW_values);

      std::vector<double> temperature_values(n_q_points);

      double local_temperature_integral = 0;

      // compute the integral quantities by quadrature
      for (const auto &cell : this->get_dof_handler().active_cell_iterators())
        if (cell->is_locally_owned())
          {
            fe_values.reinit (cell);
            fe_values[this->introspection().extractors.temperature].get_function_values (this->get_solution(),
                                                                                         temperature_values);
            for (unsigned int q=0; q<n_q_points; ++q)
              {
                local_temperature_integral += temperature_values[q]*fe_values.JxW(q);
              }
          }

      // compute min/max by simply
      // looping over the elements of the
      // solution vector. the reason is
      // that minimum and maximum are
      // usually attained at the
      // boundary, and so taking their
      // values at Gauss quadrature
      // points gives an inaccurate
      // picture of their true values
      double local_min_temperature = std::numeric_limits<double>::max();
      double local_max_temperature = -std::numeric_limits<double>::max();
      const unsigned int temperature_block = this->introspection().block_indices.temperature;
      IndexSet range = this->get_solution().block(temperature_block).locally_owned_elements();
      for (unsigned int i=0; i<range.n_elements(); ++i)
        {
          const unsigned int idx = range.nth_index_in_set(i);
          const double val =  this->get_solution().block(temperature_block)(idx);

          local_min_temperature = std::min<double> (local_min_temperature, val);
          local_max_temperature = std::max<double> (local_max_temperature, val);
        }

      const double global_temperature_integral
        = Utilities::MPI::sum (local_temperature_integral, this->get_mpi_communicator());
      double global_min_temperature = 0;
      double global_max_temperature = 0;

      // now do the reductions that are
      // min/max operations. do them in
      // one communication by multiplying
      // one value by -1
      {
        double local_values[2] = { -local_min_temperature, local_max_temperature };
        double global_values[2];

        Utilities::MPI::max (local_values, this->get_mpi_communicator(), global_values);

        global_min_temperature = -global_values[0];
        global_max_temperature = global_values[1];
      }

      double global_mean_temperature = global_temperature_integral / this->get_volume();
      statistics.add_value ("Minimal temperature (K)",
                            global_min_temperature);
      statistics.add_value ("Average temperature (K)",
                            global_mean_temperature);
      statistics.add_value ("Maximal temperature (K)",
                            global_max_temperature);
      if ((this->get_fixed_temperature_boundary_indicators().size() > 0)
          &&
          (this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
           !=
           this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators())))
        statistics.add_value ("Average nondimensional temperature (K)",
                              (global_mean_temperature - global_min_temperature) /
                              (this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
                               -
                               this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators())));

      // also make sure that the other columns filled by the this object
      // all show up with sufficient accuracy and in scientific notation
      {
        const char *columns[] = { "Minimal temperature (K)",
                                  "Average temperature (K)",
                                  "Maximal temperature (K)"
                                };
        for (unsigned int i=0; i<sizeof(columns)/sizeof(columns[0]); ++i)
          {
            statistics.set_precision (columns[i], 8);
            statistics.set_scientific (columns[i], true);
          }

        if ((this->get_fixed_temperature_boundary_indicators().size() > 0)
            &&
            (this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators())
             !=
             this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators())))
          {
            statistics.set_precision ("Average nondimensional temperature (K)", 8);
            statistics.set_scientific ("Average nondimensional temperature (K)", true);
          }
      }

      std::ostringstream output;
      output.precision(4);
      output << global_min_temperature << " K, "
             << global_temperature_integral / this->get_volume() << " K, "
             << global_max_temperature << " K";

      return std::pair<std::string, std::string> ("Temperature min/avg/max:",
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

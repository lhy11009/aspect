/*
  Copyright (C) 2011 - 2024 by the authors of the ASPECT code.

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
#include <aspect/simulator_access.h>
#include <aspect/structured_data.h>

#include <aspect/geometry_model/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/gravity_model/interface.h>

#include <aspect/material_model/interface.h>
#include <aspect/material_model/equation_of_state/interface.h>
#include <aspect/material_model/utilities.h>
#include <aspect/utilities.h>

#include <deal.II/base/exceptions.h>
#include <boost/lexical_cast.hpp>
#include <list>


namespace aspect
{
  namespace MaterialModel
  {
    namespace MaterialUtilities
    {
      namespace Lookup
      {
        double
        MaterialLookup::specific_heat(const double temperature,
                                      const double pressure) const
        {
          return value(temperature,pressure,specific_heat_values,interpolation);
        }

        double
        MaterialLookup::density(const double temperature,
                                const double pressure) const
        {
          return value(temperature,pressure,density_values,interpolation);
        }

        double
        MaterialLookup::thermal_expansivity(const double temperature,
                                            const double pressure) const
        {
          return value(temperature,pressure,thermal_expansivity_values,interpolation);
        }

        double
        MaterialLookup::seismic_Vp(const double temperature,
                                   const double pressure) const
        {
          return value(temperature,pressure,vp_values,false);
        }

        double
        MaterialLookup::seismic_Vs(const double temperature,
                                   const double pressure) const
        {
          return value(temperature,pressure,vs_values,false);
        }

        double
        MaterialLookup::enthalpy(const double temperature,
                                 const double pressure) const
        {
          return value(temperature,pressure,enthalpy_values,true);
        }

        double
        MaterialLookup::dHdT (const double temperature,
                              const double pressure) const
        {
          const double h = value(temperature,pressure,enthalpy_values,interpolation);
          const double dh = value(temperature+delta_temp,pressure,enthalpy_values,interpolation);
          return (dh - h) / delta_temp;
        }

        double
        MaterialLookup::dHdp (const double temperature,
                              const double pressure) const
        {
          const double h = value(temperature,pressure,enthalpy_values,interpolation);
          const double dh = value(temperature,pressure+delta_press,enthalpy_values,interpolation);
          return (dh - h) / delta_press;
        }

        std::array<std::pair<double, unsigned int>,2>
        MaterialLookup::enthalpy_derivatives(const std::vector<double> &temperatures,
                                             const std::vector<double> &pressures,
                                             const unsigned int n_substeps) const
        {
          Assert(temperatures.size() == pressures.size(),ExcInternalError());
          const unsigned int n_q_points = temperatures.size();
          unsigned int n_T(0), n_p(0);
          double dHdT(0.0), dHdp(0.0);

          for (unsigned int q=0; q<n_q_points; ++q)
            {
              for (unsigned int p=0; p<n_q_points; ++p)
                {
                  if (std::fabs(temperatures[q] - temperatures[p]) > 100.0 * std::numeric_limits<double>::epsilon() * std::fabs(temperatures[q] + std::numeric_limits<double>::epsilon()))
                    {
                      for (unsigned int substep = 0; substep < n_substeps; ++substep)
                        {
                          const double step_ratio = static_cast<double>(substep)/static_cast<double>(n_substeps);
                          const double step_ratio_next = static_cast<double>(substep+1)/static_cast<double>(n_substeps);

                          const double current_pressure = pressures[q]
                                                          + step_ratio
                                                          * (pressures[p]-pressures[q]);
                          const double T1_substep = temperatures[q]
                                                    + step_ratio
                                                    * (temperatures[p]-temperatures[q]);
                          const double T2_substep = temperatures[q]
                                                    + step_ratio_next
                                                    * (temperatures[p]-temperatures[q]);
                          const double enthalpy2 = enthalpy(T2_substep,current_pressure);
                          const double enthalpy1 = enthalpy(T1_substep,current_pressure);
                          dHdT += (enthalpy2-enthalpy1)/(T2_substep-T1_substep);
                          ++n_T;
                        }
                    }
                  if (std::fabs(pressures[q] - pressures[p]) > 100.0 * std::numeric_limits<double>::epsilon() * std::fabs(pressures[q] + std::numeric_limits<double>::epsilon()))
                    {
                      for (unsigned int substep = 0; substep < n_substeps; ++substep)
                        {
                          const double step_ratio = static_cast<double>(substep)/static_cast<double>(n_substeps);
                          const double step_ratio_next = static_cast<double>(substep+1)/static_cast<double>(n_substeps);
                          const double current_temperature = temperatures[q]
                                                             + step_ratio
                                                             * (temperatures[p]-temperatures[q]);
                          const double p1_substep = pressures[q]
                                                    + step_ratio
                                                    * (pressures[p]-pressures[q]);
                          const double p2_substep = pressures[q]
                                                    + step_ratio_next
                                                    * (pressures[p]-pressures[q]);
                          const double enthalpy2 = enthalpy(current_temperature,p2_substep);
                          const double enthalpy1 = enthalpy(current_temperature,p1_substep);
                          dHdp += (enthalpy2-enthalpy1)/(p2_substep-p1_substep);
                          ++n_p;
                        }
                    }
                }
            }

          if ((n_T > 0)
              && (n_p > 0))
            {
              dHdT /= n_T;
              dHdp /= n_p;
            }

          std::array<std::pair<double, unsigned int>,2> derivatives;
          derivatives[0] = std::make_pair(dHdT,n_T);
          derivatives[1] = std::make_pair(dHdp,n_p);
          return derivatives;
        }

        double
        MaterialLookup::dRhodp (const double temperature,
                                const double pressure) const
        {
          const double rho = value(temperature,pressure,density_values,interpolation);
          const double drho = value(temperature,pressure+delta_press,density_values,interpolation);
          return (drho - rho) / delta_press;
        }

        unsigned int
        MaterialLookup::dominant_phase (const double temperature,
                                        const double pressure) const
        {
          if (!has_dominant_phase_column)
            AssertThrow(false, ExcMessage("You can not ask for the column with the dominant phase if it does not exist in the data file."));
          return value(temperature, pressure, dominant_phase_indices);
        }

        bool
        MaterialLookup::has_dominant_phase() const
        {
          return has_dominant_phase_column;
        }

        std::vector<std::string>
        MaterialLookup::phase_volume_column_names() const
        {
          return phase_column_names;
        }

        double
        MaterialLookup::phase_volume_fraction(const int phase_id,
                                              const double temperature,
                                              const double pressure) const
        {
          return value(temperature,pressure,phase_volume_fractions[phase_id],interpolation);
        }


        const std::vector<std::string> &
        MaterialLookup::get_dominant_phase_names() const
        {
          return dominant_phase_names;
        }

        double
        MaterialLookup::value (const double temperature,
                               const double pressure,
                               const Table<2, double> &values,
                               const bool interpol) const
        {
          const double nT = get_nT(temperature);
          const unsigned int inT = static_cast<unsigned int>(nT);

          const double np = get_np(pressure);
          const unsigned int inp = static_cast<unsigned int>(np);

          Assert(inT<values.n_rows(), ExcMessage("Attempting to look up a temperature value with index greater than the number of rows."));
          Assert(inp<values.n_cols(), ExcMessage("Attempting to look up a pressure value with index greater than the number of columns."));

          if (!interpol)
            return values[inT][inp];
          else
            {
              // compute the coordinates of this point in the
              // reference cell between the data points
              const double xi = nT-inT;
              const double eta = np-inp;

              Assert ((0 <= xi) && (xi <= 1), ExcInternalError());
              Assert ((0 <= eta) && (eta <= 1), ExcInternalError());

              // use these coordinates for a bilinear interpolation
              return ((1-xi)*(1-eta)*values[inT][inp] +
                      xi    *(1-eta)*values[inT+1][inp] +
                      (1-xi)*eta    *values[inT][inp+1] +
                      xi    *eta    *values[inT+1][inp+1]);
            }
        }

        unsigned int
        MaterialLookup::value (const double temperature,
                               const double pressure,
                               const Table<2, unsigned int> &values) const
        {
          const double nT = get_nT(temperature);
          const unsigned int inT = static_cast<unsigned int>(nT);

          const double np = get_np(pressure);
          const unsigned int inp = static_cast<unsigned int>(np);

          Assert(inT<values.n_rows(), ExcMessage("Attempting to look up a temperature value with index greater than the number of rows."));
          Assert(inp<values.n_cols(), ExcMessage("Attempting to look up a pressure value with index greater than the number of columns."));

          return values[inT][inp];
        }

        std::array<double,2>
        MaterialLookup::get_pT_steps() const
        {
          std::array<double,2> pt_steps;
          pt_steps[0] = delta_press;
          pt_steps[1] = delta_temp;
          return pt_steps;
        }

        double
        MaterialLookup::get_nT(const double temperature) const
        {
          double bounded_temperature=std::max(min_temp, temperature);
          bounded_temperature=std::min(bounded_temperature, max_temp-delta_temp);

          return (bounded_temperature-min_temp)/delta_temp;
        }

        double
        MaterialLookup::get_np(const double pressure) const
        {
          double bounded_pressure=std::max(min_press, pressure);
          bounded_pressure=std::min(bounded_pressure, max_press-delta_press);

          return (bounded_pressure-min_press)/delta_press;
        }

        HeFESToReader::HeFESToReader(const std::string &material_filename,
                                     const std::string &derivatives_filename,
                                     const bool interpol,
                                     const MPI_Comm comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=std::numeric_limits<double>::lowest();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=std::numeric_limits<double>::lowest();
          n_temperature=0;
          n_pressure=0;

          std::string temp;

          // Read material data
          {
            // Read data from disk and distribute among processes
            std::istringstream in(Utilities::read_and_distribute_file_content(material_filename, comm));

            bool parsed_first_column = false;
            unsigned int i = 0;
            double current_pressure = 0.0;
            double old_pressure = -1.0;
            while (!in.eof())
              {
                in >> current_pressure;
                if (in.fail())
                  {
                    in.clear();
                  }

                if (!parsed_first_column)
                  {
                    if (current_pressure > old_pressure)
                      old_pressure = current_pressure;
                    else if (current_pressure <= old_pressure)
                      {
                        n_pressure = i;
                        parsed_first_column = true;
                      }
                  }

                std::getline(in, temp);
                if (in.eof())
                  break;
                ++i;
              }

            in.clear();
            in.seekg (0, in.beg);

            n_temperature = i / n_pressure;

            Assert(i == n_temperature * n_pressure,
                   ExcMessage("Material table size not consistent."));

            density_values.reinit(n_temperature,n_pressure);
            thermal_expansivity_values.reinit(n_temperature,n_pressure);
            specific_heat_values.reinit(n_temperature,n_pressure);
            vp_values.reinit(n_temperature,n_pressure);
            vs_values.reinit(n_temperature,n_pressure);
            enthalpy_values.reinit(n_temperature,n_pressure);

            i = 0;
            while (!in.eof())
              {
                double P = 0.0;
                double depth,T;
                double rho,vb,vs,vp,vsq,vpq,h;
                std::string code;
                double alpha = 0.0;
                double cp = 0.0;

                in >> P >> depth >> T;
                if (in.fail())
                  in.clear();
                // conversion from [GPa] to [Pa]
                P *= 1e9;

                min_press=std::min(P,min_press);
                min_temp=std::min(T,min_temp);
                max_temp = std::max(T,max_temp);
                max_press = std::max(P,max_press);

                in >> rho;
                if (in.fail())
                  {
                    in.clear();
                    rho = density_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                else
                  rho *= 1e3; // conversion from [g/cm^3] to [kg/m^3]

                in >> vb;
                if (in.fail())
                  in.clear();

                in >> vs;
                if (in.fail())
                  {
                    in.clear();
                    vs = vs_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                in >> vp;
                if (in.fail())
                  {
                    in.clear();
                    vp = vp_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                in >> vsq >> vpq;

                in >> h;
                if (in.fail())
                  {
                    in.clear();
                    h = enthalpy_values[(i-1)%n_temperature][(i-1)/n_temperature];
                  }
                else
                  h *= 1e6; // conversion from [kJ/g] to [J/kg]

                std::getline(in, temp);
                if (in.eof())
                  break;

                density_values[i/n_pressure][i%n_pressure]=rho;
                thermal_expansivity_values[i/n_pressure][i%n_pressure]=alpha;
                specific_heat_values[i/n_pressure][i%n_pressure]=cp;
                vp_values[i/n_pressure][i%n_pressure]=vp;
                vs_values[i/n_pressure][i%n_pressure]=vs;
                enthalpy_values[i/n_pressure][i%n_pressure]=h;

                ++i;
              }

            delta_temp = (max_temp - min_temp) / (n_temperature - 1);
            delta_press = (max_press - min_press) / (n_pressure - 1);

            AssertThrow(max_temp >= 0.0, ExcMessage("Read in of Material header failed (max_temp)."));
            AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
            AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
            AssertThrow(max_press >= 0, ExcMessage("Read in of Material header failed (max_press)."));
            AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
            AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));
          }

          // If requested read derivative data
          if (derivatives_filename != "")
            {
              std::string temp;
              // Read data from disk and distribute among processes
              std::istringstream in(Utilities::read_and_distribute_file_content(derivatives_filename, comm));

              int i = 0;
              while (!in.eof())
                {
                  double P = 0.0;
                  double depth,T;
                  double cp,alpha,alpha_eff;
                  double temp1,temp2;

                  in >> P >> depth >> T;
                  if (in.fail())
                    in.clear();


                  in >> cp;
                  if (in.fail() || (cp <= std::numeric_limits<double>::min()))
                    {
                      in.clear();
                      cp = specific_heat_values[(i-1)%n_temperature][(i-1)/n_temperature];
                    }
                  else
                    cp *= 1e3; // conversion from [J/g/K] to [J/kg/K]

                  in >> alpha >> alpha_eff;
                  if (in.fail() || (alpha_eff <= std::numeric_limits<double>::min()))
                    {
                      in.clear();
                      alpha_eff = thermal_expansivity_values[(i-1)%n_temperature][(i-1)/n_temperature];
                    }
                  else
                    {
                      alpha *= 1e-5;
                      alpha_eff *= 1e-5;
                    }

                  in >> temp1 >> temp2;
                  if (in.fail())
                    in.clear();


                  std::getline(in, temp);
                  if (in.eof())
                    break;

                  specific_heat_values[i/n_pressure][i%n_pressure]=cp;
                  thermal_expansivity_values[i/n_pressure][i%n_pressure]=alpha_eff;

                  ++i;
                }
            }
        }

        PerplexReader::PerplexReader(const std::string &filename,
                                     const bool interpol,
                                     const MPI_Comm comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=std::numeric_limits<double>::lowest();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=std::numeric_limits<double>::lowest();
          n_temperature=0;
          n_pressure=0;
          has_dominant_phase_column = false;

          std::string temp;
          // Read data from disk and distribute among processes
          std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

          // The following lines read in a PerpleX tab file in standard format
          // The first 13 lines are a header in the format:
          // |<perplex version>
          // <table filename>
          // <grid dim>
          // <grid variable 1> (usually T(K) or P(bar))
          // <min grid variable 1>
          // <delta grid variable 1>
          // <n steps grid variable 1>
          // <grid variable 2> (usually T(K) or P(bar))
          // <min grid variable 2>
          // <delta grid variable 2>
          // <n steps grid variable 2>
          // Number of property columns in the table
          // Column names

          // First line is the Perplex version number
          std::getline(in, temp); // get next line, table file name

          std::getline(in, temp); // get next line, dimension of table
          unsigned int n_variables;
          in >> n_variables;
          AssertThrow (n_variables==2, ExcMessage("The PerpleX file " + filename + " must be two dimensional (P(bar)-T(K))."));

          std::getline(in, temp); // get next line, either T(K) or P(bar)

          for (unsigned int i=0; i<2; ++i)
            {
              std::string natural_variable;
              in >> natural_variable;

              if (natural_variable == "T(K)")
                {
                  std::getline(in, temp);
                  in >> min_temp;
                  std::getline(in, temp);
                  in >> delta_temp;
                  std::getline(in, temp);
                  in >> n_temperature;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else if (natural_variable == "P(bar)")
                {
                  std::getline(in, temp);
                  in >> min_press;
                  min_press *= 1e5;  // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> delta_press;
                  delta_press *= 1e5; // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> n_pressure;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else
                {
                  AssertThrow (false, ExcMessage("The start of the PerpleX file " + filename + " does not have the expected format."));
                }
            }

          in >> n_columns;
          std::getline(in, temp); // get next line, column labels

          // here we string match to assign properties to columns
          // column i in text file -> column j in properties
          // Properties are stored in the order rho, alpha, cp, vp, vs, h
          std::vector<int> prp_indices(6, -1);
          std::vector<int> phase_column_indices;
          unsigned int dominant_phase_column_index = numbers::invalid_unsigned_int;

          // First two columns should be P(bar) and T(K).
          // Here we find the order.
          std::string column_name;
          in >> column_name;

          std::string first_natural_variable;
          if (column_name == "P(bar)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "T(K)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be T(K)."));
            }
          else if (column_name == "T(K)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "P(bar)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be P(bar)."));
            }
          else
            {
              AssertThrow(false, ExcMessage("The first column name in the PerpleX lookup file " + filename + " should be P(bar) or T(K)."));
            }

          for (unsigned int n=2; n<n_columns; ++n)
            {
              in >> column_name;
              if (column_name == "rho,kg/m3")
                prp_indices[0] = n;
              else if (column_name == "alpha,1/K")
                prp_indices[1] = n;
              else if (column_name == "cp,J/K/kg")
                prp_indices[2] = n;
              else if (column_name == "vp,km/s")
                prp_indices[3] = n;
              else if (column_name == "vs,km/s")
                prp_indices[4] = n;
              else if (column_name == "h,J/kg")
                prp_indices[5] = n;
              else if (column_name == "phase")
                {
                  has_dominant_phase_column = true;
                  dominant_phase_column_index = n;
                }
              else if (column_name.length() > 3)
                {
                  if (column_name.substr(0,13).compare("vol_fraction_") == 0)
                    {
                      if (std::find(phase_column_names.begin(),
                                    phase_column_names.end(),
                                    column_name) != phase_column_names.end())
                        {
                          AssertThrow(false,
                                      ExcMessage("The PerpleX lookup file " + filename + " must have unique column names. "
                                                 "Sometimes, the same phase is stable with >1 composition at the same "
                                                 "pressure and temperature, so you may see several columns with the same name. "
                                                 "Either combine columns with the same name, or change the names."));
                        }
                      // Populate phase_column_names with the column name
                      // and phase_column_indices with the column index in the current lookup file.
                      phase_column_indices.push_back(n);
                      phase_column_names.push_back(column_name);
                    }
                }
            }
          AssertThrow(std::all_of(prp_indices.begin(), prp_indices.end(), [](int i)
          {
            return i>=0;
          }),
          ExcMessage("The PerpleX lookup file " + filename + " must contain columns with names "
                     "rho,kg/m3, alpha,1/K, cp,J/K/kg, vp,km/s, vs,km/s and h,J/kg."));

          std::getline(in, temp); // first data line

          AssertThrow(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
          AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
          AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
          AssertThrow(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
          AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
          AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));


          max_temp = min_temp + (n_temperature-1) * delta_temp;
          max_press = min_press + (n_pressure-1) * delta_press;

          density_values.reinit(n_temperature,n_pressure);
          thermal_expansivity_values.reinit(n_temperature,n_pressure);
          specific_heat_values.reinit(n_temperature,n_pressure);
          vp_values.reinit(n_temperature,n_pressure);
          vs_values.reinit(n_temperature,n_pressure);
          enthalpy_values.reinit(n_temperature,n_pressure);

          if (has_dominant_phase_column)
            dominant_phase_indices.reinit(n_temperature,n_pressure);

          phase_volume_fractions.resize(phase_column_names.size());
          for (auto &phase_volume_fraction : phase_volume_fractions)
            phase_volume_fraction.reinit(n_temperature,n_pressure);

          unsigned int i = 0;
          std::vector<double> previous_row_values(n_columns, 0.);

          while (!in.eof())
            {
              std::vector<double> row_values(n_columns);
              std::string phase;

              for (unsigned int n=0; n<n_columns; ++n)
                {
                  if (n == dominant_phase_column_index)
                    in >> phase;
                  else
                    in >> row_values[n]; // assigned as 0 if in.fail() == True

                  // P-T grids created with PerpleX-werami sometimes contain rows
                  // filled with NaNs at extreme P-T conditions where the thermodynamic
                  // models break down. These P-T regions are typically not relevant to
                  // geodynamic modeling (they most commonly appear above
                  // mantle liquidus temperatures at low pressures).
                  // More frustratingly, PerpleX-vertex occasionally fails to find a
                  // valid mineral assemblage in small, isolated regions within the domain,
                  // and so PerpleX-werami also returns NaNs for pixels within these regions.
                  // It is recommended that the user preprocesses their input
                  // files to replace these NaNs before plugging them into ASPECT.
                  // If this lookup encounters invalid doubles it replaces them
                  // with the most recent valid double.
                  if (in.fail())
                    {
                      in.clear();
                      row_values[n] = previous_row_values[n];
                    }
                }
              previous_row_values = row_values;

              std::getline(in, temp); // read next line
              if (in.eof())
                break;

              if (std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase) == dominant_phase_names.end())
                dominant_phase_names.push_back(phase);

              // The ordering of the first two columns in the PerpleX table files
              // dictates whether the inner loop is over temperature or pressure.
              // The first column is always the inner loop.
              // The following lines populate the material property tables
              // according to that implicit loop structure.
              if (first_natural_variable == "T(K)")
                {
                  density_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[1]];
                  specific_heat_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[2]];
                  vp_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[3]];
                  vs_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[4]];
                  enthalpy_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i%n_temperature][i/n_temperature] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i%n_temperature][i/n_temperature]=row_values[phase_column_indices[n]];
                    }
                }
              else // first_natural_variable == "P(bar)"
                {
                  density_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[1]];
                  specific_heat_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[2]];
                  vp_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[3]];
                  vs_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[4]];
                  enthalpy_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i/n_pressure][i%n_pressure] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i/n_pressure][i%n_pressure]=row_values[phase_column_indices[n]];
                    }
                }
              ++i;
            }
          AssertThrow(i == n_temperature*n_pressure, ExcMessage("Material table size not consistent with header."));

        }

        PerplexReaderMorb::PerplexReaderMorb(const std::string &filename,
                                             const bool interpol,
                                             const MPI_Comm comm)
        {
          /* Initializing variables */
          interpolation = interpol;
          delta_press=numbers::signaling_nan<double>();
          min_press=std::numeric_limits<double>::max();
          max_press=std::numeric_limits<double>::lowest();
          delta_temp=numbers::signaling_nan<double>();
          min_temp=std::numeric_limits<double>::max();
          max_temp=std::numeric_limits<double>::lowest();
          n_temperature=0;
          n_pressure=0;
          has_dominant_phase_column = false;

          std::string temp;
          // Read data from disk and distribute among processes
          std::istringstream in(Utilities::read_and_distribute_file_content(filename, comm));

          // The following lines read in a PerpleX tab file in standard format
          // The first 13 lines are a header in the format:
          // |<perplex version>
          // <table filename>
          // <grid dim>
          // <grid variable 1> (usually T(K) or P(bar))
          // <min grid variable 1>
          // <delta grid variable 1>
          // <n steps grid variable 1>
          // <grid variable 2> (usually T(K) or P(bar))
          // <min grid variable 2>
          // <delta grid variable 2>
          // <n steps grid variable 2>
          // Number of property columns in the table
          // Column names

          // First line is the Perplex version number
          std::getline(in, temp); // get next line, table file name

          std::getline(in, temp); // get next line, dimension of table
          unsigned int n_variables;
          in >> n_variables;
          AssertThrow (n_variables==2, ExcMessage("The PerpleX file " + filename + " must be two dimensional (P(bar)-T(K))."));

          std::getline(in, temp); // get next line, either T(K) or P(bar)

          for (unsigned int i=0; i<2; ++i)
            {
              std::string natural_variable;
              in >> natural_variable;

              if (natural_variable == "T(K)")
                {
                  std::getline(in, temp);
                  in >> min_temp;
                  std::getline(in, temp);
                  in >> delta_temp;
                  std::getline(in, temp);
                  in >> n_temperature;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else if (natural_variable == "P(bar)")
                {
                  std::getline(in, temp);
                  in >> min_press;
                  min_press *= 1e5;  // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> delta_press;
                  delta_press *= 1e5; // conversion from [bar] to [Pa]
                  std::getline(in, temp);
                  in >> n_pressure;
                  std::getline(in, temp); // get next line, either T(K), P(bar) or number of columns
                }
              else
                {
                  AssertThrow (false, ExcMessage("The start of the PerpleX file " + filename + " does not have the expected format."));
                }
            }

          in >> n_columns;
          std::getline(in, temp); // get next line, column labels

          // here we string match to assign properties to columns
          // column i in text file -> column j in properties
          // Properties are stored in the order rho, (no alpha, cp, vp, vs, h)
          std::vector<int> prp_indices(1, -1);
          std::vector<int> phase_column_indices;
          unsigned int dominant_phase_column_index = numbers::invalid_unsigned_int;

          // First two columns should be P(bar) and T(K).
          // Here we find the order.
          std::string column_name;
          in >> column_name;

          std::string first_natural_variable;
          if (column_name == "P(bar)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "T(K)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be T(K)."));
            }
          else if (column_name == "T(K)")
            {
              first_natural_variable = column_name;
              in >> column_name;
              AssertThrow(column_name == "P(bar)", ExcMessage("The second column name in PerpleX lookup file " + filename + " should be P(bar)."));
            }
          else
            {
              AssertThrow(false, ExcMessage("The first column name in the PerpleX lookup file " + filename + " should be P(bar) or T(K)."));
            }

          for (unsigned int n=2; n<n_columns; ++n)
            {
              in >> column_name;
              if (column_name == "rho,kg/m3")
                prp_indices[0] = n;
              else if (column_name == "phase")
                {
                  has_dominant_phase_column = true;
                  dominant_phase_column_index = n;
                }
              else if (column_name.length() > 3)
                {
                  if (column_name.substr(0,13).compare("vol_fraction_") == 0)
                    {
                      if (std::find(phase_column_names.begin(),
                                    phase_column_names.end(),
                                    column_name) != phase_column_names.end())
                        {
                          AssertThrow(false,
                                      ExcMessage("The PerpleX lookup file " + filename + " must have unique column names. "
                                                 "Sometimes, the same phase is stable with >1 composition at the same "
                                                 "pressure and temperature, so you may see several columns with the same name. "
                                                 "Either combine columns with the same name, or change the names."));
                        }
                      // Populate phase_column_names with the column name
                      // and phase_column_indices with the column index in the current lookup file.
                      phase_column_indices.push_back(n);
                      phase_column_names.push_back(column_name);
                    }
                }
            }
          AssertThrow(std::all_of(prp_indices.begin(), prp_indices.end(), [](int i)
          {
            return i>=0;
          }),
          ExcMessage("The PerpleX lookup file " + filename + " must contain columns with names "
                     "rho,kg/m3, alpha,1/K, cp,J/K/kg, vp,km/s, vs,km/s and h,J/kg."));

          std::getline(in, temp); // first data line

          AssertThrow(min_temp >= 0.0, ExcMessage("Read in of Material header failed (mintemp)."));
          AssertThrow(delta_temp > 0, ExcMessage("Read in of Material header failed (delta_temp)."));
          AssertThrow(n_temperature > 0, ExcMessage("Read in of Material header failed (numtemp)."));
          AssertThrow(min_press >= 0, ExcMessage("Read in of Material header failed (min_press)."));
          AssertThrow(delta_press > 0, ExcMessage("Read in of Material header failed (delta_press)."));
          AssertThrow(n_pressure > 0, ExcMessage("Read in of Material header failed (numpress)."));


          max_temp = min_temp + (n_temperature-1) * delta_temp;
          max_press = min_press + (n_pressure-1) * delta_press;

          density_values.reinit(n_temperature,n_pressure);
          thermal_expansivity_values.reinit(n_temperature,n_pressure);
          specific_heat_values.reinit(n_temperature,n_pressure);
          vp_values.reinit(n_temperature,n_pressure);
          vs_values.reinit(n_temperature,n_pressure);
          enthalpy_values.reinit(n_temperature,n_pressure);

          if (has_dominant_phase_column)
            dominant_phase_indices.reinit(n_temperature,n_pressure);

          phase_volume_fractions.resize(phase_column_names.size());
          for (auto &phase_volume_fraction : phase_volume_fractions)
            phase_volume_fraction.reinit(n_temperature,n_pressure);

          unsigned int i = 0;
          std::vector<double> previous_row_values(n_columns, 0.);

          while (!in.eof())
            {
              std::vector<double> row_values(n_columns);
              std::string phase;

              for (unsigned int n=0; n<n_columns; ++n)
                {
                  if (n == dominant_phase_column_index)
                    in >> phase;
                  else
                    in >> row_values[n]; // assigned as 0 if in.fail() == True

                  // P-T grids created with PerpleX-werami sometimes contain rows
                  // filled with NaNs at extreme P-T conditions where the thermodynamic
                  // models break down. These P-T regions are typically not relevant to
                  // geodynamic modeling (they most commonly appear above
                  // mantle liquidus temperatures at low pressures).
                  // More frustratingly, PerpleX-vertex occasionally fails to find a
                  // valid mineral assemblage in small, isolated regions within the domain,
                  // and so PerpleX-werami also returns NaNs for pixels within these regions.
                  // It is recommended that the user preprocesses their input
                  // files to replace these NaNs before plugging them into ASPECT.
                  // If this lookup encounters invalid doubles it replaces them
                  // with the most recent valid double.
                  if (in.fail())
                    {
                      in.clear();
                      row_values[n] = previous_row_values[n];
                    }
                }
              previous_row_values = row_values;

              std::getline(in, temp); // read next line
              if (in.eof())
                break;

              if (std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase) == dominant_phase_names.end())
                dominant_phase_names.push_back(phase);

              // The ordering of the first two columns in the PerpleX table files
              // dictates whether the inner loop is over temperature or pressure.
              // The first column is always the inner loop.
              // The following lines populate the material property tables
              // according to that implicit loop structure.
              if (first_natural_variable == "T(K)")
                {
                  density_values[i%n_temperature][i/n_temperature]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i%n_temperature][i/n_temperature]= 0.0; // row_values[prp_indices[1]];
                  specific_heat_values[i%n_temperature][i/n_temperature]= 0.0; // row_values[prp_indices[2]];
                  vp_values[i%n_temperature][i/n_temperature]= 0.0; //row_values[prp_indices[3]];
                  vs_values[i%n_temperature][i/n_temperature]= 0.0; //row_values[prp_indices[4]];
                  enthalpy_values[i%n_temperature][i/n_temperature]= 0.0; //row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i%n_temperature][i/n_temperature] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i%n_temperature][i/n_temperature]=row_values[phase_column_indices[n]];
                    }
                }
              else // first_natural_variable == "P(bar)"
                {
                  density_values[i/n_pressure][i%n_pressure]=row_values[prp_indices[0]];
                  thermal_expansivity_values[i/n_pressure][i%n_pressure]=0.0;//row_values[prp_indices[1]];
                  specific_heat_values[i/n_pressure][i%n_pressure]=0.0; //row_values[prp_indices[2]];
                  vp_values[i/n_pressure][i%n_pressure]=0.0; //row_values[prp_indices[3]];
                  vs_values[i/n_pressure][i%n_pressure]=0.0; //row_values[prp_indices[4]];
                  enthalpy_values[i/n_pressure][i%n_pressure]=0.0; //row_values[prp_indices[5]];

                  if (has_dominant_phase_column)
                    {
                      std::vector<std::string>::iterator it = std::find(dominant_phase_names.begin(), dominant_phase_names.end(), phase);
                      dominant_phase_indices[i/n_pressure][i%n_pressure] = std::distance(dominant_phase_names.begin(), it);
                    }

                  for (unsigned int n=0; n<phase_volume_fractions.size(); ++n)
                    {
                      phase_volume_fractions[n][i/n_pressure][i%n_pressure]=row_values[phase_column_indices[n]];
                    }
                }
              ++i;
            }
          AssertThrow(i == n_temperature*n_pressure, ExcMessage("Material table size not consistent with header."));

        }



        void
        EntropyReader::initialize(const MPI_Comm comm,
                                  const std::string &data_directory,
                                  const std::string &material_file_name)
        {
          material_lookup = std::make_unique<Utilities::StructuredDataLookup<2>>(7,1.0);
          material_lookup->load_file(data_directory+material_file_name,
                                     comm);
        }



        double
        EntropyReader::specific_heat(const double entropy,
                                     const double pressure) const
        {
          const double specific_heat = material_lookup->get_data({entropy,pressure}, 3);
          return specific_heat;
        }



        double
        EntropyReader::density(const double entropy,
                               const double pressure) const
        {
          const double density = material_lookup->get_data({entropy,pressure}, 1);
          return density;
        }



        double
        EntropyReader::thermal_expansivity(const double entropy,
                                           const double pressure) const
        {
          const double thermal_expansivity = material_lookup->get_data({entropy,pressure}, 2);
          return thermal_expansivity;
        }



        double
        EntropyReader::temperature(const double entropy,
                                   const double pressure) const
        {
          const double temperature = material_lookup->get_data({entropy,pressure}, 0, true);
          return temperature;
        }



        double
        EntropyReader::seismic_vp(const double entropy,
                                  const double pressure) const
        {
          const double seismic_vp = material_lookup->get_data({entropy,pressure}, 4);
          return seismic_vp;
        }



        double
        EntropyReader::seismic_vs(const double entropy,
                                  const double pressure) const
        {
          const double seismic_vs = material_lookup->get_data({entropy,pressure}, 5);
          return seismic_vs;
        }



        Tensor<1, 2>
        EntropyReader::density_gradient(const double entropy,
                                        const double pressure) const
        {
          const Tensor<1, 2> density_gradient= material_lookup->get_gradients({entropy,pressure}, 1);
          return density_gradient;
        }


      }



      std::vector<double>
      compute_only_composition_fractions(const std::vector<double> &compositional_fields,
                                         const std::vector<unsigned int> &indices_to_use)
      {
        std::vector<double> composition_fractions(indices_to_use.size()+1);

        // Clip the compositional fields so they are between zero and one,
        // and sum the compositional fields for normalization purposes.
        double sum_composition = 0.0;
        std::vector<double> x_comp (indices_to_use.size());

        for (unsigned int i=0; i < x_comp.size(); ++i)
          {
            x_comp[i] = std::min(std::max(compositional_fields[indices_to_use[i]], 0.0), 1.0);
            sum_composition += x_comp[i];
          }

        // Compute background field fraction
        if (sum_composition >= 1.0)
          composition_fractions[0] = 0.0;
        else
          composition_fractions[0] = 1.0 - sum_composition;

        // Compute and possibly normalize field fractions
        for (unsigned int i=0; i < x_comp.size(); ++i)
          {
            if (sum_composition >= 1.0)
              composition_fractions[i+1] = x_comp[i]/sum_composition;
            else
              composition_fractions[i+1] = x_comp[i];
          }

        return composition_fractions;
      }



      std::vector<double>
      compute_composition_fractions(const std::vector<double> &compositional_fields,
                                    const ComponentMask &field_mask)
      {
        std::vector<double> composition_fractions(compositional_fields.size()+1);

        // Clip the compositional fields so they are between zero and one,
        // and sum the compositional fields for normalization purposes.
        double sum_composition = 0.0;
        std::vector<double> x_comp = compositional_fields;
        for (unsigned int i=0; i < x_comp.size(); ++i)
          if (field_mask[i] == true)
            {
              x_comp[i] = std::min(std::max(x_comp[i], 0.0), 1.0);
              sum_composition += x_comp[i];
            }

        // Compute background field fraction
        if (sum_composition >= 1.0)
          composition_fractions[0] = 0.0;
        else
          composition_fractions[0] = 1.0 - sum_composition;

        // Compute and possibly normalize field fractions
        for (unsigned int i=0; i < x_comp.size(); ++i)
          if (field_mask[i] == true)
            {
              if (sum_composition >= 1.0)
                composition_fractions[i+1] = x_comp[i]/sum_composition;
              else
                composition_fractions[i+1] = x_comp[i];
            }

        return composition_fractions;
      }



      std::vector<double>
      compute_volumes_from_masses(const std::vector<double> &masses,
                                  const std::vector<double> &densities,
                                  const bool return_as_fraction)
      {
        Assert(masses.size() == densities.size(),
               ExcMessage ("The mass fractions and densities vectors used for computing "
                           "volumes from masses have to have the same length! "
                           "You have provided "
                           + Utilities::int_to_string(masses.size()) +
                           " mass fractions and "
                           + Utilities::int_to_string(densities.size()) +
                           " densities."));

        const unsigned int n_fields = masses.size();
        std::vector<double> volumes(n_fields);

        if (n_fields == 1 && return_as_fraction)
          {
            volumes[0] = 1.0;
            return volumes;
          }

        double sum_volumes = 0.0;
        for (unsigned int j=0; j < n_fields; ++j)
          {
            volumes[j] = masses[j] / densities[j];
            sum_volumes += volumes[j];
          }

        if (return_as_fraction)
          {
            for (unsigned int j=0; j < n_fields; ++j)
              volumes[j] /= sum_volumes;
          }
        return volumes;
      }



      CompositionalAveragingOperation
      parse_compositional_averaging_operation (const std::string &parameter_name,
                                               const ParameterHandler &prm)
      {
        CompositionalAveragingOperation averaging_operation;
        if (prm.get (parameter_name) == "harmonic")
          averaging_operation = MaterialUtilities::harmonic;
        else if (prm.get (parameter_name) == "arithmetic")
          averaging_operation = MaterialUtilities::arithmetic;
        else if (prm.get (parameter_name) == "geometric")
          averaging_operation = MaterialUtilities::geometric;
        else if (prm.get (parameter_name) == "maximum composition")
          averaging_operation = MaterialUtilities::maximum_composition;
        else
          {
            AssertThrow(false, ExcMessage("Not a valid viscosity averaging scheme"));

            //We will never get here, but we have to return something so the compiler does not complain
            return MaterialUtilities::harmonic;
          }

        return averaging_operation;
      }



      double
      average_value (const std::vector<double> &volume_fractions,
                     const std::vector<double> &parameter_values,
                     const enum CompositionalAveragingOperation &average_type)
      {
        Assert(volume_fractions.size() == parameter_values.size(),
               ExcMessage ("The volume fractions and parameter values vectors used for averaging "
                           "have to have the same length! You have provided "
                           + Utilities::int_to_string(volume_fractions.size()) +
                           " volume fractions and "
                           + Utilities::int_to_string(parameter_values.size()) +
                           " parameter values."));

        double averaged_parameter = 0.0;

        switch (average_type)
          {
            case arithmetic:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                averaged_parameter += volume_fractions[i] * parameter_values[i];
              break;
            }
            case harmonic:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                {
                  //AssertThrow(parameter_values[i] > 0,
                  //            ExcMessage ("All parameter values must be greater than 0 for harmonic averaging!"));
                  // todo
                  if (parameter_values[i] < 0)
                    {
                      std::cout << "Warning: Utilities: average_value: there is a negative value of parameters" << std::endl;
                    }
                  averaged_parameter += volume_fractions[i]/(std::max(parameter_values[i], 1e-8));
                }
              averaged_parameter = 1.0/averaged_parameter;
              break;
            }
            case geometric:
            {
              for (unsigned int i=0; i<volume_fractions.size(); ++i)
                {
                  AssertThrow(parameter_values[i] > 0,
                              ExcMessage ("All parameter values must be greater than 0 for geometric averaging!"));
                  averaged_parameter += volume_fractions[i] * std::log(parameter_values[i]);
                }
              averaged_parameter = std::exp(averaged_parameter);
              break;
            }
            case maximum_composition:
            {
              const unsigned int idx = static_cast<unsigned int>(std::max_element( volume_fractions.begin(),
                                                                                   volume_fractions.end() )
                                                                 - volume_fractions.begin());
              averaged_parameter = parameter_values[idx];
              break;
            }
            default:
            {
              AssertThrow(false, ExcNotImplemented());
              break;
            }
          }
        return averaged_parameter;
      }



      template <int dim>
      void
      fill_averaged_equation_of_state_outputs(const EquationOfStateOutputs<dim> &eos_outputs,
                                              const std::vector<double> &mass_fractions,
                                              const std::vector<double> &volume_fractions,
                                              const unsigned int i,
                                              MaterialModelOutputs<dim> &out)
      {
        // The density, isothermal compressibility and thermal expansivity are volume-averaged
        // The specific entropy derivatives and heat capacity are mass-averaged
        out.densities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.densities, MaterialUtilities::arithmetic);
        out.compressibilities[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.compressibilities, MaterialUtilities::arithmetic);
        out.thermal_expansion_coefficients[i] = MaterialUtilities::average_value (volume_fractions, eos_outputs.thermal_expansion_coefficients, MaterialUtilities::arithmetic);
        out.entropy_derivative_pressure[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_pressure, MaterialUtilities::arithmetic);
        out.entropy_derivative_temperature[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.entropy_derivative_temperature, MaterialUtilities::arithmetic);
        out.specific_heat[i] = MaterialUtilities::average_value (mass_fractions, eos_outputs.specific_heat_capacities, MaterialUtilities::arithmetic);
      }



      double phase_average_value (const std::vector<double> &phase_function_values,
                                  const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                  const std::vector<double> &parameter_values,
                                  const unsigned int composition_index,
                                  const PhaseUtilities::PhaseAveragingOperation operation)
      {
        // Calculate base index and assign base value
        unsigned int start_phase_index = 0;
        for (unsigned int i=0; i<composition_index; ++i)
          start_phase_index += n_phase_transitions_per_composition[i] + 1;

        double averaged_parameter = parameter_values[start_phase_index];
        if (n_phase_transitions_per_composition[composition_index] > 0)
          {
            // Do averaging when there are multiple phases
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = std::log(averaged_parameter);

            for (unsigned int i=0; i<n_phase_transitions_per_composition[composition_index]; ++i)
              {
                const unsigned int phase_index = start_phase_index + i;

                Assert(phase_index+1<parameter_values.size(), ExcInternalError());
                if (operation == PhaseUtilities::logarithmic)
                  {
                    // First average by log values and then take the exponential.
                    // This is used for averaging prefactors in flow laws.
                    averaged_parameter += phase_function_values[phase_index-composition_index] * (std::log(parameter_values[phase_index+1]) - averaged_parameter);
                  }
                else if (operation == PhaseUtilities::arithmetic)
                  averaged_parameter += phase_function_values[phase_index-composition_index] * (parameter_values[phase_index+1] - averaged_parameter);

                else
                  AssertThrow(false, ExcInternalError());
              }
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = std::exp(averaged_parameter);
          }
        return averaged_parameter;
      }



      double phase_average_value1 (const std::vector<double> &phase_function_values,
                                   const std::vector<unsigned int> &n_phase_transitions_per_composition,
                                   const std::vector<double> &parameter_values,
                                   const unsigned int composition_index,
                                   const PhaseUtilities::PhaseAveragingOperation operation)
      {
        // Calculate base index and assign base value
        unsigned int start_phase_index = 0;
        for (unsigned int i=0; i<composition_index; ++i)
          start_phase_index += n_phase_transitions_per_composition[i] + 1;

        double averaged_parameter = parameter_values[start_phase_index];
        if (n_phase_transitions_per_composition[composition_index] > 0)
          {
            // Do averaging when there are multiple phases
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = log(averaged_parameter);

            for (unsigned int i=0; i<n_phase_transitions_per_composition[composition_index]; ++i)
              {
                const unsigned int phase_index = start_phase_index + i;

                Assert(phase_index+1<parameter_values.size(), ExcInternalError());
                if (operation == PhaseUtilities::logarithmic)
                  {
                    // First average by log values and then take the exponential.
                    // This is used for averaging prefactors in flow laws.
                    averaged_parameter += phase_function_values[phase_index-composition_index] * log(parameter_values[phase_index+1] / parameter_values[phase_index]);
                  }
                else if (operation == PhaseUtilities::arithmetic)
                  averaged_parameter += phase_function_values[phase_index-composition_index] * (parameter_values[phase_index+1] - parameter_values[phase_index]);

                else
                  AssertThrow(false, ExcInternalError());
              }
            if (operation == PhaseUtilities::logarithmic)
              averaged_parameter = exp(averaged_parameter);
          }
        return averaged_parameter;
      }

      template <int dim>
      PhaseFunctionInputs<dim>::PhaseFunctionInputs(const double temperature_,
                                                    const double pressure_,
                                                    const double depth_,
                                                    const double pressure_depth_derivative_,
                                                    const unsigned int phase_transition_index_)

        :
        temperature(temperature_),
        pressure(pressure_),
        depth(depth_),
        pressure_depth_derivative(pressure_depth_derivative_),
        phase_transition_index(phase_transition_index_)
      {}

      template <int dim>
      PhaseFunctionInputs1<dim>::PhaseFunctionInputs1(const double temperature_,
                                                      const double pressure_,
                                                      const double depth_,
                                                      const double pressure_depth_derivative_,
                                                      const unsigned int phase_transition_index_,
                                                      const std::vector<double> &composition)

        :
        temperature(temperature_),
        pressure(pressure_),
        depth(depth_),
        pressure_depth_derivative(pressure_depth_derivative_),
        phase_transition_index(phase_transition_index_),
        composition(std::make_shared<std::vector<double>>(composition))
      {}

      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::initialize()
      {
        AssertThrow (this->introspection().get_number_of_fields_of_type(CompositionalFieldDescription::chemical_composition)+1 == material_file_names.size(),
                     ExcMessage(" The 'discrete phase function' requires that the number of compositional fields of type chemical composition plus one (for a background field) "
                                "matches the number of lookup tables."));


        minimum_temperature = std::vector<double>(material_file_names.size());
        maximum_temperature = std::vector<double>(material_file_names.size());
        interval_temperature = std::vector<double>(material_file_names.size());
        minimum_pressure = std::vector<double>(material_file_names.size());
        maximum_pressure = std::vector<double>(material_file_names.size());
        interval_pressure = std::vector<double>(material_file_names.size());
        for (unsigned int i = 0; i < material_file_names.size(); ++i)
          {
            material_lookup
            .push_back(std::make_unique<Utilities::StructuredDataLookup<2>>(8,1.0));

            material_lookup[i]->load_file(data_directory+material_file_names[i],
                                          this->get_mpi_communicator());

            Assert(material_lookup[i]->has_equidistant_coordinates(), ExcMessage("The loaded lookup tables do not use equidistant coordinates. The class 'PhaseFunctionDiscrete' cannot currently handle non equidistant coordinates."));

            minimum_temperature[i] = material_lookup[i]->get_interpolation_point_coordinates(0).front();
            maximum_temperature[i] = material_lookup[i]->get_interpolation_point_coordinates(0).back();
            interval_temperature[i] = (maximum_temperature[i]-minimum_temperature[i])/(material_lookup[i]->get_number_of_coordinates(0)-1);
            minimum_pressure[i] = material_lookup[i]->get_interpolation_point_coordinates(1).front();
            maximum_pressure[i] = material_lookup[i]->get_interpolation_point_coordinates(1).back();
            interval_pressure[i] = (maximum_pressure[i]-minimum_pressure[i])/(material_lookup[i]->get_number_of_coordinates(1)-1);
          }
      }


      template <int dim>
      double
      PhaseFunctionDiscrete<dim>::compute_value (const PhaseFunctionInputs<dim> &in) const
      {
        // the percentage of material that has undergone the transition
        double function_value;

        // lookup the most abundant phase
        unsigned int start_phase_transition_index = 0;
        unsigned int n_comp = 0;
        for (unsigned int n_relevant_fields  = 0 ; n_relevant_fields  < this->introspection().n_chemical_composition_fields() + 1 ; n_relevant_fields ++)
          {
            if (in.phase_transition_index < start_phase_transition_index + n_phase_transitions_per_chemical_composition[n_relevant_fields])
              {
                n_comp = n_relevant_fields ;
                break;
              }
            start_phase_transition_index += n_phase_transitions_per_chemical_composition[n_relevant_fields];
          }

        const double pressure_in_bar = in.pressure/1.e5;

        Assert (in.temperature >= minimum_temperature[n_comp] && in.temperature < maximum_temperature[n_comp], ExcInternalError());
        Assert (pressure_in_bar >= minimum_pressure[n_comp] && pressure_in_bar < maximum_pressure[n_comp], ExcInternalError());

        const std::vector<double> &temperature_points = material_lookup[n_comp]->get_interpolation_point_coordinates(0);
        const std::vector<double> &pressure_points = material_lookup[n_comp]->get_interpolation_point_coordinates(1);

        // round T and p points to exact values in the table
        const unsigned int i_T = static_cast<unsigned int>(std::round((in.temperature - minimum_temperature[n_comp]) / interval_temperature[n_comp]));
        const unsigned int i_p = static_cast<unsigned int>(std::round((pressure_in_bar - minimum_pressure[n_comp]) / interval_pressure[n_comp]));

        Point<2> temperature_pressure(temperature_points[i_T], pressure_points[i_p]);

        // determine the dominant phase index
        const unsigned int dominant_phase_index = static_cast<unsigned int>(std::round(material_lookup[n_comp]->get_data(temperature_pressure, 7)));

        // match the dominant phase index to one of the transition indicators
        unsigned int matched_phase_transition_index = numbers::invalid_unsigned_int;
        for (unsigned int i = start_phase_transition_index;
             i < start_phase_transition_index + n_phase_transitions_per_chemical_composition[n_comp];
             i++)
          {
            if (transition_indicators[i] == dominant_phase_index)
              {
                matched_phase_transition_index = i;
              }
          }

        // determine the value of phase function to facilitate the exact transition
        if ((matched_phase_transition_index != numbers::invalid_unsigned_int) && in.phase_transition_index <= matched_phase_transition_index)
          function_value = 1.0;
        else
          function_value = 0.0;

        return function_value;
      }


      template <int dim>
      double
      PhaseFunctionDiscrete<dim>::compute_derivative () const
      {
        // raises an error to ensure that a phase derivative request is not made for this phase function.
        AssertThrow(false, ExcNotImplemented());
      }


      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phase_transitions () const
      {
        return transition_indicators.size();
      }



      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phases () const
      {
        return n_phases_total;
      }



      template <int dim>
      unsigned int
      PhaseFunctionDiscrete<dim>::
      n_phases_over_all_chemical_compositions () const
      {
        return n_phases_total_chemical_compositions;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phase_transitions_for_each_composition () const
      {
        return *n_phase_transitions_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phases_for_each_composition () const
      {
        return n_phases_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phase_transitions_for_each_chemical_composition () const
      {
        return n_phase_transitions_per_chemical_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunctionDiscrete<dim>::n_phases_for_each_chemical_composition () const
      {
        return n_phases_per_chemical_composition;
      }



      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Phase transition indicators", "",
                           Patterns::Anything(),
                           "A list of phase indicators in a look-up table for each phase transition. "
                           "This parameter selectively assign different rheologies to specific phases, "
                           "rather than having a unique rheology for each phase in the table. "
                           "For example, if the table has phases 0, 1, and 2, and one only want "
                           "a distinct rheology for phase 2, then only phase 2 is needed "
                           "in the list of indicator. And phases 0, 1 will just be assigned "
                           "the rheology of the base phase. ");

        prm.declare_entry ("Data directory", "$ASPECT_SOURCE_DIR/data/material-model/entropy-table/pyrtable",
                           Patterns::DirectoryName (),
                           "The path to the model data. The path may also include the special "
                           "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
                           "in which the ASPECT source files were located when ASPECT was "
                           "compiled. This interpretation allows, for example, to reference "
                           "files located in the `data/' subdirectory of ASPECT. ");

        prm.declare_entry ("Material file names", "material_table_temperature_pressure_small.txt",
                           Patterns::List (Patterns::Anything()),
                           "The file names of the material data (material "
                           "data is assumed to be in order with the ordering "
                           "of the compositional fields). Note that there are "
                           "two options on how many files need to be listed "
                           "here: 1. If only one file is provided, it is used "
                           "for the whole model domain, and compositional fields "
                           "are ignored. 2. If there is one more file name than the "
                           "number of compositional fields, then the first file is "
                           "assumed to define a `background composition' that is "
                           "modified by the compositional fields. These data files need "
                           "to have the same structure as the one necessary for equation "
                           "of state plus a new column for the phase indexes, which amounts "
                           "to 8 columns in total.");
      }



      template <int dim>
      void
      PhaseFunctionDiscrete<dim>::parse_parameters (ParameterHandler &prm)
      {
        // Retrieve the list of composition names
        const std::vector<std::string> list_of_composition_names = this->introspection().get_composition_names();

        data_directory               = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
        material_file_names          = Utilities::split_string_list(prm.get ("Material file names"));

        // Make options file for parsing maps to double arrays
        std::vector<std::string> chemical_field_names = this->introspection().chemical_composition_field_names();
        chemical_field_names.insert(chemical_field_names.begin(),"background");

        Utilities::MapParsing::Options options(chemical_field_names, "Phase transition indicators");
        options.allow_missing_keys = true;
        options.allow_multiple_values_per_key = true;
        options.store_values_per_key = true;
        options.n_values_per_key = std::vector<unsigned int>();

        const std::vector<double> transition_indicators_double          = Utilities::MapParsing::parse_map_to_double_array (prm.get("Phase transition indicators"),
                                                                          options);

        transition_indicators.reserve(transition_indicators_double.size());
        for (const double transition_indicator: transition_indicators_double)
          transition_indicators.push_back(static_cast<unsigned int>(std::round(transition_indicator)));
        n_phase_transitions_per_composition = std::make_unique<std::vector<unsigned int>>(options.n_values_per_key);

        n_phases_total = 0;
        n_phases_per_composition.clear();
        for (unsigned int n : *n_phase_transitions_per_composition)
          {
            n_phases_per_composition.push_back(n+1);
            n_phases_total += n+1;
          }

        // The background field is always the first composition
        n_phases_per_chemical_composition = {n_phases_per_composition[0]};
        n_phase_transitions_per_chemical_composition = {n_phases_per_composition[0] - 1};
        n_phases_total_chemical_compositions = n_phases_per_composition[0];
        for (auto i : this->introspection().chemical_composition_field_indices())
          {
            n_phases_per_chemical_composition.push_back(n_phases_per_composition[i+1]);
            n_phase_transitions_per_chemical_composition.push_back(n_phases_per_composition[i+1] - 1);
            n_phases_total_chemical_compositions += n_phases_per_composition[i+1];
          }
      }


      template <int dim>
      double
      PhaseFunction<dim>::compute_value (const PhaseFunctionInputs<dim> &in) const
      {
        AssertIndexRange (in.phase_transition_index, transition_temperature_lower_limits.size());
        AssertIndexRange (in.phase_transition_index, transition_temperature_upper_limits.size());

        // the percentage of material that has undergone the transition
        double function_value;
        if (in.temperature < transition_temperature_lower_limits[in.phase_transition_index] ||
            in.temperature >= transition_temperature_upper_limits[in.phase_transition_index])
          {
            // assign 0.0 if temperature is out of range
            function_value = 0.0;
          }
        else
          {
            if (use_depth_instead_of_pressure)
              {
                AssertIndexRange (in.phase_transition_index, transition_depths.size());

                // calculate the deviation from the transition point (convert temperature to depth)
                double depth_deviation = in.depth - transition_depths[in.phase_transition_index];

                if (in.pressure_depth_derivative != 0.0)
                  {
                    AssertIndexRange (in.phase_transition_index, transition_slopes.size());
                    AssertIndexRange (in.phase_transition_index, transition_temperatures.size());

                    depth_deviation -= transition_slopes[in.phase_transition_index] / in.pressure_depth_derivative
                                       * (in.temperature - transition_temperatures[in.phase_transition_index]);
                  }

                // use delta function for width = 0
                AssertIndexRange (in.phase_transition_index, transition_widths.size());
                if (transition_widths[in.phase_transition_index] == 0)
                  function_value = (depth_deviation > 0) ? 1. : 0.;
                else
                  function_value = 0.5*(1.0 + std::tanh(depth_deviation / transition_widths[in.phase_transition_index]));
              }
            else
              {
                // calculate the deviation from the transition point (convert temperature to pressure)
                AssertIndexRange (in.phase_transition_index, transition_pressures.size());
                const double pressure_deviation = in.pressure - transition_pressures[in.phase_transition_index]
                                                  - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

                // use delta function for width = 0
                AssertIndexRange (in.phase_transition_index, transition_pressure_widths.size());
                if (transition_pressure_widths[in.phase_transition_index] == 0)
                  function_value = (pressure_deviation > 0) ? 1. : 0.;
                else
                  function_value = 0.5*(1.0 + std::tanh(pressure_deviation / transition_pressure_widths[in.phase_transition_index]));
              }
          }
        return function_value;
      }

      template <int dim>
      double
      PhaseFunction<dim>::compute_value1 (const PhaseFunctionInputs1<dim> &in) const
      {
        // the percentage of material that has undergone the transition
        double function_value;
        double use_manually_method_for_spcrust = manually_method_crust[in.phase_transition_index];
        double use_manually_method_for_pyrolite = manually_method_pyrolite[in.phase_transition_index];
        double use_manually_method_for_harzburgite = manually_method_harzburgite[in.phase_transition_index];
        if ( abs(use_manually_method_for_spcrust - 1.0) < 1e-8)
          {
            function_value = eclogite_transition.compute_value_crust_1_0(in, manually_method_crust,
                                                                         transition_depths, transition_temperatures,
                                                                         transition_widths, transition_slopes);
          }
        else if ( abs(use_manually_method_for_spcrust - 1.1) < 1e-8)
          {
            function_value = eclogite_transition.compute_value_crust_1_1(in, manually_method_crust,
                                                                         transition_depths, transition_temperatures,
                                                                         transition_widths, transition_slopes);
          }
        else if ( abs(use_manually_method_for_spcrust - 1.2) < 1e-8)
          {
            function_value = eclogite_transition.compute_value_crust_1_2(in, manually_method_crust,
                                                                         transition_depths, transition_temperatures,
                                                                         transition_widths, transition_slopes);
          }
        else if ( abs(use_manually_method_for_spcrust - 1.3) < 1e-8)
          {
            function_value = eclogite_transition.compute_value_crust_1_3(in, manually_method_crust,
                                                                         transition_depths, transition_temperatures,
                                                                         transition_widths, transition_slopes);
          }
        else if ( abs(use_manually_method_for_pyrolite - 1.0) < 1e-8)
          {
            function_value = pyrolite_transition.compute_value_pyrolite_1_0(in, manually_method_pyrolite,
                                                                            transition_depths, transition_temperatures,
                                                                            transition_widths, transition_slopes);
            if (is_using_metastable_kinetics)
              {
                unsigned int metastable_index = this->introspection().compositional_index_for_name("metastable");
                if (this->is_beyond_equilibrium(in.phase_transition_index))
                  {
                    function_value *= (*in.composition)[metastable_index];
                  }
              }
          }
        else if ( abs(use_manually_method_for_pyrolite - 1.1) < 1e-8)
          {
            function_value = pyrolite_transition.compute_value_pyrolite_1_1(in, manually_method_pyrolite,
                                                                            transition_depths, transition_temperatures,
                                                                            transition_widths, transition_slopes);
          }
        else if ( abs(use_manually_method_for_harzburgite - 1.0) < 1e-8)
          {
            function_value = pyrolite_transition.compute_value_harzburgite_1_0(in, manually_method_harzburgite,
                                                                               transition_depths, transition_temperatures,
                                                                               transition_widths, transition_slopes);

            if (is_using_metastable_kinetics)
              {
                unsigned int metastable_index = this->introspection().compositional_index_for_name("metastable");
                if (this->is_beyond_equilibrium(in.phase_transition_index))
                  {
                    function_value *= (*in.composition)[metastable_index];
                  }
              }
          }
        else if (use_depth_instead_of_pressure)
          {
            // calculate the deviation from the transition point (convert temperature to depth)
            double depth_deviation = in.depth - transition_depths[in.phase_transition_index];

            if (in.pressure_depth_derivative != 0.0)
              depth_deviation -= transition_slopes[in.phase_transition_index] / in.pressure_depth_derivative
                                 * (in.temperature - transition_temperatures[in.phase_transition_index]);

            // use delta function for width = 0
            if (transition_widths[in.phase_transition_index] == 0)
              function_value = (depth_deviation > 0) ? 1. : 0.;
            else
              function_value = 0.5*(1.0 + std::tanh(depth_deviation / transition_widths[in.phase_transition_index]));

            // use lower limits and upper limits to restrict the region of phase transition
            if (in.depth < transition_depth_lower_limits[in.phase_transition_index])
              function_value = 0.0;
            else if (in.depth > transition_depth_upper_limits[in.phase_transition_index])
              function_value = 1.0;
          }
        else
          {
            // calculate the deviation from the transition point (convert temperature to pressure)
            const double pressure_deviation = in.pressure - transition_pressures[in.phase_transition_index]
                                              - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

            // use delta function for width = 0
            if (transition_pressure_widths[in.phase_transition_index] == 0)
              function_value = (pressure_deviation > 0) ? 1. : 0.;
            else
              function_value = 0.5*(1.0 + std::tanh(pressure_deviation / transition_pressure_widths[in.phase_transition_index]));

            // use lower limits and upper limits to restrict the region of phase transition
            if (in.pressure < transition_pressure_lower_limits[in.phase_transition_index])
              function_value = 0.0;
            else if (in.depth > transition_pressure_upper_limits[in.phase_transition_index])
              function_value = 1.0;
          }
        return function_value;
      }

      template <int dim>
      bool
      PhaseFunction<dim>::is_beyond_equilibrium (const unsigned int phase_transition_index) const
      {
        bool beyond_metastable_transition = false;
        unsigned int start_phase_transition_index = 0;
        unsigned int composition_index = 0;
        for (unsigned int i=0; i<n_phase_transitions_per_composition->size(); ++i)
          {
            if ((phase_transition_index >= start_phase_transition_index) && (phase_transition_index < start_phase_transition_index + (*n_phase_transitions_per_composition)[i]))
              {
                composition_index = i;
                break;
              }
            start_phase_transition_index += (*n_phase_transitions_per_composition)[i];
          }

        for (unsigned int phase_transition_index_temp = start_phase_transition_index; phase_transition_index_temp <= phase_transition_index; phase_transition_index_temp++)
          {
            if (use_metastable_kinetics[phase_transition_index_temp])
              {
                beyond_metastable_transition = true;
                break;
              }
          }
        return beyond_metastable_transition;
      }


      template <int dim>
      bool
      PhaseFunction<dim>::is_beyond_decomposition (const unsigned int phase_transition_index) const
      {
        bool beyond_metastable_decomposition = false;
        unsigned int start_phase_transition_index = 0;
        unsigned int composition_index = 0;
        for (unsigned int i=0; i<n_phase_transitions_per_composition->size(); ++i)
          {
            if ((phase_transition_index >= start_phase_transition_index) && (phase_transition_index < start_phase_transition_index + (*n_phase_transitions_per_composition)[i]))
              {
                composition_index = i;
                break;
              }
            start_phase_transition_index += (*n_phase_transitions_per_composition)[i];
          }

        for (unsigned int phase_transition_index_temp = start_phase_transition_index; phase_transition_index_temp <= phase_transition_index; phase_transition_index_temp++)
          {
            if (use_metastable_decompositions[phase_transition_index_temp])
              {
                beyond_metastable_decomposition = true;
                break;
              }
          }
        return beyond_metastable_decomposition;
      }


      template <int dim>
      bool
      PhaseFunction<dim>::is_metastable_phase (const unsigned int phase_index) const
      {
        unsigned int start_phase_index = 0;
        unsigned int composition_index = 0;
        for (unsigned int i=0; i< n_phases_per_composition.size(); ++i)
          {
            if ((start_phase_index <= phase_index) && (phase_index < start_phase_index + n_phases_per_composition[i]))
              {
                composition_index = i;
                break;
              }
            start_phase_index += n_phases_per_composition[i];
          }
        if (phase_index == start_phase_index)
          {
            // no transition yet
            return false;
          }
        else
          {
            // has transition, need the transition index that lead to this phase
            // e.g for phase 1, look into transition 0
            unsigned phase_transition_index = phase_index - composition_index - 1;
            if (this->is_beyond_equilibrium(phase_transition_index) and !this->is_beyond_decomposition(phase_transition_index))
              return true;
            else
              return false;
          }
      }

      template <int dim>
      double
      PhaseFunction<dim>::compute_derivative (const PhaseFunctionInputs<dim> &in) const
      {
        double transition_pressure;
        double pressure_width;
        double width_temp;

        // we already should have the adiabatic conditions here
        Assert (this->get_adiabatic_conditions().is_initialized(),
                ExcMessage("The adiabatic conditions need to be already initialized "
                           "to calculate the derivative of phase functions. Either call this "
                           "function after the reference conditions have been computed, or implement "
                           "a workaround for the case without reference profile."));

        // phase transition based on depth
        if (use_depth_instead_of_pressure)
          {
            const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index]);
            const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] + transition_widths[in.phase_transition_index]);
            const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] - transition_widths[in.phase_transition_index]);
            transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
            pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                    - this->get_adiabatic_conditions().pressure(transition_minus_width));
            width_temp = transition_widths[in.phase_transition_index];
          }
        // using pressure instead of depth to define the phase transition
        else
          {
            transition_pressure = transition_pressures[in.phase_transition_index];
            pressure_width = transition_pressure_widths[in.phase_transition_index];
            width_temp = transition_pressure_widths[in.phase_transition_index];
          }

        // calculate the deviation from the transition point
        const double pressure_deviation = in.pressure - transition_pressure
                                          - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

        // calculate the analytical derivative of the phase function
        if (
          (in.temperature < transition_temperature_lower_limits[in.phase_transition_index]) ||
          (in.temperature >= transition_temperature_upper_limits[in.phase_transition_index])
        )
          {
            // return 0 if temperature is out of range
            return 0;
          }
        else if (width_temp == 0)
          return 0;
        else
          return 0.5 / pressure_width * (1.0 - std::tanh(pressure_deviation / pressure_width)
                                         * std::tanh(pressure_deviation / pressure_width));
      }

      template <int dim>
      double
      PhaseFunction<dim>::compute_derivative1 (const PhaseFunctionInputs1<dim> &in) const
      {
        double transition_pressure;
        double pressure_width;
        double width_temp;

        // we already should have the adiabatic conditions here
        Assert (this->get_adiabatic_conditions().is_initialized(),
                ExcMessage("The adiabatic conditions need to be already initialized "
                           "to calculate the derivative of phase functions. Either call this "
                           "function after the reference conditions have been computed, or implement "
                           "a workaround for the case without reference profile."));

        // phase transition based on depth
        if (use_depth_instead_of_pressure)
          {
            const Point<dim,double> transition_point = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index]);
            const Point<dim,double> transition_plus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] + transition_widths[in.phase_transition_index]);
            const Point<dim,double> transition_minus_width = this->get_geometry_model().representative_point(transition_depths[in.phase_transition_index] - transition_widths[in.phase_transition_index]);
            transition_pressure = this->get_adiabatic_conditions().pressure(transition_point);
            pressure_width = 0.5 * (this->get_adiabatic_conditions().pressure(transition_plus_width)
                                    - this->get_adiabatic_conditions().pressure(transition_minus_width));
            width_temp = transition_widths[in.phase_transition_index];
          }
        // using pressure instead of depth to define the phase transition
        else
          {
            transition_pressure = transition_pressures[in.phase_transition_index];
            pressure_width = transition_pressure_widths[in.phase_transition_index];
            width_temp = transition_pressure_widths[in.phase_transition_index];
          }

        // calculate the deviation from the transition point
        const double pressure_deviation = in.pressure - transition_pressure
                                          - transition_slopes[in.phase_transition_index] * (in.temperature - transition_temperatures[in.phase_transition_index]);

        // calculate the analytical derivative of the phase function
        if (
          (in.temperature < transition_temperature_lower_limits[in.phase_transition_index]) ||
          (in.temperature >= transition_temperature_upper_limits[in.phase_transition_index])
        )
          {
            // return 0 if temperature is out of range
            return 0;
          }
        else if (width_temp == 0)
          return 0;
        else
          {

            double derivative = 0.5 / pressure_width * (1.0 - std::tanh(pressure_deviation / pressure_width)
                                                        * std::tanh(pressure_deviation / pressure_width));

            if (is_using_metastable_kinetics)
              {
                unsigned int metastable_index = this->introspection().compositional_index_for_name("metastable");
                if (this->is_beyond_equilibrium(in.phase_transition_index))
                  {
                    derivative *= (*in.composition)[metastable_index];
                  }
              }

            return derivative;
          }
      }



      template <int dim>
      unsigned int
      PhaseFunction<dim>::
      n_phase_transitions () const
      {
        if (use_depth_instead_of_pressure)
          return transition_depths.size();
        else
          return transition_pressures.size();
      }



      template <int dim>
      unsigned int
      PhaseFunction<dim>::
      n_phases () const
      {
        return n_phases_total;
      }



      template <int dim>
      unsigned int
      PhaseFunction<dim>::
      n_phases_over_all_chemical_compositions () const
      {
        return n_phases_total_chemical_compositions;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phase_transitions_for_each_composition () const
      {
        return *n_phase_transitions_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phases_for_each_composition () const
      {
        return n_phases_per_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phase_transitions_for_each_chemical_composition () const
      {
        return n_phase_transitions_per_chemical_composition;
      }



      template <int dim>
      const std::vector<unsigned int> &
      PhaseFunction<dim>::n_phases_for_each_chemical_composition () const
      {
        return n_phases_per_chemical_composition;
      }



      template <int dim>
      double
      PhaseFunction<dim>::
      get_transition_slope (const unsigned int phase_transition_index) const
      {
        return transition_slopes[phase_transition_index];
      }

      template <int dim>
      double
      PhaseFunction<dim>::
      get_compute_latent_heat(const unsigned int phase_transition_index) const
      {
        return compute_latent_heats[phase_transition_index];
      }

      template <int dim>
      double
      PhaseFunction<dim>::
      get_transition_depth (const unsigned int phase_transition_index) const
      {
        return transition_depths[phase_transition_index];
      }



      template <int dim>
      void
      PhaseFunction<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Phase transition depths", "",
                           Patterns::Anything(),
                           "A list of depths where phase transitions occur. Values must "
                           "monotonically increase. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition widths", "",
                           Patterns::Anything(),
                           "A list of widths for each phase transition, in terms of depth. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition pressures", "",
                           Patterns::Anything(),
                           "A list of pressures where phase transitions occur. Values must "
                           "monotonically increase. Define transition by depth instead of "
                           "pressure must be set to false to use this parameter. "
                           "Units: \\si{\\pascal}.");
        prm.declare_entry ("Phase transition pressure widths", "",
                           Patterns::Anything(),
                           "A list of widths for each phase transition, in terms of pressure. The phase functions "
                           "are scaled with these values, leading to a jump between phases "
                           "for a value of zero and a gradual transition for larger values. "
                           "List must have the same number of entries as Phase transition pressures. "
                           "Define transition by depth instead of pressure must be set to false "
                           "to use this parameter. "
                           "Units: \\si{\\pascal}.");
        prm.declare_entry ("Define transition by depth instead of pressure", "true",
                           Patterns::Bool (),
                           "Whether to list phase transitions by depth or pressure. If this parameter is true, "
                           "then the input file will use Phase transitions depths and Phase transition widths "
                           "to define the phase transition. If it is false, the parameter file will read in "
                           "phase transition data from Phase transition pressures and "
                           "Phase transition pressure widths.");
        prm.declare_entry ("Phase transition temperatures", "",
                           Patterns::Anything(),
                           "A list of temperatures where phase transitions occur. Higher or lower "
                           "temperatures lead to phase transition occurring in smaller or greater "
                           "depths than given in Phase transition depths, depending on the "
                           "Clapeyron slope given in Phase transition Clapeyron slopes. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition temperature upper limits",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::max()),
                           Patterns::Anything(),
                           "A list of upper temperature limits for each phase transition. Above "
                           "this temperature the respective phase transition is deactivated. The default "
                           "value means there is no upper limit for any phase transitions. "
                           "List must have the same number of entries as Phase transition depths. "
                           "When the optional temperature limits are applied, the user has to be "
                           "careful about the consistency between adjacent phases. Phase transitions "
                           "should be continuous in pressure-temperature space. "
                           "We recommend producing a phase diagram with "
                           "simple model setups to check the implementation as a starting point."
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition temperature lower limits",
                           boost::lexical_cast<std::string>(std::numeric_limits<double>::lowest()),
                           Patterns::Anything(),
                           "A list of lower temperature limits for each phase transition. Below "
                           "this temperature the respective phase transition is deactivated. The default "
                           "value means there is no lower limit for any phase transition. "
                           "List must have the same number of entries as Phase transition depths. "
                           "When the optional temperature limits are applied, the user has to be "
                           "careful about the consistency between adjacent phases. Phase transitions "
                           "should be continuous in pressure-temperature space. "
                           "We recommend producing a phase diagram with "
                           "simple model setups to check the implementation as a starting point."
                           "Units: \\si{\\kelvin}.");
        prm.declare_entry ("Phase transition Clapeyron slopes", "",
                           Patterns::Anything(),
                           "A list of Clapeyron slopes for each phase transition. A positive "
                           "Clapeyron slope indicates that the phase transition will occur in "
                           "a greater depth, if the temperature is higher than the one given in "
                           "Phase transition temperatures and in a smaller depth, if the "
                           "temperature is smaller than the one given in Phase transition temperatures. "
                           "For negative slopes the other way round. "
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\pascal\\per\\kelvin}.");
        prm.declare_entry ("Phase transition depth lower limits", "-1e16",
                           Patterns::Anything(),
                           "A list of limits for each phase transition, in terms of depth. The phase transitions "
                           "only happen at deeper region"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition depth upper limits", "1e16",
                           Patterns::Anything(),
                           "A list of limits for each phase transition, in terms of depth. The phase transitions "
                           "only happen at shallower region"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\meter}.");
        prm.declare_entry ("Phase transition pressure lower limits", "-1e16",
                           Patterns::Anything(),
                           "A list of limits for each phase transition, in terms of pressure. The phase transitions "
                           "only happen at deeper region"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\pascal}.");
        prm.declare_entry ("Phase transition pressure upper limits", "1e16",
                           Patterns::Anything(),
                           "A list of limits for each phase transition, in terms of pressure. The phase transitions "
                           "only happen at shallower region"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: \\si{\\pascal}.");
        // define the manually defined composition
        prm.declare_entry ("Manually define phase method crust", "0.0",
                           Patterns::Anything(),
                           "A list of version of method to use for each phase transition for crust"
                           "version numbers are like 1.0, 1.1, 1.2 ..."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");
        // define the manually defined composition for pyrolite
        prm.declare_entry ("Manually define phase method pyrolite", "0.0",
                           Patterns::Anything(),
                           "A list of version of method to use for each phase transition for pyrolite"
                           "version numbers are like 1.0, 1.1, 1.2 ..."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");
        // define the manually defined composition for pyrolite
        prm.declare_entry ("Manually define phase method harzburgite", "0.0",
                           Patterns::Anything(),
                           "A list of version of method to use for each phase transition for pyrolite"
                           "version numbers are like 1.0, 1.1, 1.2 ..."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");
        // compute latent heat on phases
        prm.declare_entry ("Compute latent heat", "1.0",
                           Patterns::Anything(),
                           "A list of int, indicating whether to compute latent heat on this phase transition"
                           "Entries are either 0.0 or 1.0"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");

        prm.declare_entry ("Metastable transition", "0.0",
                           Patterns::Anything(),
                           "Whether to compute phase transitions by metastable kinetics."
                           "Entries are either 0.0 or 1.0."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");

        prm.declare_entry ("Metastable decomposition", "0.0",
                           Patterns::Anything(),
                           "Whether to include decomposition reaction and mark the index of"
                           "decomposition phase in metastable kinetics."
                           "Entries are either 0.0 or 1.0."
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");

        prm.declare_entry ("Metastable transition comp", "0",
                           Patterns::Anything(),
                           "Compute phase transitions by metastable kinetics and record the metastable composition"
                           "Entries are indexes of compositoin"
                           "List must have the same number of entries as Phase transition depths. "
                           "Units: None.");

        // declare parameters for eclogite_transition
        EclogiteTransition<dim>::declare_parameters(prm);
      }



      template <int dim>
      void
      PhaseFunction<dim>::parse_parameters (ParameterHandler &prm)
      {
        std::vector<std::string> compositional_field_names = this->introspection().get_composition_names();
        // Establish that a background field is required here
        compositional_field_names.insert(compositional_field_names.begin(),"background");
        Utilities::MapParsing::Options options(compositional_field_names, "Phase transition temperatures");
        options.allow_missing_keys = true;
        options.allow_multiple_values_per_key = true;
        options.store_values_per_key = true;
        options.n_values_per_key = std::vector<unsigned int>();

        options.property_name = "Phase transition temperatures";
        transition_temperatures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        // Now we know how many entries we have for each composition
        // ensure that the number of phase transitions is consistent
        options.store_values_per_key = false;
        options.check_values_per_key = true;

        use_depth_instead_of_pressure = prm.get_bool ("Define transition by depth instead of pressure");

        if (use_depth_instead_of_pressure)
          {
            options.property_name = "Phase transition depths";
            transition_depths          = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition widths";
            transition_widths          = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition depth lower limits";
            transition_depth_lower_limits         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition depth upper limits";
            transition_depth_upper_limits         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            // parse the manually defined composition
            options.property_name = "Manually define phase method crust";
            manually_method_crust         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            // parse the manually defined composition for pyrolite
            options.property_name = "Manually define phase method pyrolite";
            manually_method_pyrolite         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            // parse the manually defined composition for harzburgite
            options.property_name = "Manually define phase method harzburgite";
            manually_method_harzburgite        = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            // parse the Compute latent heat
            options.property_name = "Compute latent heat";
            compute_latent_heats        =  Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Metastable transition";
            std::vector<double> use_metastable_kinetics_double    = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            use_metastable_kinetics.resize(use_metastable_kinetics_double.size());
            std::transform(use_metastable_kinetics_double.begin(), use_metastable_kinetics_double.end(), use_metastable_kinetics.begin(),
                           [](double val)
            {
              return static_cast<bool>(std::round(val));
            });

            options.property_name = "Metastable decomposition";
            std::vector<double>    use_metastable_decompositions_double = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            use_metastable_decompositions.resize(use_metastable_decompositions_double.size());
            std::transform(use_metastable_decompositions_double.begin(), use_metastable_decompositions_double.end(), use_metastable_decompositions.begin(),
                           [](double val)
            {
              return static_cast<bool>(std::round(val));
            });

            is_using_metastable_kinetics = std::any_of(
                                             use_metastable_kinetics.begin(),
                                             use_metastable_kinetics.end(),
                                             [](bool value)
            {
              return value;  // Lambda function to check if value is true
            }
                                           );

            // parse A value for the eclogite transition temperature
            eclogite_transition.parse_parameters(prm);
          }
        else
          {
            options.property_name = "Phase transition pressures";
            transition_pressures = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition pressure widths";
            transition_pressure_widths = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition pressure lower limits";
            transition_pressure_lower_limits         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

            options.property_name = "Phase transition pressure upper limits";
            transition_pressure_upper_limits         = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);
          }

        options.property_name = "Phase transition temperature upper limits";
        transition_temperature_upper_limits = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Phase transition temperature lower limits";
        transition_temperature_lower_limits = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        options.property_name = "Phase transition Clapeyron slopes";
        transition_slopes = Utilities::MapParsing::parse_map_to_double_array (prm.get(options.property_name), options);

        n_phase_transitions_per_composition = std::make_unique<std::vector<unsigned int>>(options.n_values_per_key);

        n_phases_total = 0;
        n_phases_per_composition.clear();
        for (unsigned int n : *n_phase_transitions_per_composition)
          {
            n_phases_per_composition.push_back(n+1);
            n_phases_total += n+1;
          }

        // The background field is always the first composition
        n_phases_per_chemical_composition = {n_phases_per_composition[0]};
        n_phase_transitions_per_chemical_composition = {n_phases_per_composition[0] - 1};
        n_phases_total_chemical_compositions = n_phases_per_composition[0];
        for (auto i : this->introspection().chemical_composition_field_indices())
          {
            n_phases_per_chemical_composition.push_back(n_phases_per_composition[i+1]);
            n_phase_transitions_per_chemical_composition.push_back(n_phases_per_composition[i+1] - 1);
            n_phases_total_chemical_compositions += n_phases_per_composition[i+1];
          }
      }

      template <int dim>
      bool
      PhaseFunction<dim>::get_is_using_metastable_kinetics() const
      {
        return is_using_metastable_kinetics;
      }

      template <int dim>
      void
      EclogiteTransition<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Eclogite transition");
        {
          // declare A value for the eclogite transition temperature
          prm.declare_entry ("Temperature for eclogite transition", "973.0", Patterns::Double (),
                             "The temperature for crustal phase transition");
          prm.declare_entry ("Temperature width for eclogite transition", "75.0", Patterns::Double (),
                             "The width of temperature for crustal phase transition");
          prm.declare_entry ("Temperature slope for eclogite transition", "1e10", Patterns::Double (),
                             "The clapeyron slope of temperature for crustal phase transition");
          prm.declare_entry ("Pressure for eclogite transition", "1.5e9", Patterns::Double (),
                             "The pressure for crustal phase transition");
          prm.declare_entry ("Pressure width for eclogite transition", "0.5e9", Patterns::Double (),
                             "The width of pressure for crustal phase transition");
          prm.declare_entry ("Pressure slope for eclogite transition", "0.0", Patterns::Double (),
                             "The pressure slope for crustal phase transition");
          prm.declare_entry ("Max pressure for eclogite transition", "5e9", Patterns::Double (),
                             "The maximum pressure for crustal phase transition."
                             "This helps to force the transition in very cold region");
          prm.declare_entry ("Max pressure width for eclogite transition", "0.5e9", Patterns::Double (),
                             "The width of maximum pressure for crustal phase transition.");
          prm.declare_entry ("Average phase functions for eclogite transition",
                             "true", Patterns::Bool (),
                             "If the phase functions from the pressure and temperature boundaries are averaged for eclogite transition");
        }
        prm.leave_subsection();
      }

      template <int dim>
      void
      EclogiteTransition<dim>::parse_parameters (ParameterHandler &prm)
      {
        prm.enter_subsection ("Eclogite transition");
        {
          crust_eclogite_transition_T     =  Utilities::string_to_double(prm.get("Temperature for eclogite transition"));
          crust_eclogite_transition_T_width     =  Utilities::string_to_double(prm.get("Temperature width for eclogite transition"));
          crust_eclogite_transition_T_slope     =  Utilities::string_to_double(prm.get("Temperature slope for eclogite transition"));
          crust_eclogite_transition_P     =  Utilities::string_to_double(prm.get("Pressure for eclogite transition"));
          crust_eclogite_transition_P_slope     =  Utilities::string_to_double(prm.get("Pressure slope for eclogite transition"));
          crust_eclogite_transition_P_width     =  Utilities::string_to_double(prm.get("Pressure width for eclogite transition"));
          crust_eclogite_transition_max_P = Utilities::string_to_double(prm.get("Max pressure for eclogite transition"));
          crust_eclogite_transition_max_P_width = Utilities::string_to_double(prm.get("Max pressure width for eclogite transition"));
          crust_eclogite_transition_PT_average = prm.get_bool("Average phase functions for eclogite transition");
        }
        prm.leave_subsection();
      }


      template <int dim>
      double
      EclogiteTransition<dim>::compute_value_crust_1_0 (const PhaseFunctionInputs1<dim> &in,
                                                        const std::vector<double> &manually_method_crust,
                                                        const std::vector<double> &transition_depths,
                                                        const std::vector<double> &transition_temperatures,
                                                        const std::vector<double> &transition_widths,
                                                        const std::vector<double> &transition_slopes) const
      {
        const double version = 1.0;
        // version 1.0
        double function_value = 0.0;
        int phase_transition_index_crust = 0;
        // composition-wise index
        while ( abs(manually_method_crust[in.phase_transition_index - phase_transition_index_crust - 1] - version) < 1e-8)
          phase_transition_index_crust++;
        // find a region in a phase diagram
        const double P0 = 1.50e9; // Pa
        std::pair<bool, double> result0 = compute_point_to_line(in, 0.0, P0, 0.0, 0.0, false, false, false);

        // define ecologite transition by temperature
        const double W1 = 75.0;
        // const double T1 = 1048.0; // K
        const double T1 = crust_eclogite_transition_T + W1;  // as what we need is the dash line
        std::pair<bool, double> result1 = compute_point_to_line(in, T1, 0.0, W1, 0.0, false, false, true);

        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_crust + 2;  // third one
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);
        // std::cout << d660 << T660 << W660 << slope660/in.pressure_depth_derivative << in.pressure << in.temperature << std::endl;

        if (result0.first && result1.first && (!result660.first))
          {
            // crustal eclogite transition
            if (phase_transition_index_crust == 0)
              {
                function_value = 0.5*(1.0 + std::tanh(result1.second/W1));
              }
            else
              function_value = 0.0;
          }
        else if ( (!result1.first) && result660.first)
          {
            // 660 for mantle
            function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        else if ( result1.first && result660.first)
          {
            // 660 for crust
            if (phase_transition_index_crust == 0)
              function_value = 1.0;
            else
              function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        else
          {
            // phase 0
            function_value = 0.0;
          }
        return function_value;
      }

      template <int dim>
      double
      EclogiteTransition<dim>::compute_value_crust_1_1 (const PhaseFunctionInputs1<dim> &in,
                                                        const std::vector<double> &manually_method_crust,
                                                        const std::vector<double> &transition_depths,
                                                        const std::vector<double> &transition_temperatures,
                                                        const std::vector<double> &transition_widths,
                                                        const std::vector<double> &transition_slopes) const
      {
        // version 1.1
        const double version = 1.1;
        double function_value = 0.0;
        int phase_transition_index_crust = 0;
        // composition-wise index
        while ( abs(manually_method_crust[in.phase_transition_index - phase_transition_index_crust - 1] - version) < 1e-8)
          phase_transition_index_crust++;
        // find a region in a phase diagram
        const double W0 = crust_eclogite_transition_P_width;
        const double P0 = crust_eclogite_transition_P + W0; // Pa
        std::pair<bool, double> result0 = compute_point_to_line(in, 0.0, P0, W0, 0.0, false, false, false);

        // define ecologite transition by temperature
        const double W1 = crust_eclogite_transition_T_width;
        // const double T1 = 1048.0; // K
        const double T1 = crust_eclogite_transition_T + W1;  // as what we need is the dash line
        std::pair<bool, double> result1 = compute_point_to_line(in, T1, 0.0, W1, 0.0, false, false, true);

        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_crust + 2;  // third one
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);
        // std::cout << d660 << T660 << W660 << slope660/in.pressure_depth_derivative << in.pressure << in.temperature << std::endl;

        if (result0.first && result1.first && (!result660.first))
          {
            // crustal eclogite transition
            // const double deviation = std::min(result0.second/W0, result1.second/W1);
            const double deviation = (result0.second/W0 + result1.second/W1) / 2.0;
            if (phase_transition_index_crust == 0)
              {
                if (true)
                  function_value = 0.5*(1.0 + std::tanh(deviation));
                else
                  {
                    if (deviation > 0.0)
                      function_value = 1.0;
                    else
                      function_value = 0.5*(2.0 + deviation);
                  }
              }
            else
              function_value = 0.0;
          }
        else if ( (!result1.first) && result660.first)
          {
            // 660 for mantle
            function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        else if ( result1.first && result660.first)
          {
            // 660 for crust
            if (phase_transition_index_crust == 0)
              function_value = 1.0;
            else
              function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        else
          {
            // phase 0
            function_value = 0.0;
          }
        return function_value;
      }

      template <int dim>
      double
      EclogiteTransition<dim>::compute_value_crust_1_2 (const PhaseFunctionInputs1<dim> &in,
                                                        const std::vector<double> &manually_method_crust,
                                                        const std::vector<double> &transition_depths,
                                                        const std::vector<double> &transition_temperatures,
                                                        const std::vector<double> &transition_widths,
                                                        const std::vector<double> &transition_slopes) const
      {
        // version 1.2
        const double version = 1.2;
        // paritial_indexes
        const int partial_index_660 = 1;
        // initiate
        double function_value = 0.0;
        int phase_transition_index_crust = 0;
        // composition-wise index
        while ( abs(manually_method_crust[in.phase_transition_index - phase_transition_index_crust - 1] - version) < 1e-8)
          phase_transition_index_crust++;
        // find a region in a phase diagram
        const double W0 = crust_eclogite_transition_P_width;
        const double P0 = crust_eclogite_transition_P + W0; // Pa
        std::pair<bool, double> result0 = compute_point_to_line(in, 0.0, P0, W0, 0.0, false, false, false);

        // define ecologite transition by temperature
        // add a slope
        double W1;
        const double T1 = crust_eclogite_transition_T + crust_eclogite_transition_T_width;  // as what we need is the dash line
        std::pair<bool, double>  result1;
        if (abs(crust_eclogite_transition_T_slope) > 1e9)
          {
            // vertical
            W1 = crust_eclogite_transition_T_width;
            result1 = compute_point_to_line(in, T1, 0.0, W1, 0.0, false, false, true);
          }
        else
          {
            // with a slope, pinpoint at (T1, P0), as W1 is a width by temperature, it is multiplied with slope
            W1 = crust_eclogite_transition_T_width * abs(crust_eclogite_transition_T_slope);
            result1 = compute_point_to_line(in, T1, crust_eclogite_transition_P, W1,
                                            crust_eclogite_transition_T_slope, false, false, false);
          }

        // line 2: maximux pressure on basaltic composition
        const double P2 = crust_eclogite_transition_max_P;
        const double W2 = crust_eclogite_transition_max_P_width;
        std::pair<bool, double> result2 = compute_point_to_line(in, 0.0, P2, W2, 0.0, false, false, false);

        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_crust + partial_index_660;  // second one
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);

        if (result0.first && result1.first && (!result660.first))
          {
            // crustal eclogite transition
            // double deviation = (result0.second/W0 + result1.second/W1) / 2.0;
            // deviation = std::max(result2.second/W2, deviation);
            double deviation = average_deviation(result0.second/W0, std::max(result1.second/W1, result2.second/W2), 2.0);
            if (phase_transition_index_crust == 0)
              {
                function_value = 0.5*(1.0 + std::tanh(deviation));
              }
            else
              function_value = 0.0;
          }
        else if ( (!result1.first) && result2.first && (!result660.first))
          {
            // crustal eclogite transition: area 2 (line 0 and line 1)
            const double deviation = result2.second/W2;
            if (phase_transition_index_crust == 0)
              {
                if (true)
                  function_value = 0.5*(1.0 + std::tanh(deviation));
                else
                  {
                    if (deviation > 0.0)
                      function_value = 1.0;
                    else
                      function_value = 0.5*(2.0 + deviation);
                  }
              }
            else
              function_value = 0.0;
          }
        else if ( result660.first)
          {
            // 660 for crust
            if (phase_transition_index_crust == 0)
              function_value = 1.0;
            else
              function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        else
          {
            // phase 0
            function_value = 0.0;
          }
        return function_value;
      }


      template <int dim>
      double
      EclogiteTransition<dim>::compute_value_crust_1_3 (const PhaseFunctionInputs1<dim> &in,
                                                        const std::vector<double> &manually_method_crust,
                                                        const std::vector<double> &transition_depths,
                                                        const std::vector<double> &transition_temperatures,
                                                        const std::vector<double> &transition_widths,
                                                        const std::vector<double> &transition_slopes) const
      {
        // version 1.3
        const double version = 1.3;
        // paritial_indexes
        const int partial_index_660 = 1;
        const int partial_index_720 = 2;
        // initiate
        double function_value = 0.0;
        int phase_transition_index_crust = 0;
        // composition-wise index
        while ( abs(manually_method_crust[in.phase_transition_index - phase_transition_index_crust - 1] - version) < 1e-8)
          phase_transition_index_crust++;
        // find a region in a phase diagram
        // define the boundary by pressure
        const double W0 = crust_eclogite_transition_P_width;
        const double P0 = crust_eclogite_transition_P + W0; // Pa
        std::pair<bool, double> result0 = compute_point_to_line(in, 1150.0, P0, W0,
                                                                crust_eclogite_transition_P_slope, false, false, false);

        // define ecologite transition by temperature
        // add a slope
        double W1;
        const double T1 = crust_eclogite_transition_T + crust_eclogite_transition_T_width;  // as what we need is the dash line
        std::pair<bool, double>  result1;
        if (abs(crust_eclogite_transition_T_slope) > 1e9)
          {
            // vertical
            W1 = crust_eclogite_transition_T_width;
            result1 = compute_point_to_line(in, T1, 0.0, W1, 0.0, false, false, true);
          }
        else
          {
            // with a slope, pinpoint at (T1, P0), as W1 is a width by temperature, it is multiplied with slope
            W1 = crust_eclogite_transition_T_width * abs(crust_eclogite_transition_T_slope);
            result1 = compute_point_to_line(in, T1, crust_eclogite_transition_P, W1,
                                            crust_eclogite_transition_T_slope, false, false, false);
          }

        // line 2: maximux pressure on basaltic composition
        const double P2 = crust_eclogite_transition_max_P;
        const double W2 = crust_eclogite_transition_max_P_width;
        std::pair<bool, double> result2 = compute_point_to_line(in, 0.0, P2, W2, 0.0, false, false, false);

        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_crust + partial_index_660;  // second one
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);

        const int phase_transition_index_720 = in.phase_transition_index - phase_transition_index_crust + partial_index_720;  // third one
        const double d720 = transition_depths[phase_transition_index_720];
        const double T720 = transition_temperatures[phase_transition_index_720];
        const double W720 = transition_widths[phase_transition_index_720];
        const double slope720 = transition_slopes[phase_transition_index_720];
        std::pair<bool, double> result720 = compute_point_to_line(in, T720, d720, W720, slope720/in.pressure_depth_derivative, true, false, false);

        if (result0.first && result1.first && (!result660.first))
          {
            // crustal eclogite transition
            // double deviation = (result0.second/W0 + result1.second/W1) / 2.0;
            // deviation = std::max(result2.second/W2, deviation);
            double deviation;
            if (crust_eclogite_transition_PT_average)
              {
                deviation = average_deviation(result0.second/W0, std::max(result1.second/W1, result2.second/W2), 2.0);
              }
            else
              {
                deviation = std::min(result0.second/W0, std::max(result1.second/W1, result2.second/W2));
              }
            if (phase_transition_index_crust == 0)
              {
                function_value = 0.5*(1.0 + std::tanh(deviation));
              }
            else
              function_value = 0.0;
          }
        else if ( (!result1.first) && result2.first && (!result660.first))
          {
            // crustal eclogite transition: area 2 (line 0 and line 1)
            const double deviation = result2.second/W2;
            if (phase_transition_index_crust == 0)
              {
                if (true)
                  function_value = 0.5*(1.0 + std::tanh(deviation));
                else
                  {
                    if (deviation > 0.0)
                      function_value = 1.0;
                    else
                      function_value = 0.5*(2.0 + deviation);
                  }
              }
            else
              function_value = 0.0;
          }
        else if (result660.first && !result720.first)
          {
            // 660 for crust
            if (phase_transition_index_crust < partial_index_660)
              function_value = 1.0;
            else if (phase_transition_index_crust == partial_index_660)
              function_value = 0.5*(1.0 + std::tanh(result660.second/W660));
            else
              function_value = 0.0;
          }
        else if (result720.first)
          {
            // 720 for crust
            if (phase_transition_index_crust < partial_index_720)
              function_value = 1.0;
            else if (phase_transition_index_crust == partial_index_720)
              function_value = 0.5*(1.0 + std::tanh(result720.second/W720));
            else
              function_value = 0.0;
          }
        else
          {
            // phase 0
            function_value = 0.0;
          }
        return function_value;
      }


      template <int dim>
      double
      PyroliteTransition<dim>::compute_value_pyrolite_1_0 (const PhaseFunctionInputs1<dim> &in,
                                                           const std::vector<double> &manually_method_pyrolite,
                                                           const std::vector<double> &transition_depths,
                                                           const std::vector<double> &transition_temperatures,
                                                           const std::vector<double> &transition_widths,
                                                           const std::vector<double> &transition_slopes) const
      {
        // version 1.0
        const double version = 1.0;

        // partial indexes of transitions
        const int partial_index_410 = 0;
        const int partial_index_520 = 1;
        const int partial_index_560 = 2;
        const int partial_index_660 = 3;
        const int partial_index_660_gt = 4;
        const int partial_index_660_gt1 = 5;
        const int partial_index_660_gt_combined = 6;

        // initiate varibles
        double function_value = 0.0;
        int phase_transition_index_pyrolite = 0;

        // composition-wise index

        // loop to get the local index relative to the 0th pyrolite phase
        // debug
        while ( in.phase_transition_index - phase_transition_index_pyrolite != 0)
          {
            // see if we reach the start of the pyrolite phases
            // as for the 0th phase tran in the pyrolite phases, this loop is false initially
            if (abs(manually_method_pyrolite[in.phase_transition_index - phase_transition_index_pyrolite - 1] - version) > 1e-8)
              break;
            // add one to the relative index within the pyrolite phases if we haven't
            phase_transition_index_pyrolite++;
          }

        // 410
        const int phase_transition_index_410 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_410;
        const double d410 = transition_depths[phase_transition_index_410];
        const double T410 = transition_temperatures[phase_transition_index_410];
        const double W410 = transition_widths[phase_transition_index_410];
        const double slope410 = transition_slopes[phase_transition_index_410];
        std::pair<bool, double> result410 = compute_point_to_line(in, T410, d410, W410, slope410/in.pressure_depth_derivative, true, false, false);

        // 520
        const int phase_transition_index_520 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_520;
        const double d520 = transition_depths[phase_transition_index_520];
        const double T520 = transition_temperatures[phase_transition_index_520];
        const double W520 = transition_widths[phase_transition_index_520];
        const double slope520 = transition_slopes[phase_transition_index_520];
        std::pair<bool, double> result520 = compute_point_to_line(in, T520, d520, W520, slope520/in.pressure_depth_derivative, true, false, false);

        // 560
        const int phase_transition_index_560 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_560;
        const double d560 = transition_depths[phase_transition_index_560];
        const double T560 = transition_temperatures[phase_transition_index_560];
        const double W560 = transition_widths[phase_transition_index_560];
        const double slope560 = transition_slopes[phase_transition_index_560];
        std::pair<bool, double> result560 = compute_point_to_line(in, T560, d560, W560, slope560/in.pressure_depth_derivative, true, false, false);

        // 660
        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660;
        //const int phase_transition_index_660 = 0;
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 0
        const int phase_transition_index_660_gt = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt;
        //const int phase_transition_index_660 = 0;
        const double d660_gt = transition_depths[phase_transition_index_660_gt];
        const double T660_gt = transition_temperatures[phase_transition_index_660_gt];
        const double W660_gt = transition_widths[phase_transition_index_660_gt];
        const double slope660_gt = transition_slopes[phase_transition_index_660_gt];
        std::pair<bool, double> result660_gt = compute_point_to_line(in, T660_gt, d660_gt, W660_gt, slope660_gt/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 1
        const int phase_transition_index_660_gt1 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt1;
        //const int phase_transition_index_660 = 0;
        const double d660_gt1 = transition_depths[phase_transition_index_660_gt1];
        const double T660_gt1 = transition_temperatures[phase_transition_index_660_gt1];
        const double W660_gt1 = transition_widths[phase_transition_index_660_gt1];
        const double slope660_gt1 = transition_slopes[phase_transition_index_660_gt1];
        std::pair<bool, double> result660_gt1 = compute_point_to_line(in, T660_gt1, d660_gt1, W660_gt1, slope660_gt1/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, combined
        const int phase_transition_index_660_gt_combined = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt_combined;
        //const int phase_transition_index_660 = 0;
        const double d660_gt_combined = transition_depths[phase_transition_index_660_gt_combined];
        const double T660_gt_combined = transition_temperatures[phase_transition_index_660_gt_combined];
        const double W660_gt_combined = transition_widths[phase_transition_index_660_gt_combined];
        const double slope660_gt_combined = transition_slopes[phase_transition_index_660_gt_combined];
        std::pair<bool, double> result660_gt_combined = compute_point_to_line(in, T660_gt_combined, d660_gt_combined, W660_gt_combined, slope660_gt_combined/in.pressure_depth_derivative, true, false, false);

        if (result410.first)
          {
            // 410 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_410)
              function_value += 0.5*(1.0 + std::tanh(result410.second/W410));
          }
        if (result520.first)
          {
            // 520 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_520)
              function_value += 0.5*(1.0 + std::tanh(result520.second/W520));
          }
        if (result560.first)
          {
            // 560 for pyrolite, Gt -> CaPv + Gt
            if (phase_transition_index_pyrolite == partial_index_560)
              function_value += 0.5*(1.0 + std::tanh(result560.second/W560));
          }
        if (result660.first)
          {
            // 660 for pyrolite, rw -> brg + fp
            if (phase_transition_index_pyrolite == partial_index_660)
              function_value += 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        if (result660_gt.first && in.temperature < T660_gt)
          {
            // 660 for pyrolite, gt -> il
            if (phase_transition_index_pyrolite == partial_index_660_gt)
              function_value += 0.5*(1.0 + std::tanh(result660_gt.second/W660_gt));
          }
        if (result660_gt1.first && in.temperature < T660_gt1)
          {
            // 660 for pyrolite, il -> brg
            if (phase_transition_index_pyrolite == partial_index_660_gt1)
              function_value += 0.5*(1.0 + std::tanh(result660_gt1.second/W660_gt1));
          }
        if (result660_gt_combined.first && in.temperature >= T660_gt_combined)
          {
            // 660 for pyrolite combined, at higher temperature, gt -> brg
            if (phase_transition_index_pyrolite == partial_index_660_gt_combined)
              function_value += 0.5*(1.0 + std::tanh(result660_gt_combined.second/W660_gt_combined));
          }
        return function_value;
      }

      template <int dim>
      double
      PyroliteTransition<dim>::compute_value_pyrolite_1_1 (const PhaseFunctionInputs1<dim> &in,
                                                           const std::vector<double> &manually_method_pyrolite,
                                                           const std::vector<double> &transition_depths,
                                                           const std::vector<double> &transition_temperatures,
                                                           const std::vector<double> &transition_widths,
                                                           const std::vector<double> &transition_slopes) const
      {
        // version 1.1
        // in this version, I have adapted 2 parts for the gt transition around 670 km instead of the previous
        // 3 parts implementation.
        const double version = 1.1;

        // partial indexes of transitions
        const int partial_index_410 = 0;
        const int partial_index_520 = 1;
        const int partial_index_560 = 2;
        const int partial_index_660 = 3;
        const int partial_index_660_gt = 4;
        const int partial_index_660_gt1 = 5;
        const int partial_index_660_gt_combined = 6;

        // initiate varibles
        double function_value = 0.0;
        int phase_transition_index_pyrolite = 0;

        // composition-wise index

        // loop to get the local index relative to the 0th pyrolite phase
        // debug
        while ( in.phase_transition_index - phase_transition_index_pyrolite != 0)
          {
            // see if we reach the start of the pyrolite phases
            // as for the 0th phase tran in the pyrolite phases, this loop is false initially
            if (abs(manually_method_pyrolite[in.phase_transition_index - phase_transition_index_pyrolite - 1] - version) > 1e-8)
              break;
            // add one to the relative index within the pyrolite phases if we haven't
            phase_transition_index_pyrolite++;
          }

        // 410
        const int phase_transition_index_410 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_410;
        const double d410 = transition_depths[phase_transition_index_410];
        const double T410 = transition_temperatures[phase_transition_index_410];
        const double W410 = transition_widths[phase_transition_index_410];
        const double slope410 = transition_slopes[phase_transition_index_410];
        std::pair<bool, double> result410 = compute_point_to_line(in, T410, d410, W410, slope410/in.pressure_depth_derivative, true, false, false);

        // 520
        const int phase_transition_index_520 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_520;
        const double d520 = transition_depths[phase_transition_index_520];
        const double T520 = transition_temperatures[phase_transition_index_520];
        const double W520 = transition_widths[phase_transition_index_520];
        const double slope520 = transition_slopes[phase_transition_index_520];
        std::pair<bool, double> result520 = compute_point_to_line(in, T520, d520, W520, slope520/in.pressure_depth_derivative, true, false, false);

        // 560
        const int phase_transition_index_560 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_560;
        const double d560 = transition_depths[phase_transition_index_560];
        const double T560 = transition_temperatures[phase_transition_index_560];
        const double W560 = transition_widths[phase_transition_index_560];
        const double slope560 = transition_slopes[phase_transition_index_560];
        std::pair<bool, double> result560 = compute_point_to_line(in, T560, d560, W560, slope560/in.pressure_depth_derivative, true, false, false);

        // 660
        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660;
        //const int phase_transition_index_660 = 0;
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 0
        const int phase_transition_index_660_gt = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt;
        //const int phase_transition_index_660 = 0;
        const double d660_gt = transition_depths[phase_transition_index_660_gt];
        const double T660_gt = transition_temperatures[phase_transition_index_660_gt];
        const double W660_gt = transition_widths[phase_transition_index_660_gt];
        const double slope660_gt = transition_slopes[phase_transition_index_660_gt];
        std::pair<bool, double> result660_gt = compute_point_to_line(in, T660_gt, d660_gt, W660_gt, slope660_gt/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 1
        const int phase_transition_index_660_gt1 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt1;
        //const int phase_transition_index_660 = 0;
        const double d660_gt1 = transition_depths[phase_transition_index_660_gt1];
        const double T660_gt1 = transition_temperatures[phase_transition_index_660_gt1];
        const double W660_gt1 = transition_widths[phase_transition_index_660_gt1];
        const double slope660_gt1 = transition_slopes[phase_transition_index_660_gt1];
        std::pair<bool, double> result660_gt1 = compute_point_to_line(in, T660_gt1, d660_gt1, W660_gt1, slope660_gt1/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, combined
        const int phase_transition_index_660_gt_combined = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt_combined;
        //const int phase_transition_index_660 = 0;
        const double d660_gt_combined = transition_depths[phase_transition_index_660_gt_combined];
        const double T660_gt_combined = transition_temperatures[phase_transition_index_660_gt_combined];
        const double W660_gt_combined = transition_widths[phase_transition_index_660_gt_combined];
        const double slope660_gt_combined = transition_slopes[phase_transition_index_660_gt_combined];
        std::pair<bool, double> result660_gt_combined = compute_point_to_line(in, T660_gt_combined, d660_gt_combined, W660_gt_combined, slope660_gt_combined/in.pressure_depth_derivative, true, false, false);

        if (result410.first)
          {
            // 410 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_410)
              function_value += 0.5*(1.0 + std::tanh(result410.second/W410));
          }
        if (result520.first)
          {
            // 520 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_520)
              function_value += 0.5*(1.0 + std::tanh(result520.second/W520));
          }
        if (result560.first)
          {
            // 560 for pyrolite, Gt -> CaPv + Gt
            if (phase_transition_index_pyrolite == partial_index_560)
              function_value += 0.5*(1.0 + std::tanh(result560.second/W560));
          }
        if (result660.first)
          {
            // 660 for pyrolite, rw -> brg + fp
            if (phase_transition_index_pyrolite == partial_index_660)
              function_value += 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        if (result660_gt.first && in.temperature < T660_gt)
          {
            // 660 for pyrolite, part 1, colder part
            if (phase_transition_index_pyrolite == partial_index_660_gt)
              function_value += 0.5*(1.0 + std::tanh(result660_gt.second/W660_gt));
          }
        if (result660_gt1.first && in.temperature > T660_gt1)
          {
            // 660 for pyrolite, part 2, hotter part
            if (phase_transition_index_pyrolite == partial_index_660_gt1)
              function_value += 0.5*(1.0 + std::tanh(result660_gt1.second/W660_gt1));
          }
        if (result660_gt_combined.first)
          {
            // 660 for pyrolite combined, at higher temperature, gt -> brg
            if (phase_transition_index_pyrolite == partial_index_660_gt_combined)
              function_value += 0.5*(1.0 + std::tanh(result660_gt_combined.second/W660_gt_combined));
          }
        return function_value;
      }

      template <int dim>
      double
      PyroliteTransition<dim>::compute_value_harzburgite_1_0 (const PhaseFunctionInputs1<dim> &in,
                                                              const std::vector<double> &manually_method_harzburgite,
                                                              const std::vector<double> &transition_depths,
                                                              const std::vector<double> &transition_temperatures,
                                                              const std::vector<double> &transition_widths,
                                                              const std::vector<double> &transition_slopes) const
      {
        // version 1.0
        const double version = 1.0;

        // partial indexes of transitions
        const int partial_index_410 = 0;
        const int partial_index_520 = 1;
        const int partial_index_560 = 2;
        const int partial_index_660 = 3;
        const int partial_index_660_gt = 4;
        const int partial_index_660_gt1 = 5;
        const int partial_index_660_gt_combined = 6;

        // initiate varibles
        double function_value = 0.0;
        int phase_transition_index_pyrolite = 0;

        // composition-wise index

        // loop to get the local index relative to the 0th pyrolite phase
        // debug
        while ( in.phase_transition_index - phase_transition_index_pyrolite != 0)
          {
            // see if we reach the start of the pyrolite phases
            // as for the 0th phase tran in the pyrolite phases, this loop is false initially
            if (abs(manually_method_harzburgite[in.phase_transition_index - phase_transition_index_pyrolite - 1] - version) > 1e-8)
              break;
            // add one to the relative index within the pyrolite phases if we haven't
            phase_transition_index_pyrolite++;
          }

        // 410
        const int phase_transition_index_410 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_410;
        const double d410 = transition_depths[phase_transition_index_410];
        const double T410 = transition_temperatures[phase_transition_index_410];
        const double W410 = transition_widths[phase_transition_index_410];
        const double slope410 = transition_slopes[phase_transition_index_410];
        std::pair<bool, double> result410 = compute_point_to_line(in, T410, d410, W410, slope410/in.pressure_depth_derivative, true, false, false);

        // 520
        const int phase_transition_index_520 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_520;
        const double d520 = transition_depths[phase_transition_index_520];
        const double T520 = transition_temperatures[phase_transition_index_520];
        const double W520 = transition_widths[phase_transition_index_520];
        const double slope520 = transition_slopes[phase_transition_index_520];
        std::pair<bool, double> result520 = compute_point_to_line(in, T520, d520, W520, slope520/in.pressure_depth_derivative, true, false, false);

        // 560
        const int phase_transition_index_560 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_560;
        const double d560 = transition_depths[phase_transition_index_560];
        const double T560 = transition_temperatures[phase_transition_index_560];
        const double W560 = transition_widths[phase_transition_index_560];
        const double slope560 = transition_slopes[phase_transition_index_560];
        std::pair<bool, double> result560 = compute_point_to_line(in, T560, d560, W560, slope560/in.pressure_depth_derivative, true, false, false);

        // 660
        const int phase_transition_index_660 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660;
        //const int phase_transition_index_660 = 0;
        const double d660 = transition_depths[phase_transition_index_660];
        const double T660 = transition_temperatures[phase_transition_index_660];
        const double W660 = transition_widths[phase_transition_index_660];
        const double slope660 = transition_slopes[phase_transition_index_660];
        std::pair<bool, double> result660 = compute_point_to_line(in, T660, d660, W660, slope660/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 0
        const int phase_transition_index_660_gt = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt;
        //const int phase_transition_index_660 = 0;
        const double d660_gt = transition_depths[phase_transition_index_660_gt];
        const double T660_gt = transition_temperatures[phase_transition_index_660_gt];
        const double W660_gt = transition_widths[phase_transition_index_660_gt];
        const double slope660_gt = transition_slopes[phase_transition_index_660_gt];
        std::pair<bool, double> result660_gt = compute_point_to_line(in, T660_gt, d660_gt, W660_gt, slope660_gt/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, part 1
        const int phase_transition_index_660_gt1 = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt1;
        //const int phase_transition_index_660 = 0;
        const double d660_gt1 = transition_depths[phase_transition_index_660_gt1];
        const double T660_gt1 = transition_temperatures[phase_transition_index_660_gt1];
        const double W660_gt1 = transition_widths[phase_transition_index_660_gt1];
        const double slope660_gt1 = transition_slopes[phase_transition_index_660_gt1];
        std::pair<bool, double> result660_gt1 = compute_point_to_line(in, T660_gt1, d660_gt1, W660_gt1, slope660_gt1/in.pressure_depth_derivative, true, false, false);

        // 660 for gt, combined
        const int phase_transition_index_660_gt_combined = in.phase_transition_index - phase_transition_index_pyrolite + partial_index_660_gt_combined;
        //const int phase_transition_index_660 = 0;
        const double d660_gt_combined = transition_depths[phase_transition_index_660_gt_combined];
        const double T660_gt_combined = transition_temperatures[phase_transition_index_660_gt_combined];
        const double W660_gt_combined = transition_widths[phase_transition_index_660_gt_combined];
        const double slope660_gt_combined = transition_slopes[phase_transition_index_660_gt_combined];
        std::pair<bool, double> result660_gt_combined = compute_point_to_line(in, T660_gt_combined, d660_gt_combined, W660_gt_combined, slope660_gt_combined/in.pressure_depth_derivative, true, false, false);

        if (result410.first)
          {
            // 410 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_410)
              function_value += 0.5*(1.0 + std::tanh(result410.second/W410));
          }
        if (result520.first)
          {
            // 520 for pyrolite
            if (phase_transition_index_pyrolite == partial_index_520)
              function_value += 0.5*(1.0 + std::tanh(result520.second/W520));
          }
        if (result560.first)
          {
            // 560 for pyrolite, Gt -> CaPv + Gt
            if (phase_transition_index_pyrolite == partial_index_560)
              function_value += 0.5*(1.0 + std::tanh(result560.second/W560));
          }
        if (result660.first)
          {
            // 660 for pyrolite, rw -> brg + fp
            if (phase_transition_index_pyrolite == partial_index_660)
              function_value += 0.5*(1.0 + std::tanh(result660.second/W660));
          }
        if (result660_gt.first && in.temperature < T660_gt)
          {
            // 660 for pyrolite, gt -> il
            if (phase_transition_index_pyrolite == partial_index_660_gt)
              function_value += 0.5*(1.0 + std::tanh(result660_gt.second/W660_gt));
          }
        if (result660_gt1.first && in.temperature < T660_gt1)
          {
            // 660 for pyrolite, il -> brg
            if (phase_transition_index_pyrolite == partial_index_660_gt1)
              function_value += 0.5*(1.0 + std::tanh(result660_gt1.second/W660_gt1));
          }
        if (result660_gt_combined.first && in.temperature > T660_gt_combined)
          {
            // 660 for pyrolite combined, at higher temperature, gt -> brg
            if (phase_transition_index_pyrolite == partial_index_660_gt_combined)
              function_value += 0.5*(1.0 + std::tanh(result660_gt_combined.second/W660_gt_combined));
          }
        return function_value;
      }

      template <int dim>
      std::pair<bool, double>
      compute_point_to_line (const PhaseFunctionInputs1<dim> &in,
                             const double T, const double P, const double W, const double slope,
                             bool by_depth, bool is_negative, bool is_vertical)
      {
        // In this approach, we define a transition as a solid line and a range.
        // The solid line is a rigid boundary for the new phase.
        // While the range is a width of transition.
        double deviation;
        bool is_in;
        if (is_vertical)
          deviation = in.temperature - T;
        else
          {
            if (by_depth)
              deviation = in.depth - P - slope * (in.temperature - T);
            else
              deviation = in.pressure - P - slope * (in.temperature - T);
          }
        // In this approach, a transition must has a direction in defination
        // We need an opposite direction when the transition defined from higher pressure to lower pressure.
        if (is_negative)
          deviation *= -1.0;
        // Deviation must be smaller than 2*W, value for function would be 0.04 there.
        is_in = (deviation > -2.0 * W);
        return std::make_pair(is_in, deviation);
      }


      double average_deviation(double x1, double x2, double pinpoint)
      {
        double deviation = 0.0;
        if (x1 < pinpoint && x2 < pinpoint)
          {
            deviation = pinpoint - sqrt(pow(pinpoint - x1, 2.0) + pow(pinpoint - x2, 2.0));
          }
        else if (x1 < pinpoint && x2 > pinpoint)
          {
            deviation = x1;
          }
        else if (x1 > pinpoint && x2 < pinpoint)
          {
            deviation = x2;
          }
        else
          {
            deviation = std::min(x1, x2);
          }
        return deviation;
      }

      // RK4 solver
      // Constructor
      IRK4Solver::IRK4Solver()
      {
        double sqrt3 = std::sqrt(3.0) / 6.0;

        A =
        {
          {0.25, 0.25 - sqrt3},
          {0.25 + sqrt3, 0.25}
        };

        b = {0.5, 0.5};
        c = {0.5 - sqrt3, 0.5 + sqrt3};
      }

      // Solve method with debug option
      std::pair<std::vector<double>, std::vector<std::vector<double>>>
      IRK4Solver::solve(const std::function<std::vector<double>(double, const std::vector<double>&)> &f,
                        const std::vector<double> &y0, const std::pair<double, double> &t_span, double h, bool debug)
      {
        double t0 = t_span.first;
        double t_end = t_span.second;

        std::vector<double> t_values = {t0};
        std::vector<std::vector<double>> y_values = {y0};

        std::vector<double> y = y0;
        double t = t0;

        const double epsilon = 1e-6;

        while (t + h < t_end * (1.0 + epsilon))
          {

            //increase t
            t += h;

            size_t n = y.size();
            size_t stages = b.size();

            // Residual function for nonlinear solver
            auto residual = [&](const std::vector<std::vector<double>> &Y) -> std::vector<std::vector<double>>
            {
              std::vector<std::vector<double>> res(stages, std::vector<double>(n, 0.0));
              for (size_t i = 0; i < stages; ++i)
                {
                  for (size_t k = 0; k < n; ++k)
                    {
                      res[i][k] = Y[i][k] - y[k];
                      for (size_t j = 0; j < stages; ++j)
                        {
                          auto f_eval = f(t + c[j] * h, Y[j]);
                          res[i][k] -= h * A[i][j] * f_eval[k];
                        }
                    }
                }
              return res;
            };

            // Initial guess for stages
            std::vector<std::vector<double>> Y(stages, y);

            // Solve for Y using simple fixed-point iteration (replaceable with a better solver)
            size_t max_iterations = 100;
            double tolerance = 1e-6;
            for (size_t iter = 0; iter < max_iterations; ++iter)
              {
                auto res = residual(Y);

                // Debug: Print residuals, relative residuals, and temporary solutions if debug is true
                if (debug)
                  {
                    std::cout << "Iteration " << iter << ", t = " << t << ":\n";
                    for (size_t i = 0; i < stages; ++i)
                      {
                        std::cout << "  Stage " << i << " Residual: ";
                        for (size_t k = 0; k < n; ++k)
                          {
                            std::cout << res[i][k] << " ";
                          }
                        std::cout << "\n";

                        std::cout << "  Stage " << i << " Relative Residual: ";
                        for (size_t k = 0; k < n; ++k)
                          {
                            double relative_residual = std::abs(res[i][k]) / (std::abs(Y[i][k]) + 1e-12);
                            std::cout << relative_residual << " ";
                          }
                        std::cout << "\n";

                        std::cout << "  Stage " << i << " Solution: ";
                        for (size_t k = 0; k < n; ++k)
                          {
                            std::cout << Y[i][k] << " ";
                          }
                        std::cout << "\n";
                      }
                  }

                // Update Y values based on residuals
                bool converged = true;
                for (size_t i = 0; i < stages; ++i)
                  {
                    for (size_t k = 0; k < n; ++k)
                      {
                        double relative_residual = std::abs(res[i][k]) / (std::abs(Y[i][k]) + 1e-12); // Add small value to avoid division by zero
                        Y[i][k] -= res[i][k];
                        if (relative_residual > tolerance)
                          {
                            converged = false;
                          }
                      }
                  }

                if (converged) break;
                if (iter == max_iterations - 1)
                  {
                    throw std::runtime_error("Nonlinear solver failed to converge");
                  }
              }

            // Update solution using stages
            std::vector<double> dy(n, 0.0);
            for (size_t k = 0; k < n; ++k)
              {
                for (size_t i = 0; i < stages; ++i)
                  {
                    dy[k] += b[i] * f(t + c[i] * h, Y[i])[k];
                  }
              }

            for (size_t k = 0; k < n; ++k)
              {
                y[k] += h * dy[k];
              }

            t_values.push_back(t);
            y_values.push_back(y);

          }

        return {t_values, y_values};
      }

// Constructor
      PTKinetics::PTKinetics():
        A(std::exp(-18.0)), n(3.2),
        dS(7.7), dV(3.16e-6),
        dHa(274e3), Vstar(3.3e-6),
        fs(1e-3), Vm(4.05e-5),
        gamma(0.6), K0(3.65e38),
        nucleation_type(0)
      {}

      PTKinetics::PTKinetics(const double A_, const double n_, const double dS_, const double dV_,
                             const double dHa_, const double Vstar_, const double fs_, const double Vm_
                             , const double gamma_, const double K0_, const int nucleation_type)
        : A(A_), n(n_), dS(dS_), dV(dV_),
          dHa(dHa_), Vstar(Vstar_), fs(fs_), Vm(Vm_),
          gamma(gamma_), K0(K0_), nucleation_type(nucleation_type)
      {}

// Growth rate (without ΔGr term)
      double PTKinetics::growth_rate_P1(double P, double T, double Coh) const
      {
        return A * std::pow(Coh, n) * std::exp(-(dHa + P * Vstar) / (R * T));
      }

// Full growth rate (with ΔGr term)
      double PTKinetics::growth_rate(double P, double T, double Peq, double Teq, double Coh) const
      {
        assert(P > Peq && "P must be greater than Peq");
        double dGr = dV * (P - Peq);
        return growth_rate_P1(P, T, Coh) * T * (1.0 - std::exp(-dGr / (R * T)));
      }

// Nucleation rate
      double PTKinetics::nucleation_rate(double P, double T, double Peq, double Teq) const
      {
        assert(P > Peq && "P must be >= Peq");
        double dGr = dV * (P - Peq);
        double delta_G_hom = (16.0 * M_PI * fs * std::pow(Vm, 2) * std::pow(gamma, 3)) / (3.0 * std::pow(dGr, 2));
        double Q_a = dHa + P * Vstar;
        return K0 * T * std::exp(-delta_G_hom / (k * T)) * std::exp(-Q_a / (R * T));
      }

      // Function to calculate the Avrami number using corrected Equation (19)
      double calculate_avrami_number(double I_max, double Y_max, double kappa, double D)
      {
        // Compute the Avrami number
        return std::pow(D * D / kappa, 4) * I_max * std::pow(Y_max, 3);
      }

// Function to calculate dimensionless time (sigma_s)
      double calculate_sigma_s(double I_PT, double Y_PT, double d_0, double kappa, double D)
      {
        // Check for invalid inputs
        if (I_PT == 0 || Y_PT == 0 || d_0 == 0)
          {
            return std::numeric_limits<double>::infinity(); // Return infinity for invalid inputs
          }

        // Compute sigma_s using the given formula
        double sigma_s = (kappa / std::pow(D, 2)) * std::pow((I_PT * std::pow(Y_PT, 2) * d_0) / 6.7, -1.0 / 3.0);
        return sigma_s;
      }

// Function to calculate dimensionless time (sigma_s) for vector inputs
      double calculate_sigma_s(const std::vector<double> &I_array, const std::vector<double> &Y_array, double d_0, double kappa, double D)
      {
        // Validate input sizes
        if (I_array.size() != Y_array.size())
          {
            throw std::invalid_argument("I_array and Y_array must have the same size");
          }

        // Compute the average of I_array and Y_array
        double I_PT = std::accumulate(I_array.begin(), I_array.end(), 0.0) / I_array.size();
        double Y_PT = std::accumulate(Y_array.begin(), Y_array.end(), 0.0) / Y_array.size();

        // Call the scalar version of calculate_sigma_s
        return calculate_sigma_s(I_PT, Y_PT, d_0, kappa, D);
      }


// Function to solve for the extended volume after site saturation
      double solve_extended_volume_post_saturation(double Y, double s, double kappa, double D, double d0)
      {
        // Calculate the extended volume based on the given parameters
        double X3 = (6.7 * std::pow(D, 2) / (d0 * kappa)) * Y * s;
        return X3;
      }

// Overloaded version for vector inputs
      std::vector<double> solve_extended_volume_post_saturation(const double Y, const std::vector<double> &s, double kappa, double D, double d0)
      {
        // Validate input sizes

        // Compute extended volume for each pair of Y and s
        std::vector<double> X3(s.size());
        for (size_t i = 0; i < s.size(); ++i)
          {
            X3[i] = solve_extended_volume_post_saturation(Y, s[i], kappa, D, d0);
          }
        return X3;
      }

// Basic Ode and solvers
// Define the ODE system based on the modified Equation (18)
      std::vector<double> ode_system(double s, const std::vector<double> &X, double Av,
                                     const std::function<double(double)> &Y_prime,
                                     const std::function<double(double)> &I_prime)
      {
        // Extract state variables
        double X0 = X[0];
        double X1 = X[1];
        double X2 = X[2];
        double X3 = X[3];

        // Avrami factor
        double Av_factor = std::pow(Av, 0.25);

        // Calculate the derivatives based on the Avrami equation
        double dX3 = Av_factor * (4.0 * Y_prime(s) * X2);
        double dX2 = Av_factor * (M_PI * Y_prime(s) * X1);
        double dX1 = Av_factor * (2.0 * Y_prime(s) * X0);
        double dX0 = Av_factor * I_prime(s);

        // Return the derivatives as a vector
        return {dX0, dX1, dX2, dX3};
      }

// Solve the Modified Equation (18)
      std::vector<std::vector<double>> solve_modified_equations_eq18(
        double Av,
        const std::function<double(double)> &Y_prime_func,
        const std::function<double(double)> &I_prime_func,
        const std::pair<double, double> &s_span,
        const std::vector<double> &X_ini,
        int n_span,
        bool debug)
      {

        // Debugging output
        if (debug)
          {
            std::cout << "solve_modified_equations_eq18" << std::endl << "X0 = [";
            for (const auto &x : X_ini)
              {
                std::cout << x << " ";
              }
            std::cout << "]" << std::endl;

            std::cout << "s_span = (" << s_span.first << ", " << s_span.second << ")" << std::endl;
            std::cout << "Av = " << Av << std::endl;
            std::cout << "Y_prime_func(s_span.first) = " << Y_prime_func(s_span.first) << ", "
                      << "Y_prime_func(s_span.second) = " << Y_prime_func(s_span.second) << std::endl;
            std::cout << "I_prime_func(s_span.first) = " << I_prime_func(s_span.first) << ", "
                      << "I_prime_func(s_span.second) = " << I_prime_func(s_span.second) << std::endl;
          }

        // Define the ODE system for IRK4Solver
        auto odes = [&](double s, const std::vector<double> &X) -> std::vector<double>
        {
          return ode_system(s, X, Av, Y_prime_func, I_prime_func);
        };

        // Time step size
        double h = (s_span.second - s_span.first) / n_span;

        // Create the solver
        IRK4Solver solver;

        // Solve the ODE system
        auto [t_values, X_values] = solver.solve(odes, X_ini, s_span, h, debug);

        if (debug)
          {
            std::cout << "      h = " << h << std::endl;
            std::cout << "      X_values.size() = " << X_values.size() << std::endl;
            size_t _size = t_values.size();
            size_t _size1 = X_values[0].size();
            for (size_t i = 0; i < _size; ++i)
              {
                std::cout << "      t[" << i << "] = " << t_values[i] << std::endl;
                for (size_t j = 0; j < _size1; ++j)
                  {
                    std::cout << "      X["<< i << "][" << j << "] = " << X_values[i][j] << std::endl;
                  }
              }
          }

        // Safeguard for very small values
        const double threshold = 1e-12; // Define a small threshold
        for (auto &row : X_values)
          {
            for (auto &value : row)
              {
                if (std::abs(value) < threshold)
                  {
                    value = 0.0; // Set to zero to avoid computational issues
                  }
              }
          }

        // Combine results for return
        std::vector<std::vector<double>> solution(X_values.size(), std::vector<double>(X_ini.size(), 0.0));
        for (size_t i = 0; i < X_values.size(); ++i)
          {
            solution[i] = X_values[i];
          }

        return solution;
      }

      template <int dim>
      void
      MO_KINETICS<dim>::initialize()
      {
        PT_eq = {0.0, 0.0, 0.0}; // Initialize P, T, cl to 0.0
        n_col = 8;
      }

      template <int dim>
      void
      MO_KINETICS<dim>::declare_parameters (ParameterHandler &prm)
      {
        prm.declare_entry ("Grain size", "1e-2", Patterns::Double (0.),
                           "The fixed initial grain size of the material. "
                          );

        // exp(-18.0) = 1.52299797e-8
        prm.declare_entry ("Prefactor for boundary grain growth", "1.52299797e-8", Patterns::Double (0.),
                           "The value of the prefactor for grain boundary growth"
                          );

        prm.declare_entry ("OH exponent for boundary grain growth", "3.2", Patterns::Double (0.),
                           "The value of stress exponent number for grain boundary growth"
                          );

        prm.declare_entry ("Volume fraction change for metastable transition", "3.16e-6", Patterns::Double (0.),
                           "The value of volume fraction change for metastable transition"
                          );

        prm.declare_entry ("Activation enthalpy for boundary grain growth", "274e3", Patterns::Double (0.),
                           "The value of activation enthalpy for boundary grain growth"
                          );

        prm.declare_entry ("Activation volume for boundary grain growth", "3.3e-6", Patterns::Double (0.),
                           "The value of activation volume for boundary grain growth"
                          );

        prm.declare_entry ("Shape factor for gibbs affinity", "6e-4", Patterns::Double (0.),
                           "The value of shape factor for gibbs affinity"
                          );

        prm.declare_entry ("Molar volume of olivine", "4.05e-5", Patterns::Double (0.),
                           "The value of olivine molar volume"
                          );

        prm.declare_entry ("Surface energy for gibbs affinity", "0.46", Patterns::Double (0.),
                           "The value of Surface energy for gibbs affinity"
                          );

        prm.declare_entry ("Prefactor for nucleation", "1e30", Patterns::Double (0.),
                           "The value of prefactor for nucleation"
                          );

        prm.declare_entry ("Phase transition pressure at equilibrium", "13.5e9", Patterns::Double (0.),
                           "The value of pressure at equilibrium phase transition"
                          );

        prm.declare_entry ("Phase transition temperature at equilibrium", "1740.0", Patterns::Double (0.),
                           "The value of temperature at equilibrium phase transition"
                          );

        prm.declare_entry ("Phase transition Claypeyron slope at equilibrium", "2e6", Patterns::Double (0.),
                           "The value of Clayperon slope at equilibrium phase transition"
                          );

        prm.declare_entry ("Type of nucleation", "1", Patterns::Integer(0),
                           "The type of nucleation: 0 - volumetric; 1 - surface"
                          );
      }

      template <int dim>
      void
      MO_KINETICS<dim>::parse_parameters (ParameterHandler &prm)
      {
        d0 = prm.get_double("Grain size");
        A = prm.get_double("Prefactor for boundary grain growth");
        n = prm.get_double("OH exponent for boundary grain growth");
        dV = prm.get_double("Volume fraction change for metastable transition");
        dHa = prm.get_double("Activation enthalpy for boundary grain growth");
        Vstar = prm.get_double("Activation volume for boundary grain growth");
        fs = prm.get_double("Shape factor for gibbs affinity");
        Vm = prm.get_double("Molar volume of olivine");
        gamma = prm.get_double("Surface energy for gibbs affinity");
        K0 = prm.get_double("Prefactor for nucleation");
        nucleation_type = prm.get_integer("Type of nucleation");

        const double P0_eq = prm.get_double("Phase transition pressure at equilibrium");
        const double T0_eq = prm.get_double("Phase transition temperature at equilibrium");
        const double cl_eq = prm.get_double("Phase transition Claypeyron slope at equilibrium");
        PT_eq = {P0_eq, T0_eq, cl_eq};
      }

      template <int dim>
      void MO_KINETICS<dim>::linkAndSetKineticsModel()
      {
        kinetics = std::make_shared<PTKinetics>(
                     A, n,
                     dS, dV,             // dS, dv
                     dHa, Vstar,                // dHa, Vstar
                     fs, Vm,               //fs, Vm
                     gamma, K0,            // gamma, K0
                     nucleation_type // nucleation_type
                   );
        assert(kinetics->get_nucleation_type() == nucleation_type);
      }

// Fix the kinetics model
      template <int dim>
      void MO_KINETICS<dim>::setKineticsFixed(double P, double T, double Coh)
      {
        assert(!PT_eq.empty() && "PT_eq must be set before calling setKineticsFixed!");

        double P_eq = computeEqP(T); // compute equilibrium values
        double T_eq = computeEqT(P);

        // fix grain growth and nucleation functions to specific P, T condition
        // convert nucleation rate to volumetric nucleation rate
        if (kinetics->get_nucleation_type() == 0)
          {
            f_nu = 1.0;
          }
        else if (kinetics->get_nucleation_type() == 1)
          {
            f_nu = 6.7 / d0;
          }
        else
          {
            throw std::runtime_error("This nucleation type is not implemented yet.");
          }
        Y_func = [this, P, T, P_eq, T_eq, Coh](double t)
        {
          return kinetics->growth_rate(P, T, P_eq, T_eq, Coh);
        };
        I_func = [this, P, T, P_eq, T_eq](double t)
        {
          return f_nu*kinetics->nucleation_rate(P, T, P_eq, T_eq);
        };
      }

// Set phase transition equilibrium
      template <int dim>
      void MO_KINETICS<dim>::setPTEq(double P0, double T0, double cl)
      {
        PT_eq = {P0, T0, cl}; // Update vector with new values
      }

      template <int dim>
      double MO_KINETICS<dim>::computeEqP(double T)
      {
        return PT_eq[0] + PT_eq[2] * (T - PT_eq[1]);
      }

      template <int dim>
      double MO_KINETICS<dim>::computeEqT(double P)
      {
        return PT_eq[1] + (P - PT_eq[0]) / PT_eq[2];
      }

      template <int dim>
      double MO_KINETICS<dim>::growth_rate(double P, double T, double Coh)
      {
        double P_eq = computeEqP(T); // compute equilibrium values
        double T_eq = computeEqT(P);

        // return 0 if equilibrium phase transition condition is not met
        if (P > P_eq)
          {
            return kinetics->growth_rate(P, T, P_eq, T_eq, Coh);
          }
        else
          {
            return 0.0;
          }
      }

      template <int dim>
      double MO_KINETICS<dim>::nucleation_rate(double P, double T)
      {
        double P_eq = computeEqP(T); // compute equilibrium values
        double T_eq = computeEqT(P);

        // return 0 if equilibrium phase transition condition is not met
        if (P > P_eq)
          {
            return f_nu*kinetics->nucleation_rate(P, T, P_eq, T_eq);
          }
        else
          {
            return 0.0;
          }
      }

      template <int dim>
      std::pair<std::vector<std::vector<double>>, std::vector<bool>> MO_KINETICS<dim>::solveModifiedEquation(
        const std::pair<double, double> &t_span,
        const std::vector<double> &X_ini,
        bool is_saturated,
        int n_span,
        bool debug)
      {

        // Ensure previous steps are valid
        assert(Y_func && I_func && "Kinetics functions must be set before solving!");

        // Compute scaling variables
        double I_max = std::max(1e-50, I_func(t_span.first));
        double Y_max = std::max(1e-50, Y_func(t_span.first));
        double Av = calculate_avrami_number(I_max, Y_max, kappa, D);

        // Define non-dimensionalized versions of Y_func and I_func
        auto Y_prime_func = [this, Y_max](double s)
        {
          return Y_func(s * t_scale) / Y_max; // Scale time by t_scale
        };

        auto I_prime_func = [this, I_max](double s)
        {
          return I_func(s * t_scale) / I_max; // Scale time by t_scale
        };


        if (debug)
          {
            std::cout << "solveModifiedEquation: I_max = " << I_max << ", Y_max = " << Y_max << ", Av = " << Av << std::endl;
          }

        // Compute scaling factors
        X_scale_array =
        {
          std::pow(I_max, 0.75) *std::pow(Y_max, -0.75),
          std::pow(I_max, 0.5) *std::pow(Y_max, -0.5),
          std::pow(I_max, 0.25) *std::pow(Y_max, -0.25),
          1.0
        };

        // Non-dimensionalize the initial solution
        std::vector<double> X_ini_nd(4, 0.0);
        for (size_t k = 0; k < X_ini.size(); ++k)
          {
            X_ini_nd[k] = X_ini[k] / X_scale_array[k];
          }

        // Non-dimensionalize the time span
        std::pair<double, double> s_span = {t_span.first / t_scale, t_span.second / t_scale};
        std::vector<double> s_values(n_span+1);
        for (int i = 0; i <= n_span; ++i)
          {
            s_values[i] = s_span.first + i * (s_span.second - s_span.first) / n_span;
          }

        std::vector<double> I_array(n_span+1), Y_array(n_span+1);
        for (int i = 0; i <= n_span; ++i)
          {
            I_array[i] = I_func(s_values[i] * t_scale);
            Y_array[i] = Y_func(s_values[i] * t_scale);
          }

        // Compute saturation condition
        double s_saturation = calculate_sigma_s(I_array, Y_array, d0, kappa, D);

        if (debug)
          {
            std::cout << "I_array[0] = " << I_array[0] << std::endl;
            std::cout << "Y_array[0] = " << Y_array[0] << std::endl;
            std::cout << "solveModifiedEquation: s_saturation = " << s_saturation << "\n";
            std::cout << "solveModifiedEquation: t_saturation = " << s_saturation *t_scale << "\n";
            std::cout << "solveModifiedEquation: is_saturated = " << is_saturated << "\n";
          }

        // Initialize result containers
        std::vector<std::vector<double>> X_array(4, std::vector<double>(n_span+1, 0.0));
        std::vector<bool> is_saturated_array(n_span+1, false);

        if (!is_saturated)
          {
            auto it = std::find_if(s_values.begin(), s_values.end(), [s_saturation, s_values](double s)
            {
              return s > s_values[0] + s_saturation;
            });
            if (it != s_values.end() && std::distance(s_values.begin(), it) > 1)
              {
                // Pre-saturation & saturation
                int i0 = std::distance(s_values.begin(), it);
                std::pair<double, double> s_span_us = {s_values[0], s_values[i0]};
                auto solution_nd = solve_modified_equations_eq18(Av, Y_prime_func, I_prime_func, s_span_us, X_ini_nd, i0, debug);


                // Scale and store pre-saturation results
                for (size_t j = 0; j < solution_nd.size(); ++j)
                  {
                    for (size_t k = 0; k < solution_nd[j].size(); ++k)
                      {
                        X_array[k][j] = solution_nd[j][k] * X_scale_array[k];
                      }
                  }
                // Pre-assign values post-saturation
                for (size_t j = solution_nd.size(); j < static_cast<size_t>(n_span + 1); ++j)
                  {
                    for (size_t k = 0; k < solution_nd[solution_nd.size()-1].size(); ++k)
                      {
                        X_array[k][j] = solution_nd[solution_nd.size()-1][k] * X_scale_array[k];
                      }
                  }

                // Post-saturation, increment from values derived by the analytical solution
                auto post_saturation = solve_extended_volume_post_saturation(Y_max, s_values, kappa, D, d0);
                auto post_saturation_ini = solve_extended_volume_post_saturation(Y_max, s_values[i0], kappa, D, d0);
                for (int i = i0; i <= n_span; ++i)
                  {
                    X_array[3][i] = X_array[3][i0] + post_saturation[i] - post_saturation_ini;
                  }
                std::fill(is_saturated_array.begin() + i0, is_saturated_array.end(), true);
              }
            else if (it != s_values.end() && std::distance(s_values.begin(), it) <= 1)
              {
                // assign the initial values
                for (size_t k = 0; k < X_ini.size(); ++k)
                  {
                    X_array[k][0] = X_ini[k];
                  }

                // solve the nucleation in range of s_values[0], s_saturation
                std::pair<double, double> s_span_us = {s_values[0], s_values[0] + s_saturation};

                auto solution_nd = solve_modified_equations_eq18(Av, Y_prime_func, I_prime_func, s_span_us, X_ini_nd, n_span, debug);

                for (size_t j = 1; j < static_cast<size_t>(n_span + 1); ++j)
                  {
                    for (size_t k = 0; k < solution_nd[j].size(); ++k)
                      {
                        X_array[k][j] = solution_nd[solution_nd.size()-1][k] * X_scale_array[k];
                      }
                  }

                // saturation at the 1st sub-step
                // Post-saturation, increment from values derived by the analytical solution
                auto post_saturation = solve_extended_volume_post_saturation(Y_max, s_values, kappa, D, d0);
                auto post_saturation_ini = solve_extended_volume_post_saturation(Y_max, s_values[1], kappa, D, d0);
                for (size_t j = 1; j < static_cast<size_t>(n_span + 1); ++j)
                  {
                    X_array[3][j] = X_array[3][1] + post_saturation[j] - post_saturation_ini;
                  }

                std::fill(is_saturated_array.begin()+1, is_saturated_array.end(), true);
              }
            else
              {
                // solve the unsaturated condition
                auto solution_nd = solve_modified_equations_eq18(Av, Y_prime_func, I_prime_func, s_span, X_ini_nd, n_span, debug);

                // Debug: Print X_scale_array if debug is true
                if (debug)
                  {
                    std::cout << "X_scale_array: ";
                    for (const auto &scale : X_scale_array)
                      {
                        std::cout << scale << " ";
                      }
                    std::cout << "\n";
                  }

                // assign the values from the solution
                for (size_t j = 0; j < solution_nd.size(); ++j)
                  {
                    for (size_t k = 0; k < solution_nd[j].size(); ++k)
                      {
                        // Debug print for each value
                        if (debug)
                          {
                            std::cout << "solution_nd[" << j << "][" << k << "] = "  << solution_nd[j][k] << ", X_scale_array[" << k << "] = " << X_scale_array[k] << std::endl;
                          }
                        X_array[k][j] = solution_nd[j][k] * X_scale_array[k];
                      }
                  }
              }
          }
        else
          {
            // assign the initial values
            for (size_t j = 0; j < static_cast<size_t>(n_span + 1); ++j)
              {
                for (size_t k = 0; k < X_ini.size(); ++k)
                  {
                    X_array[k][j] = X_ini[k];
                  }
              }

            // Full saturation, increment from values derived by the analytical solution
            auto post_saturation = solve_extended_volume_post_saturation(Y_max, s_values, kappa, D, d0);
            auto post_saturation_ini = solve_extended_volume_post_saturation(Y_max, s_values[0], kappa, D, d0);
            for (size_t j = 0; j < static_cast<size_t>(n_span + 1); ++j)
              {
                X_array[3][j] = X_ini[3] + post_saturation[j] - post_saturation_ini;
              }
            std::fill(is_saturated_array.begin(), is_saturated_array.end(), true);
          }

        return {X_array, is_saturated_array};
      }

      template <int dim>
      std::vector<std::vector<double>> MO_KINETICS<dim>::solve(const double P, const double T, const double t_min, const double t_max,
                                                                const int n_t, const int n_span, const bool debug, const std::vector<double> X_ini, const bool is_saturated_ini)
      {

        // Initialize variables
        // Check if X has exactly 4 elements
        if (X_ini.size() != 4)
          {
            throw std::invalid_argument("Error: X must have exactly 4 elements.");
          }

        std::vector<std::vector<double>> results(n_t * n_span+1, std::vector<double>(n_col, 0.0));

        // Compute equilibrium pressure
        double Peq = computeEqP(T);

        // Set initial values
        double last_derivative = 0.0;
        std::vector<double> X(X_ini) ;
        bool is_saturated = is_saturated_ini;
        // Loop over time steps
        for (int i_t = 0; i_t < n_t; ++i_t)
          {
            if (debug)
              {
                std::cout << "i_t: " << i_t << std::endl;
              }

            // Define the time span for the current step
            double t_piece_min = t_min + (t_max - t_min) / n_t * i_t;
            double t_piece_max = t_min + (t_max - t_min) / n_t * (i_t + 1);
            std::pair<double, double> t_span = {t_piece_min, t_piece_max};

            std::vector<std::vector<double>> X_array(4, std::vector<double>(n_span+1, 0.0));
            std::vector<bool> is_saturated_array(n_span+1, false);

            if (P > Peq)
              {
                // Solve the kinetics if equilibrium condition is met
                auto solution = solveModifiedEquation(t_span, X, is_saturated, n_span, false);
                X_array = solution.first;
                is_saturated_array = solution.second;
                X = {X_array[0].back(), X_array[1].back(), X_array[2].back(), X_array[3].back()};
                is_saturated = is_saturated_array.back();
              }
            else
              {
                // Assign trivial values if equilibrium condition is not met
                X = {0.0, 0.0, 0.0, 0.0};
                is_saturated = false;
              }

            // Compute results for each step in the current time span
            std::vector<double> V_array(n_span+1);
            const double t_interval = (t_piece_max - t_piece_min) / n_span;
            for (int j = 0; j <= n_span; ++j)
              {
                // check all values are positive
                for (size_t i_x = 0; i_x < 4; ++i_x)
                  {
                    // local cpp
                    if (X_array[i_x][j] < 0.0)
                      {
                        const int foo0 = 0;
                        const int foo1 = 1;
                        /*
                        throw std::runtime_error(
                          "Error: X_array contains negative value at index " + std::to_string(i_x) +
                          " (value = " + std::to_string(X_array[i_x][j]) + ").\n"
                          "Context: P = " + std::to_string(P) +
                          ", Peq = " + std::to_string(Peq) +
                          ", T = " + std::to_string(T) +
                          ".\nAll values in X_array must be > 0.");
                        */
                      }
                    AssertThrow (X_array[i_x][j] >= 0.0,
                                 ExcMessage("X_array contains negative value at index " + std::to_string(i_x) + ", " + std::to_string(j) +
                                            " (value = " + std::to_string(X_array[i_x][j]) + ").\n"
                                            "Context: P = " + std::to_string(P) +
                                            ", Peq = " + std::to_string(Peq) +
                                            ", T = " + std::to_string(T) + ".\nAll values must be > 0."));
                  }

                double threshold = 50.0; // Define a threshold for large values
                if (X_array[3][j] > threshold)
                  {
                    // If X_array[3][j] is too large, directly set V_array[j] to 1
                    V_array[j] = 1.0;
                  }
                else
                  {
                    // Otherwise, compute the exponential term
                    V_array[j] = 1.0 - std::exp(-X_array[3][j]);
                  }

                // debug
                if (P > Peq)
                  {
                    const int foo0 = 0;
                    const int foo1 = 1;
                  }

                // note i_t * n_span + j would have the same value
                // for the last node the piece and the first node in the consequtive piece
                // compute the derivative and recode the last value to append the next "first" value
                const double derivative = (j==0)? last_derivative: (V_array[j] - V_array[j-1])/t_interval;
                if (j == n_span)
                  last_derivative = derivative;

                results[i_t * n_span + j] =
                {
                  t_piece_min + t_interval * j,
                  X_array[0][j],
                  X_array[1][j],
                  X_array[2][j],
                  X_array[3][j],
                  V_array[j],
                  static_cast<double>(is_saturated_array[j]),
                  derivative
                };
              }
          }
        return results;
      }


    }
  }
}


// explicit instantiations
namespace aspect
{
  namespace MaterialModel
  {
    namespace MaterialUtilities
    {
#define INSTANTIATE(dim) \
  template void fill_averaged_equation_of_state_outputs<dim> (const EquationOfStateOutputs<dim> &, \
                                                              const std::vector<double> &, \
                                                              const std::vector<double> &, \
                                                              const unsigned int, \
                                                              MaterialModelOutputs<dim> &); \
  template struct PhaseFunctionInputs<dim>; \
  template struct PhaseFunctionInputs1<dim>; \
  template class PhaseFunctionDiscrete<dim>; \
  template class PhaseFunction<dim>;\
  template class MO_KINETICS<dim>;

      ASPECT_INSTANTIATE(INSTANTIATE)

#undef INSTANTIATE
    }
  }
}

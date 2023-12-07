/*
  Copyright (C) 2011 - 2019 by the authors of the ASPECT code.

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

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/fe/fe_values.h>
#include <aspect/global.h>
#include <aspect/simulator_signals.h>
#include <aspect/utilities.h>
#include <aspect/geometry_model/interface.h>
              


namespace aspect
{
  using namespace dealii;

  // Global variables (to be set by parameters)
  bool prescribe_internal_temperatures;
  Utilities::Coordinates::CoordinateSystem coordinate_system_indicator_function;
  Utilities::Coordinates::CoordinateSystem coordinate_system_function;
  std::string temperature_model, geometry_model_name;
  // these are defined for the plate model
  double top_temperature, subducting_plate_velocity, 
  overiding_plate_age, x_extent, y_extent, area_width,
  potential_mantle_temperature, overiding_plate_area_width,
  theta_M, theta_M_rad, outer_radius;

  // Because we do not initially know what dimension we're in, we need
  // function parser objects for both 2d and 3d.
  Functions::ParsedFunction<2> prescribed_temperature_indicator_function_2d (1);
  Functions::ParsedFunction<3> prescribed_temperature_indicator_function_3d (1);
  Functions::ParsedFunction<2> prescribed_temperature_function_2d (1);
  Functions::ParsedFunction<3> prescribed_temperature_function_3d (1);

  double plate_model_T(const double c1, const double c2)
  {
      // Use plate model to compute the temperature, migrated from the world builder
      // in order to generate consistent result with the world builder.
      const double x = c1;
      const double y = c2;
      const double thermal_diffusivity = 1e-6;
      const double thermal_expansion_coefficient = 3e-5;
      const double gravity_norm = 10.0;
      const double specific_heat = 1250.0;
      const double max_depth = 150e3 ; // same with the wb file
      const int sommation_number = 100; // same as in the World Builder
      const double distance_ridge = x;
      const double depth = y_extent - y;
      const double age = distance_ridge / subducting_plate_velocity;
      const double bottom_temperature_local =  potential_mantle_temperature *
                                                std::exp(((thermal_expansion_coefficient* gravity_norm) /
                                                           specific_heat) * depth);
  
      double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);
  
      for (int i = 1; i<sommation_number+1; ++i)
        {
          // suming over the "sommation_number"
          // use a spreading ridge around the left corner and a constant age around the right corner 
          if (x < area_width)
            temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp((((subducting_plate_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                     std::sqrt(((subducting_plate_velocity*subducting_plate_velocity*max_depth*max_depth) /
                                                (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * M_PI * M_PI)) *
                                    ((subducting_plate_velocity * age) / max_depth)));
          else if (x > x_extent - area_width)
            temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp(-1.0 * i * i * M_PI * M_PI * thermal_diffusivity * overiding_plate_age / (max_depth * max_depth)));
  
        }
        return temperature;  
  }
  
  double plate_model_T1_2d(const Point<2> &p, const double depth)
  {
      // Use plate model to compute the temperature, migrated from the world builder
      // in order to generate consistent result with the world builder.
      // Fix the overiding plate side by creating a minor ridge there.
      const double x = p[0];
      const double thermal_diffusivity = 1e-6;
      const double thermal_expansion_coefficient = 3e-5;
      const double gravity_norm = 10.0;
      const double specific_heat = 1250.0;
      const double max_depth = 150e3 ; // same with the wb file
      const int sommation_number = 100; // same as in the World Builder
      const double overriding_plate_pseudo_velocity = overiding_plate_area_width / overiding_plate_age;
      const double bottom_temperature_local =  potential_mantle_temperature *
                                                std::exp(((thermal_expansion_coefficient* gravity_norm) /
                                                           specific_heat) * depth);
  
      double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);

      double age;
  
      for (int i = 1; i<sommation_number+1; ++i)
        {
          // suming over the "sommation_number"
          // use a spreading ridge around the left corner and a pseudo spreading ridge around the right corner 
          if (x < area_width)
          {
            // for the subducting plate
            age = x / subducting_plate_velocity;
            temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp((((subducting_plate_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                     std::sqrt(((subducting_plate_velocity*subducting_plate_velocity*max_depth*max_depth) /
                                                (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * M_PI * M_PI)) *
                                    ((subducting_plate_velocity * age) / max_depth)));
          }
          else if (x > x_extent - overiding_plate_area_width)
            {
              // for the overiding plate
              age = (x_extent - x) / overiding_plate_area_width * overiding_plate_age;
              temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp((((overriding_plate_pseudo_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                     std::sqrt(((overriding_plate_pseudo_velocity*overriding_plate_pseudo_velocity*max_depth*max_depth) /
                                                (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * M_PI * M_PI)) *
                                    ((overriding_plate_pseudo_velocity * age) / max_depth)));
            }
        }
        return temperature;  
  }
  
  double plate_model_T1_2d_chunk(const Point<2> &p, const double depth)
  {
      // Use plate model to compute the temperature, migrated from the world builder
      // in order to generate consistent result with the world builder.
      // Fix the overiding plate side by creating a minor ridge there.
      const double theta = std::atan2(p[1], p[0]);
      const double thermal_diffusivity = 1e-6;
      const double thermal_expansion_coefficient = 3e-5;
      const double gravity_norm = 10.0;
      const double specific_heat = 1250.0;
      const double max_depth = 150e3 ; // same with the wb file
      const int sommation_number = 100; // same as in the World Builder
      const double overriding_plate_pseudo_velocity = overiding_plate_area_width / overiding_plate_age;
      const double bottom_temperature_local =  potential_mantle_temperature *
                                                std::exp(((thermal_expansion_coefficient* gravity_norm) /
                                                           specific_heat) * depth);
  
      double temperature = top_temperature + (bottom_temperature_local - top_temperature) * (depth / max_depth);
  
      for (int i = 1; i<sommation_number+1; ++i)
        {
          // suming over the "sommation_number"
          // use a spreading ridge around the left corner and a pseudo spreading ridge around the right corner 
          if (outer_radius * theta < area_width)
          {
            const double dist = outer_radius * theta;
            // for the subducting plate
            const double age = dist / subducting_plate_velocity;
            temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp((((subducting_plate_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                     std::sqrt(((subducting_plate_velocity*subducting_plate_velocity*max_depth*max_depth) /
                                                (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * M_PI * M_PI)) *
                                    ((subducting_plate_velocity * age) / max_depth)));
          }
          else if (outer_radius * (theta_M_rad - theta) < overiding_plate_area_width)
            {
              const double dist = outer_radius * (theta_M_rad - theta);
              // for the overiding plate
              const double age = dist / overiding_plate_area_width * overiding_plate_age;
              temperature = temperature + (bottom_temperature_local - top_temperature) *
                          ((2 / (double(i) * M_PI)) * std::sin((double(i) * M_PI * depth) / max_depth) *
                           std::exp((((overriding_plate_pseudo_velocity * max_depth)/(2 * thermal_diffusivity)) -
                                     std::sqrt(((overriding_plate_pseudo_velocity*overriding_plate_pseudo_velocity*max_depth*max_depth) /
                                                (4*thermal_diffusivity*thermal_diffusivity)) + double(i) * double(i) * M_PI * M_PI)) *
                                    ((overriding_plate_pseudo_velocity * age) / max_depth)));
            }
        }
        return temperature;  
  }

  /**
   * Declare additional parameters.
   */
  void declare_parameters(const unsigned int dim,
                          ParameterHandler &prm)
  {
    prm.declare_entry ("Prescribe internal temperatures", "false",
                       Patterns::Bool (),
                       "Whether or not to use any prescribed internal temperatures. "
                       "Locations in which to prescribe temperatures are defined "
                       "in section ``Prescribed temperatures/Indicator function'' "
                       "and the temperatures are defined in section ``Prescribed "
                       "temperatures/Temperature function''. Indicators are evaluated "
                       "at the center of each cell, and all DOFs associated with "
                       "the specified temperature component at the indicated cells "
                       "are constrained."
                      );
    // Declare temperature & Geometry information for the plate model
    prm.declare_entry("Adiabatic surface temperature", "1673.0",
                     Patterns::Double (0.),
                     "Potential temperature"
                     );
    prm.enter_subsection ("Geometry");
    {
      prm.declare_entry ("Model name", "unspecified",
                        Patterns::Anything(),
                        "Geometry models");
      prm.enter_subsection ("Box");
      {
        prm.declare_entry("X extent", "6.7830e+06",
                          Patterns::Double (0.),
                          "X Extent of the Box"
                          );
        prm.declare_entry("Y extent", "2.8900e6",
                          Patterns::Double (0.),
                          "Y Extent of the Box"
                          );
      }
      prm.leave_subsection();

      prm.enter_subsection ("Chunk");
      {
          prm.declare_entry ("Chunk outer radius", "1.",
                             Patterns::Double (0.),
                             "Radius at the top surface of the chunk. Units: \\si{\\meter}.");
          
          prm.declare_entry ("Chunk maximum longitude", "0.",
                             Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                             "Minimum longitude of the chunk. Units: degrees.");
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();

    // Declare variables defined in the "Prescribed temperatures" subsection
    prm.enter_subsection ("Prescribed temperatures");
    {
      prm.enter_subsection ("Indicator function");
      {
        prm.declare_entry ("Coordinate system", "cartesian",
                           Patterns::Selection ("cartesian|spherical|depth"),
                           "A selection that determines the assumed coordinate "
                           "system for the function variables. Allowed values "
                           "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                           "are interpreted as r,phi or r,phi,theta in 2D/3D "
                           "respectively with theta being the polar angle. `depth' "
                           "will create a function, in which only the first "
                           "parameter is non-zero, which is interpreted to "
                           "be the depth of the point.");
        if (dim == 2)
          Functions::ParsedFunction<2>::declare_parameters (prm, 1);
        else
          Functions::ParsedFunction<3>::declare_parameters (prm, 1);
      }
      prm.leave_subsection ();

      prm.declare_entry("Model name", "function", 
                        Patterns::Selection ("function|plate model|plate model 1"),
                        "A selection that determines the model to use "
                        "for the indicated area."
                        );
      // variables used in the plate model
      prm.enter_subsection ("Plate model");
      {
        prm.declare_entry("Subducting plate velocity", "1.5855e-09",
                          Patterns::Double (0.),
                          "Velocity of the subducting plate"
                          );
        prm.declare_entry("Overiding plate age", "1.2614e+15",
                          Patterns::Double (0.),
                          "Velocity of the overiding plate"
                          );
        prm.declare_entry("Top temperature", "273.0",
                          Patterns::Double (0.),
                          "Temperature of the top boundary"
                          );
        prm.declare_entry("Area width", "2.75e5",
                          Patterns::Double (0.),
                          "Width of the area"
                          );
      }
      prm.leave_subsection();

      prm.enter_subsection ("Plate model 1");
      {
        prm.declare_entry("Subducting plate velocity", "1.5855e-09",
                          Patterns::Double (0.),
                          "Velocity of the subducting plate"
                          );
        prm.declare_entry("Overiding plate age", "1.2614e+15",
                          Patterns::Double (0.),
                          "Velocity of the overiding plate"
                          );
        prm.declare_entry("Top temperature", "273.0",
                          Patterns::Double (0.),
                          "Temperature of the top boundary"
                          );
        prm.declare_entry("Area width", "2.75e5",
                          Patterns::Double (0.),
                          "Width of the area"
                          );
        prm.declare_entry("Overiding area width", "2.75e5",
                          Patterns::Double (0.),
                          "Width of the area below the overiding plate."
                          );
      }
      prm.leave_subsection();

      prm.enter_subsection ("Temperature function");
      {
        prm.declare_entry ("Coordinate system", "cartesian",
                           Patterns::Selection ("cartesian|spherical|depth"),
                           "A selection that determines the assumed coordinate "
                           "system for the function variables. Allowed values "
                           "are `cartesian', `spherical', and `depth'. `spherical' coordinates "
                           "are interpreted as r,phi or r,phi,theta in 2D/3D "
                           "respectively with theta being the polar angle. `depth' "
                           "will create a function, in which only the first "
                           "parameter is non-zero, which is interpreted to "
                           "be the depth of the point.");
        if (dim == 2)
          Functions::ParsedFunction<2>::declare_parameters (prm, 1);
        else
          Functions::ParsedFunction<3>::declare_parameters (prm, 1);
      }
      prm.leave_subsection ();
    }
    prm.leave_subsection ();
  }

  template <int dim>
  void parse_parameters(const Parameters<dim> &,
                        ParameterHandler &prm)
  {
    prescribe_internal_temperatures = prm.get_bool ("Prescribe internal temperatures");
    // Get temperature & Geometry information for the plate model
    potential_mantle_temperature = prm.get_double("Adiabatic surface temperature");
    prm.enter_subsection ("Geometry model");
    {
      geometry_model_name = prm.get ("Model name");
      prm.enter_subsection ("Box");
      {
        x_extent = prm.get_double("X extent");
        y_extent = prm.get_double("Y extent");
      }
      prm.leave_subsection();

      prm.enter_subsection ("Chunk");
      {
          outer_radius = prm.get_double("Chunk outer radius");
          theta_M = prm.get_double("Chunk maximum longitude");
          theta_M_rad = theta_M * M_PI / 180.0;
      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
    // Get parameters defined for the Prescribed temperatures plugin
    prm.enter_subsection ("Prescribed temperatures");
    {
      prm.enter_subsection("Indicator function");
      {
        coordinate_system_indicator_function = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        try
          {
            if (dim == 2)
              prescribed_temperature_indicator_function_2d.parse_parameters (prm);
            else
              prescribed_temperature_indicator_function_3d.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed temperatures.Indicator function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
      }
      prm.leave_subsection();

      temperature_model = prm.get("Model name");
      if (temperature_model == "plate model" || temperature_model == "plate model 1")
        AssertThrow(dim == 2,
                    ExcMessage("Plate model only works in 2d right now"));

      if (temperature_model == "plate model")
      {
          // plate model
          prm.enter_subsection ("Plate model");
          {
            area_width = prm.get_double("Area width");
            subducting_plate_velocity = prm.get_double("Subducting plate velocity");
            overiding_plate_age = prm.get_double("Overiding plate age");
            top_temperature = prm.get_double("Top temperature");
          }
          prm.leave_subsection();
      }
      else if (temperature_model == "plate model 1")
      {
          // plate model 1
          prm.enter_subsection ("Plate model 1");
          {
            area_width = prm.get_double("Area width");
            subducting_plate_velocity = prm.get_double("Subducting plate velocity");
            overiding_plate_age = prm.get_double("Overiding plate age");
            overiding_plate_area_width = prm.get_double("Overiding area width");
            top_temperature = prm.get_double("Top temperature");
          }
          prm.leave_subsection();
      }

      prm.enter_subsection("Temperature function");
      {
        coordinate_system_function = Utilities::Coordinates::string_to_coordinate_system(prm.get("Coordinate system"));
        try
          {
            if (dim == 2)
              prescribed_temperature_function_2d.parse_parameters (prm);
            else
              prescribed_temperature_function_3d.parse_parameters (prm);
          }
        catch (...)
          {
            std::cerr << "ERROR: FunctionParser failed to parse\n"
                      << "\t'Prescribed temperatures.Temperature function'\n"
                      << "with expression\n"
                      << "\t'" << prm.get("Function expression") << "'";
            throw;
          }
      }
      prm.leave_subsection();
    }
    prm.leave_subsection ();
  }

  /**
   * This function retrieves the unit support points (in the unit cell) for the current element.
   * The DGP element used when 'set Use locally conservative discretization = true' does not
   * have support points. If these elements are in use, a fictitious support point at the cell
   * center is returned for each shape function that corresponds to the pressure variable,
   * whereas the support points for the temperature are correct; the fictitious points don't matter
   * because we only use this function when interpolating the temperature variable, and ignore the
   * evaluation at the pressure support points.
   */
  template <int dim>
  std::vector< Point<dim> >
  get_unit_support_points_for_temperature(const SimulatorAccess<dim> &simulator_access)
  {
    std::vector< Point<dim> > unit_support_points;
    if ( !simulator_access.get_parameters().use_locally_conservative_discretization )
      {
        return simulator_access.get_fe().get_unit_support_points();
      }
    else
      {
        //special case for discontinuous pressure elements, which lack unit support points
        std::vector< Point<dim> > unit_support_points;
        const unsigned int dofs_per_cell = simulator_access.get_fe().dofs_per_cell;
        for (unsigned int dof=0; dof < dofs_per_cell; ++dof)
          {
            // base will hold element, base_index holds node/shape function within that element
            const unsigned int
            base       = simulator_access.get_fe().system_to_base_index(dof).first.first,
            base_index = simulator_access.get_fe().system_to_base_index(dof).second;
            // get the unit support points for the relevant element
            std::vector< Point<dim> > my_support_points = simulator_access.get_fe().base_element(base).get_unit_support_points();
            if ( my_support_points.size() == 0 )
              {
                //manufacture a support point, arbitrarily at cell center
                if ( dim==2 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5) );
                if ( dim==3 )
                  unit_support_points.push_back( Point<dim> (0.5,0.5,0.5) );
              }
            else
              {
                unit_support_points.push_back(my_support_points[base_index]);
              }
          }
        return unit_support_points;
      }
  }

  namespace
  {
    template <int dim>
    const Point<2> &as_2d(const Point<dim> &p);

    template <>
    const Point<2> &as_2d(const Point<2> &p)
    {
      return p;
    }

    template <int dim>
    const Point<3> &as_3d(const Point<dim> &p);

    template <>
    const Point<3> &as_3d(const Point<3> &p)
    {
      return p;
    }

  }


  /**
   * This function is called by a signal which is triggered after the other constraints
   * have been calculated. This enables us to define additional constraints in the mass
   * matrix on any arbitrary degree of freedom in the model space.
   */
  template <int dim>
  void constrain_internal_temperatures (const SimulatorAccess<dim> &simulator_access,
                                      AffineConstraints<double> &current_constraints)
  {
    if (prescribe_internal_temperatures)
      {
        const std::vector< Point<dim> > points = get_unit_support_points_for_temperature(simulator_access);
        const Quadrature<dim> quadrature (points);
        FEValues<dim> fe_values (simulator_access.get_fe(), quadrature, update_quadrature_points);
        typename DoFHandler<dim>::active_cell_iterator cell;

        // Loop over all cells
        for (cell = simulator_access.get_dof_handler().begin_active();
             cell != simulator_access.get_dof_handler().end();
             ++cell)
          if (! cell->is_artificial())
            {
              fe_values.reinit (cell);
              std::vector<unsigned int> local_dof_indices(simulator_access.get_fe().dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);

              for (unsigned int q=0; q<quadrature.size(); q++)
                // If it's okay to constrain this DOF
                if (current_constraints.can_store_line(local_dof_indices[q]) &&
                    !current_constraints.is_constrained(local_dof_indices[q]))
                  {
                    // Get the temperature component index
                    const unsigned int c_idx =
                      simulator_access.get_fe().system_to_component_index(q).first;

                    // If we're on one of the temperature DOFs
                    if (c_idx == simulator_access.introspection().component_indices.temperature)
                      {
                        // we get time passed as seconds (always) but may want
                        // to reinterpret it in years
                        double time = simulator_access.get_time();
                        if (simulator_access.convert_output_to_years())
                          time /= year_in_seconds;

                        prescribed_temperature_indicator_function_2d.set_time (time);
                        prescribed_temperature_indicator_function_3d.set_time (time);
                        prescribed_temperature_function_2d.set_time (time);
                        prescribed_temperature_function_3d.set_time (time);

                        const Point<dim> p = fe_values.quadrature_point(q);

                        // Because we defined and parsed our parameter
                        // file differently for 2d and 3d we need to
                        // be sure to query the correct object for
                        // function values. The function parser
                        // objects expect points of a certain
                        // dimension, but Point p will be compiled for
                        // both 2d and 3d, so we need to do some trickery
                        // to make this compile.
                        double indicator, u_i;
                        if (dim == 2)
                          {
                            Utilities::NaturalCoordinate<dim> indicator_function_point =
                            simulator_access.get_geometry_model().cartesian_to_other_coordinates(p, coordinate_system_indicator_function);
                            Utilities::NaturalCoordinate<dim> function_point =
                            simulator_access.get_geometry_model().cartesian_to_other_coordinates(p, coordinate_system_function);
                            indicator = prescribed_temperature_indicator_function_2d.value
                                        (as_2d(Utilities::convert_array_to_point<dim>(indicator_function_point.get_coordinates())),
                                         0);
                            if (temperature_model == "function"){
                              // use the defined function
                              u_i       = prescribed_temperature_function_2d.value
                                          (as_2d(Utilities::convert_array_to_point<dim>(function_point.get_coordinates())),
                                           0);
                            }
                            else if (temperature_model == "plate model"){
                              // use the plate model
                              const Point<2> p = as_2d(Utilities::convert_array_to_point<dim>(function_point.get_coordinates()));
                              u_i       = plate_model_T(p[0], p[1]);
                            }
                            else if (temperature_model == "plate model 1"){
                              // use the plate model1
                              const Point<2> p2d = as_2d(Utilities::convert_array_to_point<dim>(function_point.get_coordinates()));
                              const double depth = simulator_access.get_geometry_model().depth(p);
                              if (indicator > 0.5)
                              {
                                if (geometry_model_name == "box")
                                  u_i       = plate_model_T1_2d(p2d, depth);
                                else if (geometry_model_name == "chunk") 
                                  u_i       = plate_model_T1_2d_chunk(p2d, depth);
                                else
                                {
                                  AssertThrow (false, ExcMessage ("The prescribed_temperature plugin only works"
                                                  "with the geometry of box or chunk:\n\n"));
                                }
                              }
                            }
                          }
                        else
                          {
                            Utilities::NaturalCoordinate<dim> indicator_function_point =
                            simulator_access.get_geometry_model().cartesian_to_other_coordinates(p, coordinate_system_indicator_function);
                            Utilities::NaturalCoordinate<dim> function_point =
                            simulator_access.get_geometry_model().cartesian_to_other_coordinates(p, coordinate_system_function);
                            indicator = prescribed_temperature_indicator_function_3d.value
                                        (as_3d(Utilities::convert_array_to_point<dim>(indicator_function_point.get_coordinates())),
                                         0);
                            u_i       = prescribed_temperature_function_3d.value
                                        (as_3d(Utilities::convert_array_to_point<dim>(function_point.get_coordinates())),
                                         0);
                          }

                        if (indicator > 0.5)
                          {
                            // Add a constraint of the form dof[q] = u_i
                            // to the list of constraints.
                            current_constraints.add_line (local_dof_indices[q]);
                            current_constraints.set_inhomogeneity (local_dof_indices[q], u_i);
                          }
                      }
                  }
            }
      }
  }

  // Connect declare_parameters and parse_parameters to appropriate signals.
  void parameter_connector ()
  {
    SimulatorSignals<2>::declare_additional_parameters.connect (&declare_parameters);
    SimulatorSignals<3>::declare_additional_parameters.connect (&declare_parameters);

    SimulatorSignals<2>::parse_additional_parameters.connect (&parse_parameters<2>);
    SimulatorSignals<3>::parse_additional_parameters.connect (&parse_parameters<3>);
  }

  // Connect constraints function to correct signal.
  template <int dim>
  void signal_connector (SimulatorSignals<dim> &signals)
  {
    signals.post_constraints_creation.connect (&constrain_internal_temperatures<dim>);
  }

  // Tell Aspect to send signals to the connector functions
  ASPECT_REGISTER_SIGNALS_PARAMETER_CONNECTOR(parameter_connector)
  ASPECT_REGISTER_SIGNALS_CONNECTOR(signal_connector<2>, signal_connector<3>)
}

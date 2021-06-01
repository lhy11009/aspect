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


#include <aspect/global.h>
#include "subduction_temperature2d_ellipse.h"
#include "spline.h"
#include <aspect/material_model/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/boundary_temperature/interface.h>
#include <aspect/geometry_model/box.h>
#include <aspect/geometry_model/two_merged_boxes.h>
#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <cmath>
#include <algorithm>
#include <functional>

#include <boost/lexical_cast.hpp>
#include <boost/math/tools/minima.hpp>

namespace aspect
{
  namespace InitialTemperature
  {

    template <int dim>
    Subduction2T<dim>::Subduction2T ()
    {}

    template <int dim>
    void
    Subduction2T<dim>::initialize ()
    {
      Utilities::AsciiDataInitial<dim>::initialize(1);
    }

// set up initial temperature
    template <int dim>
    double
    Subduction2T<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      const double SM=2.700e+05;
      const double Ro=6.371e6;
      std::vector<double> phases_below_lith_depth(3, 0.0);
      std::vector<double> phases_below_lith_temperature(3, 0.0);
      phases_below_lith_depth[0] = 410.0;
      phases_below_lith_depth[1] = 520.0;
      phases_below_lith_depth[2] = 670.0;
      phases_below_lith_temperature[0] = 1718.94072224;
      phases_below_lith_temperature[1] = 1746.82723665;
      phases_below_lith_temperature[2] = 1782.31922429;
      const std::array<double,dim> spherical_coordinates =
            aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      const double r = spherical_coordinates[0];
      const double phi = spherical_coordinates[1];
      double temperature;
      if(Ro-r < SM)
      {
        temperature = slab_model(r, phi);
        // for extended boussinesq, I use this along with "adiabatic" temperature profile
        if (extended_boussinesq == true)
          temperature = std::max(0.0001, temperature - adiabatic_surface_temperature_gradient * (Ro - r));
      }
      else
      {
        // temperature = katrina_mantle_temperature(r);
        // temperature = 0.0;
        temperature = adiabatic_surface_temperature;
      }
      if (extended_boussinesq == false)
        temperature = temperature_substract_adiabatic(r, temperature);
      return temperature;
    }

    /**
     * consist model temperature with mantle phase transition
    */
    /*
    template <int dim>
    double
    Subduction2T<dim>::
    below_lith_model (const double r, 
                      const double mantle_temperature, 
                      const double lith_depth,
                      const std::vector<double>& phases_below_lith_depth,
                      const std::vector<double>& phases_below_lith_temperature) const
    {
    }
    */
    
    template <int dim>
    double
    Subduction2T<dim>::
    slab_model (const double r, const double phi) const
    {
      const double Ro=6.371e6;
      const double PHIC=0.628319;
      const double PHIM=1.06465;
      const double AC=2.000e+05;
      const double BC=4.000e+05;
      const double ST=2.000e+05;
      const double SM=2.700e+05;
      const double VSUB=1.58549e-09;
      const double AGEOP=1.26144e15;
      const double TS=2.730e+02;
      const double TM=1.673e+03; 
      const double K=1.000e-06;

      const double tolerance = 1e-6;  // tolerance for the algorithm looking for shortest distance

      double temperature;


      // a thickness of thermal boundary on the upper slab boundary
      const double Ht_slab_in = thermal_boundary_width_factor_in * sqrt(K*BC/VSUB);
      const double Ht_slab_out = thermal_boundary_width_factor_out * sqrt(K*BC/VSUB);

      // a condition to determine if the point is within the ellipse envelop of
      // the initial slab
      double m = ellipse_equation_right(AC, BC, Ro*(phi-PHIC), r+BC-Ro);
      bool is_in_ellipse = (m <= 1);

      if (phi<=0.0)
      {
        // ghost domain
        temperature = TM;
      }
      else if ((phi<=PHIC)&&(phi>0))
      {
        // subducting plate
        temperature = half_space_cooling(TM, TS, abs(Ro-r), K, Ro*phi/VSUB);
      }
      else if (((Ro-r)<=ST)
                && is_in_ellipse
                &&(phi>PHIC))
                {
        // slab
        // temperature of the plate at the same depth
        const double plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, Ro*phi/VSUB);
        
        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, BC/AC, 
                                                                        Ro*(phi-PHIC)/AC, (r+BC-Ro)/AC, tolerance); // distance from focus to slab surface
        const double depth_in_slab = sqrt(ellipse_distance_sqr(AC, BC, Ro*(phi-PHIC), 
                                                               r+BC-Ro, results.first));
        const double slab_temperature = half_space_cooling(TM, TS, depth_in_slab, K, Ro*phi/VSUB);

        // cooling temperature as the smaller of these two

        const double cooling_temperature = std::min(plate_temperature, slab_temperature);

        // a perturbation on the upper surface
        // This is implemented as tranform from internal temperature of the slab to surface temperature.
        // The surface temperature is taken as the average of TS and temperature of overiding plate.
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        // method for perturbation
        double perturbation_fucntion_value = sin(M_PI/2.0*depth_in_slab/Ht_slab_in);
        const double perturbation_temperature = (depth_in_slab < Ht_slab_in)?
                                          (over_plate_temperature -cooling_temperature)
                                          * (1- perturbation_fucntion_value) / 2.0:
                                          0.0;
        
        // Final temperature is the sum of the two parts.
        temperature = cooling_temperature + perturbation_temperature;
                }
      else if (((Ro-r)>ST)
               &&((Ro-r)<=SM)
               && is_in_ellipse)
               {
        // temperature below the slab tip
        // A half-space cooling temperature is assigned once more but
        // an intermediate temperature is taken for the surface temperature.
        // depth of this interface
        // This have the same expression to depth_in_slab, except that this is actually below the initial slab

        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, BC/AC, 
                                                                        Ro*(phi-PHIC)/AC, (r+BC-Ro)/AC, tolerance); // distance from focus to slab surface
        const double depth_in_slab = sqrt(ellipse_distance_sqr(AC, BC, Ro*(phi-PHIC), 
                                                               r+BC-Ro, results.first));

        // Intermediate temperature on this interface
        const double T_interface = TS+(TM-TS)*(abs(Ro-r)-ST)/(SM-ST);

        const double cooling_temperature = half_space_cooling(TM, T_interface, depth_in_slab, K, Ro*phi/VSUB);
        
        // temperature overiding plate
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        
        // a perturbation on the upper surface
        // a thickness of thermal boundary on top
        // This is implemented as tranform from internal temperature of the slab to surface temperature.
        // The surface temperature is taken as the average of the intermediate temperature,
        // and temperature of overiding plate.
        const double perturbation_temperature = (depth_in_slab < Ht_slab_in)?
                                          (over_plate_temperature -cooling_temperature)
                                          * (1-sin(M_PI/2.0*depth_in_slab/Ht_slab_in)) / 2.0:
                                          0.0;

        // Final temperature is the sum of the two parts.
        temperature = cooling_temperature + perturbation_temperature;
               }
      else if (PHIM-phi>0)
      {
        // overiding plate
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        
        // This formula includes absolute value just to make sure it is consistent with the previous 'depth_in_slab'
        // the depth out of the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, BC/AC, 
                                                                        Ro*(phi-PHIC)/AC, (r+BC-Ro)/AC, tolerance); // distance from focus to slab surface
        const double depth_out_slab = sqrt(ellipse_distance_sqr(AC, BC, Ro*(phi-PHIC), 
                                                               r+BC-Ro, results.first));
        
        // a perturbation on the upper surface of the slab
        // This is implemented as tranform from temperature of the overiding plate to surface temperature.
        // Out of the upper boundary layer of the slab, this is 0.0
        double perturbation_temperature = 0.0;
        if (depth_out_slab<Ht_slab_out)
        { 
          if ((Ro-r)<ST)
          {
            // Surface temperature is taken as the average of TS and temperature of overiding plate.
            perturbation_temperature = (TS - over_plate_temperature)
                                        * (1-sin(M_PI/2.0*depth_out_slab/Ht_slab_out)) / 2.0;
          }
          else if ((Ro-r)<SM)
          {
            // Surface temperature is taken as the average of an intermediate temperature
            // and temperature of overiding plate.
            // Intermediate temperature on this interface
            const double T_interface = TS+(TM-TS)*(abs(Ro-r)-ST)/(SM-ST);

            perturbation_temperature = (T_interface - over_plate_temperature)
                                        * (1-sin(M_PI/2.0*depth_out_slab/Ht_slab_out)) / 2.0;
          }
        }
        temperature = over_plate_temperature + perturbation_temperature;
      }
      else
      {
        // mantle
        temperature = TM;
      }
      return temperature;
    }

    template <int dim>
    double
    Subduction2T<dim>::
    katrina_mantle_temperature(const double r) const
    {
      const double Ro = 6371e3;
      const double Dlith = 270e3;
      const double Dcmb = 2890e3;
      std::vector<double> phases_below_lith_depth(5, 0.0);
      std::vector<double> phases_below_lith_temperature(5, 0.0);
      std::vector<double> quary_depths(101, 0.0);
      phases_below_lith_depth[0] = Dlith;
      phases_below_lith_depth[1] = 410.0e3;
      phases_below_lith_depth[2] = 520.0e3;
      phases_below_lith_depth[3] = 670.0e3;
      phases_below_lith_depth[4] = Dcmb;
      phases_below_lith_temperature[0] = 1673.0;
      phases_below_lith_temperature[1] = 1718.94072224;
      phases_below_lith_temperature[2] = 1746.82723665;
      phases_below_lith_temperature[3] = 1782.31922429;
      phases_below_lith_temperature[4] = 2347.66;
      tk::spline s;
      s.set_points(phases_below_lith_depth, phases_below_lith_temperature);
      double temperature = s(Ro-r, true);
      return temperature;
    } 

    template <int dim>
    double
    Subduction2T<dim>::
    temperature_substract_adiabatic(const double r, const double temperature) const
    {
      const double Ro = 6371e3;
      const double Dadiabatic = 190e3;
      const double adiabatic_gradient = 0.25/1e3;
      if(r<Ro-Dadiabatic)
        return temperature - adiabatic_gradient * (Ro-r - Dadiabatic);
      else
        return temperature;
    }

    template <int dim>
    void
    Subduction2T<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry ("Adiabatic surface temperature", "0.",
                         Patterns::Double(),
                         "In order to make the problem in the first time step easier to "
                         "solve, we need a reasonable guess for the temperature and pressure. "
                         "To obtain it, we use an adiabatic pressure and temperature field. "
                         "This parameter describes what the `adiabatic' temperature would "
                         "be at the surface of the domain (i.e. at depth zero). Note "
                         "that this value need not coincide with the boundary condition "
                         "posed at this point. Rather, the boundary condition may differ "
                         "significantly from the adiabatic value, and then typically "
                         "induce a thermal boundary layer."
                         "\n\n"
                         "For more information, see the section in the manual that discusses "
                         "the general mathematical model.");
      prm.enter_subsection ("Initial temperature model");
      { 
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d_Vs_YT16.txt");
        prm.enter_subsection("Subduction 2d temperature ellipse");
        {
         // TODO:
         prm.declare_entry ("Thermal boundary width factor in", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
         prm.declare_entry ("Thermal boundary width factor out", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
         prm.declare_entry ("Adiabatic surface temperature gradient", "0.0",
                            Patterns::Double (),
                            "This controls the surface temperature gradient");
         prm.declare_entry("Extended boussinesq", "false", Patterns::Bool(),
                            "Initial temperature for extended boussinesq approximation");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Subduction2T<dim>::parse_parameters (ParameterHandler &prm)
    {
      adiabatic_surface_temperature   = prm.get_double ("Adiabatic surface temperature");
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Subduction 2d temperature ellipse");
        {
          // todo
          thermal_boundary_width_factor_in = prm.get_double("Thermal boundary width factor in");
          thermal_boundary_width_factor_out = prm.get_double("Thermal boundary width factor out");
          adiabatic_surface_temperature_gradient = prm.get_double("Adiabatic surface temperature gradient");
          extended_boussinesq = prm.get_bool("Extended boussinesq");
        }
        prm.leave_subsection();

        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
    
    template <int dim>
    double
    Adiabatic1<dim>::
    initial_temperature (const Point<dim> &position) const
    {
      // convert input ages to seconds
      const double age_top =    (this->convert_output_to_years() ? age_top_boundary_layer * year_in_seconds
                                 : age_top_boundary_layer);
      const double age_bottom = (this->convert_output_to_years() ? age_bottom_boundary_layer * year_in_seconds
                                 : age_bottom_boundary_layer);

      // First, get the temperature of the adiabatic profile at a representative
      // point at the top and bottom boundary of the model
      // if adiabatic heating is switched off, assume a constant profile
      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const Point<dim> bottom_point = this->get_geometry_model().representative_point(this->get_geometry_model().maximal_depth());
      const double adiabatic_surface_temperature = this->get_adiabatic_conditions().temperature(surface_point);
      const double adiabatic_bottom_temperature = (this->include_adiabatic_heating())
                                                  ?
                                                  this->get_adiabatic_conditions().temperature(bottom_point)
                                                  :
                                                  adiabatic_surface_temperature;

      // then, get the temperature at the top and bottom boundary of the model
      // if no boundary temperature is prescribed simply use the adiabatic.
      // This implementation assumes that the top and bottom boundaries have
      // prescribed temperatures and minimal_temperature() returns the value
      // at the surface and maximal_temperature() the value at the bottom.
      const double T_surface = (this->has_boundary_temperature()
                                ?
                                this->get_boundary_temperature_manager().minimal_temperature(
                                  this->get_fixed_temperature_boundary_indicators())
                                :
                                adiabatic_surface_temperature);
      const double T_bottom = (this->has_boundary_temperature()
                               ?
                               this->get_boundary_temperature_manager().maximal_temperature(
                                 this->get_fixed_temperature_boundary_indicators())
                               :
                               adiabatic_bottom_temperature);

      // get a representative profile of the compositional fields as an input
      // for the material model
      const double depth = this->get_geometry_model().depth(position);

      // look up material properties
      MaterialModel::MaterialModelInputs<dim> in(1, this->n_compositional_fields());
      MaterialModel::MaterialModelOutputs<dim> out(1, this->n_compositional_fields());
      in.position[0]=position;
      in.temperature[0]=this->get_adiabatic_conditions().temperature(position);
      in.pressure[0]=this->get_adiabatic_conditions().pressure(position);
      in.velocity[0]= Tensor<1,dim> ();
      for (unsigned int c=0; c<this->n_compositional_fields(); ++c)
        in.composition[0][c] = function->value(Point<1>(depth),c);
      in.strain_rate.resize(0); // adiabat has strain=0.
      this->get_material_model().evaluate(in, out);

      const double kappa = ( (this->get_parameters().formulation_temperature_equation ==
                              Parameters<dim>::Formulation::TemperatureEquation::reference_density_profile)
                             ?
                             out.thermal_conductivities[0] /
                             (this->get_adiabatic_conditions().density(in.position[0]) * out.specific_heat[0])
                             :
                             out.thermal_conductivities[0] / (out.densities[0] * out.specific_heat[0])
                           );

      // analytical solution for the thermal boundary layer from half-space cooling model
      const double surface_cooling_temperature = age_top > 0.0 ?
                                                 (T_surface - adiabatic_surface_temperature) *
                                                 erfc(this->get_geometry_model().depth(position) /
                                                      (2 * sqrt(kappa * age_top)))
                                                 : 0.0;
      const double bottom_heating_temperature = (age_bottom > 0.0 && this->get_adiabatic_conditions().is_initialized()) ?
                                                (T_bottom - adiabatic_bottom_temperature + subadiabaticity)
                                                * erfc((this->get_geometry_model().maximal_depth()
                                                        - this->get_geometry_model().depth(position)) /
                                                       (2 * sqrt(kappa * age_bottom)))
                                                : 0.0;

      // set the initial temperature perturbation
      // first: get the center of the perturbation, then check the distance to the
      // evaluation point. the center is supposed to lie at the center of the bottom
      // surface.
      Point<dim> mid_point;
      if (perturbation_position == "center")
        {
          if (Plugins::plugin_type_matches<const GeometryModel::SphericalShell<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::SphericalShell<dim> &shell_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::SphericalShell<dim>> (this->get_geometry_model());

              const double inner_radius = shell_geometry_model.inner_radius();
              const double half_opening_angle = numbers::PI/180.0 * 0.5 * shell_geometry_model.opening_angle();
              if (dim==2)
                {
                  // choose the center of the perturbation at half angle along the inner radius
                  mid_point(0) = inner_radius * std::sin(half_opening_angle),
                  mid_point(1) = inner_radius * std::cos(half_opening_angle);
                }
              else if (dim==3)
                {
                  // if the opening angle is 90 degrees (an eighth of a full spherical
                  // shell, then choose the point on the inner surface along the first
                  // diagonal
                  if (shell_geometry_model.opening_angle() == 90)
                    {
                      mid_point(0) = inner_radius*std::sqrt(1./3),
                      mid_point(1) = inner_radius*std::sqrt(1./3),
                      mid_point(2) = inner_radius*std::sqrt(1./3);
                    }
                  else
                    {
                      // otherwise do the same as in 2d
                      mid_point(0) = inner_radius * std::sin(half_opening_angle) * std::cos(half_opening_angle),
                      mid_point(1) = inner_radius * std::sin(half_opening_angle) * std::sin(half_opening_angle),
                      mid_point(2) = inner_radius * std::cos(half_opening_angle);
                    }
                }
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::Chunk<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::Chunk<dim> &chunk_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Chunk<dim>> (this->get_geometry_model());

              const double inner_radius = chunk_geometry_model.inner_radius();

              const double west_longitude = chunk_geometry_model.west_longitude(); // in radians
              const double longitude_range = chunk_geometry_model.longitude_range(); // in radians
              const double longitude_midpoint = west_longitude + 0.5 * longitude_range;

              if (dim==2)
                {
                  // choose the center of the perturbation at half angle along the inner radius
                  mid_point(0) = inner_radius * std::cos(longitude_midpoint),
                  mid_point(1) = inner_radius * std::sin(longitude_midpoint);
                }
              else if (dim==3)
                {

                  const double south_latitude = chunk_geometry_model.south_latitude(); // in radians
                  const double latitude_range = chunk_geometry_model.latitude_range(); // in radians
                  const double latitude_midpoint = south_latitude + 0.5 * latitude_range;
                  mid_point(0) = inner_radius * std::cos(latitude_midpoint) * std::cos(longitude_midpoint);
                  mid_point(1) = inner_radius * std::cos(latitude_midpoint) * std::sin(longitude_midpoint);
                  mid_point(2) = inner_radius * std::sin(latitude_midpoint);
                }
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::Box<dim>>(this->get_geometry_model()))
            {
              const GeometryModel::Box<dim> &box_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::Box<dim>> (this->get_geometry_model());

              // for the box geometry, choose a point at the center of the bottom face.
              // (note that the loop only runs over the first dim-1 coordinates, leaving
              // the depth variable at zero)
              mid_point = box_geometry_model.get_origin();
              for (unsigned int i=0; i<dim-1; ++i)
                mid_point(i) += 0.5 * box_geometry_model.get_extents()[i];
            }
          else if (Plugins::plugin_type_matches<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model()))
            {
              const GeometryModel::TwoMergedBoxes<dim> &two_merged_boxes_geometry_model =
                Plugins::get_plugin_as_type<const GeometryModel::TwoMergedBoxes<dim>> (this->get_geometry_model());

              // for the box geometry, choose a point at the center of the bottom face.
              // (note that the loop only runs over the first dim-1 coordinates, leaving
              // the depth variable at zero)
              mid_point = two_merged_boxes_geometry_model.get_origin();
              for (unsigned int i=0; i<dim-1; ++i)
                mid_point(i) += 0.5 * two_merged_boxes_geometry_model.get_extents()[i];
            }
          else
            AssertThrow (false,
                         ExcMessage ("Not a valid geometry model for the initial temperature model"
                                     "adiabatic."));
        }

      const double perturbation = (mid_point.distance(position) < radius) ? amplitude
                                  : 0.0;


      // add the subadiabaticity
      const double zero_depth = 0.174;
      const double nondimensional_depth = (this->get_geometry_model().depth(position) / this->get_geometry_model().maximal_depth() - zero_depth)
                                          / (1.0 - zero_depth);
      double subadiabatic_T = 0.0;
      if (nondimensional_depth > 0)
        subadiabatic_T = -subadiabaticity * nondimensional_depth * nondimensional_depth;

      // If adiabatic heating is disabled, apply all perturbations to
      // constant adiabatic surface temperature instead of adiabatic profile.
      const double temperature_profile = (this->include_adiabatic_heating())
                                         ?
                                         this->get_adiabatic_conditions().temperature(position)
                                         :
                                         adiabatic_surface_temperature;

      // return sum of the adiabatic profile, the boundary layer temperatures and the initial
      // temperature perturbation.
      return temperature_profile + surface_cooling_temperature
             + (perturbation > 0.0 ? std::max(bottom_heating_temperature + subadiabatic_T,perturbation)
                : bottom_heating_temperature + subadiabatic_T) - adiabatic_surface_temperature;
    }


    template <int dim>
    void
    Adiabatic1<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
      prm.declare_entry ("Adiabatic surface temperature", "0.",
                         Patterns::Double(),
                         "In order to make the problem in the first time step easier to "
                         "solve, we need a reasonable guess for the temperature and pressure. "
                         "To obtain it, we use an adiabatic pressure and temperature field. "
                         "This parameter describes what the `adiabatic' temperature would "
                         "be at the surface of the domain (i.e. at depth zero). Note "
                         "that this value need not coincide with the boundary condition "
                         "posed at this point. Rather, the boundary condition may differ "
                         "significantly from the adiabatic value, and then typically "
                         "induce a thermal boundary layer."
                         "\n\n"
                         "For more information, see the section in the manual that discusses "
                         "the general mathematical model.");
        prm.enter_subsection("Subduction adiabatic");
        {
          prm.declare_entry ("Age top boundary layer", "0.",
                             Patterns::Double (0.),
                             "The age of the upper thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Age bottom boundary layer", "0.",
                             Patterns::Double (0.),
                             "The age of the lower thermal boundary layer, used for the calculation "
                             "of the half-space cooling model temperature. Units: years if the "
                             "'Use years in output instead of seconds' parameter is set; "
                             "seconds otherwise.");
          prm.declare_entry ("Radius", "0.",
                             Patterns::Double (0.),
                             "The Radius (in m) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain.");
          prm.declare_entry ("Amplitude", "0.",
                             Patterns::Double (0.),
                             "The amplitude (in K) of the initial spherical temperature perturbation "
                             "at the bottom of the model domain. This perturbation will be added to "
                             "the adiabatic temperature profile, but not to the bottom thermal "
                             "boundary layer. Instead, the maximum of the perturbation and the bottom "
                             "boundary layer temperature will be used.");
          prm.declare_entry ("Position", "center",
                             Patterns::Selection ("center"),
                             "Where the initial temperature perturbation should be placed. If `center' is "
                             "given, then the perturbation will be centered along a `midpoint' of some "
                             "sort of the bottom boundary. For example, in the case of a box geometry, "
                             "this is the center of the bottom face; in the case of a spherical shell "
                             "geometry, it is along the inner surface halfway between the bounding "
                             "radial lines.");
          prm.declare_entry ("Subadiabaticity", "0.",
                             Patterns::Double (0.),
                             "If this value is larger than 0, the initial temperature profile will "
                             "not be adiabatic, but subadiabatic. This value gives the maximal "
                             "deviation from adiabaticity. Set to 0 for an adiabatic temperature "
                             "profile. Units: \\si{\\kelvin}.\n\n"
                             "The function object in the Function subsection "
                             "represents the compositional fields that will be used as a reference "
                             "profile for calculating the thermal diffusivity. "
                             "This function is one-dimensional and depends only on depth. The format of this "
                             "functions follows the syntax understood by the "
                             "muparser library, see Section~\\ref{sec:muparser-format}.");
          prm.enter_subsection("Function");
          {
            Functions::ParsedFunction<1>::declare_parameters (prm, 1);
          }
          prm.leave_subsection();
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }


    template <int dim>
    void
    Adiabatic1<dim>::parse_parameters (ParameterHandler &prm)
    {
      // we need to get the number of compositional fields here to
      // initialize the function parser. unfortunately, we can't get it
      // via SimulatorAccess from the simulator itself because at the
      // current point the SimulatorAccess hasn't been initialized
      // yet. so get it from the parameter file directly.
      adiabatic_surface_temperature   = prm.get_double ("Adiabatic surface temperature");
      prm.enter_subsection ("Compositional fields");
      const unsigned int n_compositional_fields = prm.get_integer ("Number of fields");
      prm.leave_subsection ();

      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Subduction adiabatic");
        {
          age_top_boundary_layer = prm.get_double ("Age top boundary layer");
          age_bottom_boundary_layer = prm.get_double ("Age bottom boundary layer");
          radius = prm.get_double ("Radius");
          amplitude = prm.get_double ("Amplitude");
          perturbation_position = prm.get("Position");
          subadiabaticity = prm.get_double ("Subadiabaticity");
          if (n_compositional_fields > 0)
            {
              prm.enter_subsection("Function");
              try
                {
                  function
                    = std_cxx14::make_unique<Functions::ParsedFunction<1>>(n_compositional_fields);
                  function->parse_parameters (prm);
                }
              catch (...)
                {
                  std::cerr << "ERROR: FunctionParser failed to parse\n"
                            << "\t'Initial temperature model.Adiabatic.Function'\n"
                            << "with expression\n"
                            << "\t'" << prm.get("Function expression") << "'"
                            << "More information about the cause of the parse error \n"
                            << "is shown below.\n";
                  throw;
                }

              prm.leave_subsection();
            }
        }
        prm.leave_subsection ();
      }
      prm.leave_subsection ();
    }

  }


  namespace InitialComposition
  {
    template <int dim>
    double
    Subduction2C<dim>::
    initial_composition (const Point<dim> &position, const unsigned int n_comp) const
    {
      // point position
      const std::array<double,dim> spherical_coordinates =
            aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      const double r = spherical_coordinates[0];
      const double phi = spherical_coordinates[1];

      // determine composition
      bool is_spcrust = false;
      bool is_spharz = false;
      bool is_opcrust = false;
      bool is_opharz = false;
      
      const double Ro = 6.371e6;
      const double DCS = 7.5e3;  // crustal thickness
      const double DHS = 3.520e4;  // depth of the harzburgite layer
      const double PHIC = 0.628319; // trench angular position
      const double ST=2.000e+05; // depth of the slab
      const double AC=2.000e+05; // shape of the ellipse
      const double BC=4.000e+05;
      const double tolerance = 1e-6;  // tolerance for the algorithm looking for shortest distance
      
      // a condition to determine if the point is within the ellipse envelop of
      // the initial slab
      double m = ellipse_equation_right(AC, BC, Ro*(phi-PHIC), r+BC-Ro);
      bool is_in_ellipse = (m <= 1);
      double depth_in_slab = 0.0;
      
      //spcrust
      is_spcrust = ((phi<PHIC)&&(Ro-r<=DCS));
      bool is_in_slab = (((Ro-r)<=ST) && is_in_ellipse && (phi>PHIC));
      // compute depth in slab when query point is in the ellipse
      if (is_in_slab)
      {
        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, BC/AC, 
                                                                        Ro*(phi-PHIC)/AC, (r+BC-Ro)/AC, tolerance); // distance from focus to slab surface
        depth_in_slab = sqrt(ellipse_distance_sqr(AC, BC, Ro*(phi-PHIC), 
                                                  r+BC-Ro, results.first));
      }
      is_spcrust = (is_spcrust || 
                    (is_in_slab && (depth_in_slab <= DCS)));  // crustal slab
      //spharz
      is_spharz = ((phi<PHIC)&&(Ro-r>DCS)&&(Ro-r<=DHS));
      is_spharz = (is_spharz || 
                    (is_in_slab && (depth_in_slab > DCS) && (depth_in_slab <= DHS)));  // harzburgitic slab
      //opcrust
      is_opcrust = ((phi>=PHIC)&&(Ro-r<=DCS)&&(!is_in_ellipse));
      //opharz
      is_opharz = ((phi>=PHIC)&&(Ro-r>DCS)&&(Ro-r<=DHS)&&(!is_in_ellipse));

      // assign composition
      if (is_spcrust && n_comp == 0)
        return 1.0;
      else if (is_spharz && n_comp == 1)
        return 1.0;
      else if (is_opcrust && n_comp == 2)
        return 1.0;
      else if (is_opharz && n_comp == 3)
        return 1.0;
      else
        return 0.0;
    }
    
    template <int dim>
    void
    Subduction2C<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Subduction 2d composition ellipse");
        {

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    
    template <int dim>
    void
    Subduction2C<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Subduction 2d composition ellipse");
        {
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}
    
double
ellipse_equation_right(const double AC, const double BC, const double x, const double y)
{
  // compute the value on the right side of a ellipse equation
  double m = sqrt(pow(x/AC, 2.0) + pow(y/BC, 2.0));
  return m;
}

double
ellipse_spherical_dr(const double AC, const double BC, const double x, const double y)
{
  // distance from focus to slab surface
  double a, b ,c, consf, rpc;
  if (BC > AC)
  {
    a = BC;
    b = AC;
    c = sqrt(a*a - b*b);
    rpc = sqrt(x*x + (y-c)*(y-c));
    consf = (y-c)/rpc; // angle of the focus-intersection with the semi-long axis
  }
  else
  {
    a = AC;
    b = BC;
    c = sqrt(a*a - b*b);
    rpc = sqrt((x-c)*(x-c) + y*y);
    consf = (x-c)/rpc; // angle of the focus-intersection with the semi-long axis
  }
  double e = sqrt(1 - b*b /(a*a));
  double p = a * (1 - e*e);
  c = a * e;
  double rout = p / (1 + e*consf);
  double d = abs(rout - rpc);
  return d;
}
    
    
double
half_space_cooling (const double internal_temperature, 
                    const double surface_temperature,
                    const double depth, 
                    const double thermal_diffusivity,
                    const double age)
{
  double temperature = internal_temperature - 
                       (internal_temperature - surface_temperature) * erfc(depth / (2 * sqrt(thermal_diffusivity * age)));
  return temperature;
}


std::pair<double, int> 
ellipse_distance_sqr_shortest1(const double AC, const double BC,
                               const double xp, const double yp, const double tolerance)
{
    double theta0 = 0.0;
    double theta1 = M_PI / 2.0;
    double distance_sqr_div2 = ellipse_distance_sqr_div(AC, BC, xp, 
                                                       yp, theta0);
    
    int i = 0;
    while(abs(distance_sqr_div2 / (AC * AC)) > tolerance)
    {
       double theta2 = (theta0 + theta1) / 2.0; 
       distance_sqr_div2 = ellipse_distance_sqr_div(AC, BC, xp, 
                                                    yp, theta2);
        if(distance_sqr_div2 > 0.0)
            theta1 = theta2;
        else if(distance_sqr_div2 < 0.0)
            theta0 = theta2;
        else
        {
            // 0.0, break and return
            theta0 = theta2;
            break;
        }
        i++;
    }
    // return both the result and the number of iteration 
    return std::make_pair(theta0, i);
}
    
double
ellipse_distance_sqr_div(const double AC, const double BC, const double xp, 
                                const double yp, const double theta)
{
    const double distance_sqr_div = (BC*BC - AC*AC) * sin(2.0*theta)
                                    + 2.0 * AC * xp * sin(theta) 
                                    - 2.0 * BC * yp * cos(theta);
    return distance_sqr_div;
}
    
    
double
ellipse_distance_sqr(const double AC, const double BC, const double xp, 
                    const double yp, const double theta)
{
    const double distance_sqr = pow(AC * cos(theta) - xp, 2.0) +
                                 pow(BC * sin(theta) - yp, 2.0);
    return distance_sqr;
}


// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Subduction2T,
                                              "subduction 2d temperature ellipse",
                                              "Implementation of an initial temperature for 2d subduction model")
  }
  
  namespace InitialComposition
  {
    ASPECT_REGISTER_INITIAL_COMPOSITION_MODEL(Subduction2C,
                                              "subduction 2d composition ellipse",
                                              "Implementation of an initial composition model for 2d subduction model")
  }
  
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Adiabatic1,
                                              "subduction adiabatic",
                                              "Temperature is prescribed as an adiabatic, which is consistent to use with slab temperature "
                                              "profile with upper and lower thermal boundary layers, "
                                              "whose ages are given as input parameters.")
  }
}

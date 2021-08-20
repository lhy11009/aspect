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
      // convert to spherical coordinates
      const std::array<double,dim> spherical_coordinates =
            aspect::Utilities::Coordinates::cartesian_to_spherical_coordinates(position);
      const double r = spherical_coordinates[0];
      const double phi = spherical_coordinates[1];
      double temperature = 0.0;
      if (extended_boussinesq)
      {
        if(outer_radius-r < depth_of_lith)
        {
          // lithosphere: apply a differential temperature
          const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
          temperature = slab_model(r, phi) - this->get_adiabatic_conditions().temperature(surface_point);
        }
        else
        {
          // mantle: use the adiabatic temperature
          temperature = 0.0;
        }
      }
      else // this is manually defined
        temperature = temperature_substract_adiabatic(r, temperature);
      return temperature;
    }
    
    template <int dim>
    double
    Subduction2T<dim>::
    slab_model (const double r, const double phi) const
    {
      // todo
      const double tr_rad= trench_longitude * M_PI / 180.0;
      const double max_rad= maximum_longtitude * M_PI / 180.0;
      const double year = 3600 * 24 * 365.0;
      const double v_sub_s = v_sub_plate / year;
      const double age_op_s = age_op_plate * year;
      const Point<dim> surface_point = this->get_geometry_model().representative_point(0.0);
      const double TM = this->get_adiabatic_conditions().temperature(surface_point);
      const double K=1.000e-06;

      const double tolerance = 1e-6;  // tolerance for the algorithm looking for shortest distance

      double temperature;


      // a thickness of thermal boundary on the upper slab boundary
      const double Ht_slab_in = thermal_boundary_width_factor_in * sqrt(K*ellipse_r_axis/v_sub_s);
      const double Ht_slab_out = thermal_boundary_width_factor_out * sqrt(K*ellipse_r_axis/v_sub_s);

      // a condition to determine if the point is within the ellipse envelop of
      // the initial slab
      double m = ellipse_equation_right(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), r+ellipse_r_axis-outer_radius);
      bool is_in_ellipse = (m <= 1);

      if (phi<=0.0)
      {
        // ghost domain
        temperature = TM;
      }
      else if ((phi<=tr_rad)&&(phi>0))
      {
        // subducting plate
        temperature = half_space_cooling(TM, outer_temperature, abs(outer_radius-r), K, outer_radius*phi/v_sub_s);
      }
      else if (((outer_radius-r)<=depth_of_slab)
                && is_in_ellipse
                &&(phi>tr_rad))
                {
        // slab
        // temperature of the plate at the same depth
        const double plate_temperature = half_space_cooling(TM, outer_temperature, abs(outer_radius-r), K, outer_radius*phi/v_sub_s);
        
        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, ellipse_r_axis/ellipse_phi_axis, 
                                                                        outer_radius*(phi-tr_rad)/ellipse_phi_axis, (r+ellipse_r_axis-outer_radius)/ellipse_phi_axis, tolerance); // distance from focus to slab surface
        const double depth_in_slab = sqrt(ellipse_distance_sqr(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), 
                                                               r+ellipse_r_axis-outer_radius, results.first));
        const double slab_temperature = half_space_cooling(TM, outer_temperature, depth_in_slab, K, outer_radius*phi/v_sub_s);

        // cooling temperature as the smaller of these two

        const double cooling_temperature = std::min(plate_temperature, slab_temperature);

        // a perturbation on the upper surface
        // This is implemented as tranform from internal temperature of the slab to surface temperature.
        // The surface temperature is taken as the average of outer_temperature and temperature of overiding plate.
        const double over_plate_temperature = half_space_cooling(TM, outer_temperature, abs(outer_radius-r), K, age_op_s);
        // method for perturbation
        double perturbation_fucntion_value = sin(M_PI/2.0*depth_in_slab/Ht_slab_in);
        const double perturbation_temperature = (depth_in_slab < Ht_slab_in)?
                                          (over_plate_temperature -cooling_temperature)
                                          * (1- perturbation_fucntion_value) / 2.0:
                                          0.0;
        
        // Final temperature is the sum of the two parts.
        temperature = cooling_temperature + perturbation_temperature;
                }
      else if (((outer_radius-r)>depth_of_slab)
               &&((outer_radius-r)<=depth_of_lith)
               && is_in_ellipse)
               {
        // temperature below the slab tip
        // A half-space cooling temperature is assigned once more but
        // an intermediate temperature is taken for the surface temperature.
        // depth of this interface
        // This have the same expression to depth_in_slab, except that this is actually below the initial slab

        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, ellipse_r_axis/ellipse_phi_axis, 
                                                                        outer_radius*(phi-tr_rad)/ellipse_phi_axis, (r+ellipse_r_axis-outer_radius)/ellipse_phi_axis, tolerance); // distance from focus to slab surface
        const double depth_in_slab = sqrt(ellipse_distance_sqr(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), 
                                                               r+ellipse_r_axis-outer_radius, results.first));

        // Intermediate temperature on this interface
        const double T_interface = outer_temperature+(TM-outer_temperature)*(abs(outer_radius-r)-depth_of_slab)/(depth_of_lith-depth_of_slab);

        const double cooling_temperature = half_space_cooling(TM, T_interface, depth_in_slab, K, outer_radius*phi/v_sub_s);
        
        // temperature overiding plate
        const double over_plate_temperature = half_space_cooling(TM, outer_temperature, abs(outer_radius-r), K, age_op_s);
        
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
      else if (max_rad-phi>0)
      {
        // overiding plate
        const double over_plate_temperature = half_space_cooling(TM, outer_temperature, abs(outer_radius-r), K, age_op_s);
        
        // This formula includes absolute value just to make sure it is consistent with the previous 'depth_in_slab'
        // the depth out of the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, ellipse_r_axis/ellipse_phi_axis, 
                                                                        outer_radius*(phi-tr_rad)/ellipse_phi_axis, (r+ellipse_r_axis-outer_radius)/ellipse_phi_axis, tolerance); // distance from focus to slab surface
        const double depth_out_slab = sqrt(ellipse_distance_sqr(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), 
                                                               r+ellipse_r_axis-outer_radius, results.first));
        
        // a perturbation on the upper surface of the slab
        // This is implemented as tranform from temperature of the overiding plate to surface temperature.
        // Out of the upper boundary layer of the slab, this is 0.0
        double perturbation_temperature = 0.0;
        if (depth_out_slab<Ht_slab_out)
        { 
          if ((outer_radius-r)<depth_of_slab)
          {
            // Surface temperature is taken as the average of TS and temperature of overiding plate.
            perturbation_temperature = (outer_temperature - over_plate_temperature)
                                        * (1-sin(M_PI/2.0*depth_out_slab/Ht_slab_out)) / 2.0;
          }
          else if ((outer_radius-r)<depth_of_lith)
          {
            // Surface temperature is taken as the average of an intermediate temperature
            // and temperature of overiding plate.
            // Intermediate temperature on this interface
            const double T_interface = outer_temperature+(TM-outer_temperature)*(abs(outer_radius-r)-depth_of_slab)/(depth_of_lith-depth_of_slab);

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
      double temperature = s(outer_radius-r, true);
      return temperature;
    } 

    template <int dim>
    double
    Subduction2T<dim>::
    temperature_substract_adiabatic(const double r, const double temperature) const
    {
      const double outer_radius = 6371e3;
      const double Dadiabatic = 190e3;
      const double adiabatic_gradient = 0.25/1e3;
      if(r<outer_radius-Dadiabatic)
        return temperature - adiabatic_gradient * (outer_radius-r - Dadiabatic);
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
      prm.enter_subsection("Boundary temperature model");
      {
          prm.enter_subsection("Spherical constant");
          {
            prm.declare_entry ("Outer temperature", "0.",
                               Patterns::Double (),
                               "Temperature at the outer boundary (lithosphere water/air). Units: \\si{\\kelvin}.");
          }
          prm.leave_subsection();
      }
      prm.leave_subsection();
      prm.enter_subsection("Geometry model");
      {
          prm.enter_subsection("Chunk");
          {
              prm.declare_entry ("Chunk outer radius", "1.",
                                 Patterns::Double (0.),
                                 "Radius at the top surface of the chunk. Units: \\si{\\meter}.");
              prm.declare_entry ("Chunk maximum longitude", "1.",
                                Patterns::Double (-180., 360.), // enables crossing of either hemisphere
                                "Maximum longitude of the chunk. Units: degrees.");
          }
          prm.leave_subsection();
      }
      prm.leave_subsection();
      prm.enter_subsection ("Initial temperature model");
      { 
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d_Vs_YT16.txt");
        prm.enter_subsection("Subduction 2d temperature ellipse");
        {
         prm.declare_entry ("Thermal boundary width factor in", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
         prm.declare_entry ("Thermal boundary width factor out", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
         prm.declare_entry("Extended boussinesq", "false", Patterns::Bool(),
                            "Initial temperature for extended boussinesq approximation");
         // todo
         prm.declare_entry ("Depth of lithsphere", "270e3",
                            Patterns::Double (),
                            "This parameter controls the depth of the lithosphere");
         prm.declare_entry ("Depth of slab", "200e3",
                            Patterns::Double (),
                            "This parameter controls the depth of the slab tip");
         prm.declare_entry ("Trench longitude", "36.0",
                            Patterns::Double (),
                            "This parameter controls the position of trench");
         prm.declare_entry ("Ellipse phi axis", "2.000e+05",
                            Patterns::Double (),
                            "This parameter controls the long axis of ellipse");
         prm.declare_entry ("Ellipse r axis", "4.000e+05",
                            Patterns::Double (),
                            "This parameter controls the short axis of ellipse");
         prm.declare_entry ("Velocity subducting plate", "0.05",
                            Patterns::Double (),
                            "This parameter controls the velocity of the subducting plate (m/year)");
         prm.declare_entry ("Age overiding plate", "40e6",
                            Patterns::Double (),
                            "This parameter controls the age of the overideing plate (year)");
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
      prm.enter_subsection("Boundary temperature model");
      {
          prm.enter_subsection("Spherical constant");
          {
            outer_temperature = prm.get_double ("Outer temperature");
          }
          prm.leave_subsection();
      }
      prm.leave_subsection();
      
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          outer_radius = prm.get_double ("Chunk outer radius");
          maximum_longtitude = prm.get_double ("Chunk maximum longitude");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Subduction 2d temperature ellipse");
        {
          thermal_boundary_width_factor_in = prm.get_double("Thermal boundary width factor in");
          thermal_boundary_width_factor_out = prm.get_double("Thermal boundary width factor out");
          extended_boussinesq = prm.get_bool("Extended boussinesq");
          // todo
          depth_of_lith = prm.get_double("Depth of lithsphere");
          depth_of_slab = prm.get_double("Depth of slab");
          trench_longitude = prm.get_double("Trench longitude");
          ellipse_phi_axis = prm.get_double("Ellipse phi axis");
          ellipse_r_axis = prm.get_double("Ellipse r axis");
          v_sub_plate = prm.get_double("Velocity subducting plate");
          age_op_plate = prm.get_double("Age overiding plate");
        }
        prm.leave_subsection();

        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
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

      // todo 
      const double tr_rad= trench_longitude * M_PI / 180.0;
      const double tolerance = 1e-6;  // tolerance for the algorithm looking for shortest distance
      
      // a condition to determine if the point is within the ellipse envelop of
      // the initial slab
      double m = ellipse_equation_right(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), r+ellipse_r_axis-outer_radius);
      bool is_in_ellipse = (m <= 1);
      double depth_in_slab = 0.0;
      
      //spcrust
      is_spcrust = ((phi<tr_rad)&&(outer_radius-r<=crustal_thickness));
      bool is_in_slab = (((outer_radius-r)<=depth_of_slab) && is_in_ellipse && (phi>tr_rad));
      // compute depth in slab when query point is in the ellipse
      if (is_in_slab)
      {
        // the depth withing the initial slab
        std::pair<double, int> results = ellipse_distance_sqr_shortest1(1.0, ellipse_r_axis/ellipse_phi_axis, 
                                                                        outer_radius*(phi-tr_rad)/ellipse_phi_axis, (r+ellipse_r_axis-outer_radius)/ellipse_phi_axis, tolerance); // distance from focus to slab surface
        depth_in_slab = sqrt(ellipse_distance_sqr(ellipse_phi_axis, ellipse_r_axis, outer_radius*(phi-tr_rad), 
                                                  r+ellipse_r_axis-outer_radius, results.first));
      }
      is_spcrust = (is_spcrust || 
                    (is_in_slab && (depth_in_slab <= crustal_thickness)));  // crustal slab
      //spharz
      is_spharz = ((phi<tr_rad)&&(outer_radius-r>crustal_thickness)&&(outer_radius-r<=harz_thickness));
      is_spharz = (is_spharz || 
                    (is_in_slab && (depth_in_slab > crustal_thickness) && (depth_in_slab <= harz_thickness)));  // harzburgitic slab
      //opcrust
      is_opcrust = ((phi>=tr_rad)&&(outer_radius-r<=crustal_thickness)&&(!is_in_ellipse));
      //opharz
      is_opharz = ((phi>=tr_rad)&&(outer_radius-r>crustal_thickness)&&(outer_radius-r<=harz_thickness)&&(!is_in_ellipse));

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
      prm.enter_subsection("Geometry model");
      {
          prm.enter_subsection("Chunk");
          {
              prm.declare_entry ("Chunk outer radius", "1.",
                                 Patterns::Double (0.),
                                 "Radius at the top surface of the chunk. Units: \\si{\\meter}.");
          }
          prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Subduction 2d composition ellipse");
        {
          // todo
          prm.declare_entry ("Depth of slab", "200e3",
                              Patterns::Double (),
                              "This parameter controls the depth of the slab tip");
          prm.declare_entry ("Trench longitude", "36.0",
                              Patterns::Double (),
                              "This parameter controls the position of trench");
          prm.declare_entry ("Crustal thickness", "7.5e3",
                              Patterns::Double (),
                              "This parameter controls the crustal thickness");
          prm.declare_entry ("Harzburgite thickness", "3.520e4",
                              Patterns::Double (),
                              "This parameter controls the harzburgite layer thickness");
          prm.declare_entry ("Ellipse phi axis", "2.000e+05",
                             Patterns::Double (),
                             "This parameter controls the long axis of ellipse");
          prm.declare_entry ("Ellipse r axis", "4.000e+05",
                             Patterns::Double (),
                             "This parameter controls the short axis of ellipse");

        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
    
    template <int dim>
    void
    Subduction2C<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry model");
      {
        prm.enter_subsection("Chunk");
        {
          outer_radius = prm.get_double ("Chunk outer radius");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();

      prm.enter_subsection("Initial composition model");
      {
        prm.enter_subsection("Subduction 2d composition ellipse");
        {
          // todo
          depth_of_slab = prm.get_double("Depth of slab");
          trench_longitude = prm.get_double("Trench longitude");
          crustal_thickness = prm.get_double("Crustal thickness");
          harz_thickness = prm.get_double("Harzburgite thickness");
          ellipse_phi_axis = prm.get_double("Ellipse phi axis");
          ellipse_r_axis = prm.get_double("Ellipse r axis");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }
  }
}
    
double
ellipse_equation_right(const double ellipse_phi_axis, const double ellipse_r_axis, const double x, const double y)
{
  // compute the value on the right side of a ellipse equation
  double m = sqrt(pow(x/ellipse_phi_axis, 2.0) + pow(y/ellipse_r_axis, 2.0));
  return m;
}

double
ellipse_spherical_dr(const double ellipse_phi_axis, const double ellipse_r_axis, const double x, const double y)
{
  // distance from focus to slab surface
  double a, b ,c, consf, rpc;
  if (ellipse_r_axis > ellipse_phi_axis)
  {
    a = ellipse_r_axis;
    b = ellipse_phi_axis;
    c = sqrt(a*a - b*b);
    rpc = sqrt(x*x + (y-c)*(y-c));
    consf = (y-c)/rpc; // angle of the focus-intersection with the semi-long axis
  }
  else
  {
    a = ellipse_phi_axis;
    b = ellipse_r_axis;
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
ellipse_distance_sqr_shortest1(const double ellipse_phi_axis, const double ellipse_r_axis,
                               const double xp, const double yp, const double tolerance)
{
    double theta0 = 0.0;
    double theta1 = M_PI / 2.0;
    double distance_sqr_div2 = ellipse_distance_sqr_div(ellipse_phi_axis, ellipse_r_axis, xp, 
                                                       yp, theta0);
    
    int i = 0;
    while(abs(distance_sqr_div2 / (ellipse_phi_axis * ellipse_phi_axis)) > tolerance)
    {
       double theta2 = (theta0 + theta1) / 2.0; 
       distance_sqr_div2 = ellipse_distance_sqr_div(ellipse_phi_axis, ellipse_r_axis, xp, 
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
ellipse_distance_sqr_div(const double ellipse_phi_axis, const double ellipse_r_axis, const double xp, 
                                const double yp, const double theta)
{
    const double distance_sqr_div = (ellipse_r_axis*ellipse_r_axis - ellipse_phi_axis*ellipse_phi_axis) * sin(2.0*theta)
                                    + 2.0 * ellipse_phi_axis * xp * sin(theta) 
                                    - 2.0 * ellipse_r_axis * yp * cos(theta);
    return distance_sqr_div;
}
    
    
double
ellipse_distance_sqr(const double ellipse_phi_axis, const double ellipse_r_axis, const double xp, 
                    const double yp, const double theta)
{
    const double distance_sqr = pow(ellipse_phi_axis * cos(theta) - xp, 2.0) +
                                 pow(ellipse_r_axis * sin(theta) - yp, 2.0);
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
}

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
#include "subduction_temperature2d.h"
#include "spline.h"
#include <aspect/material_model/interface.h>
#include <aspect/initial_composition/interface.h>
#include <aspect/adiabatic_conditions/interface.h>
#include <aspect/simulator_access.h>
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
      const bool extended_boussinesq = false;
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
      Point<dim> polar_position = cartesian_to_polar(position);
      const double r = polar_position[0];
      const double phi = polar_position[1];
      double temperature;
      if(Ro-r < SM) 
        temperature = slab_model(r, phi);
      else
      {
        temperature = katrina_mantle_temperature(r);
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
        const double depth_in_slab = ellipse_spherical_dr(AC, BC, Ro*(phi-PHIC), r+BC-Ro); // distance from focus to slab surface
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
        const double depth_in_slab = ellipse_spherical_dr(AC, BC, Ro*(phi-PHIC), r+BC-Ro); // distance from focus to slab surface

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
        const double depth_out_slab = ellipse_spherical_dr(AC, BC, Ro*(phi-PHIC), r+BC-Ro); // distance from focus to slab surface
        
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
    ellipse_equation_right(const double AC, const double BC, const double x, const double y) const
    {
      // compute the value on the right side of a ellipse equation
      double m = sqrt(pow(x/AC, 2.0) + pow(y/BC, 2.0));
      return m;
    }
    
    template <int dim>
    double
    Subduction2T<dim>::
    ellipse_spherical_dr(const double AC, const double BC, const double x, const double y) const
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
    double
    Subduction2T<dim>::
    half_space_cooling (const double internal_temperature, 
                        const double surface_temperature,
                        const double depth, 
                        const double thermal_diffusivity,
                        const double age) const
    {
      double temperature = internal_temperature - 
                           (internal_temperature - surface_temperature) * erfc(depth / (2 * sqrt(thermal_diffusivity * age)));
      return temperature;
    }

    template <int dim>
    void
    Subduction2T<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      { 
        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d_Vs_YT16.txt");
        prm.enter_subsection("Subduction 2d temperature");
        {
         // todo
         prm.declare_entry ("Thermal boundary width factor in", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
         prm.declare_entry ("Thermal boundary width factor out", "1.0",
                            Patterns::Double (),
                            "This parameter controls the width of the slab thermal boundary");
        }
        prm.leave_subsection();
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Subduction2T<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        prm.enter_subsection("Subduction 2d temperature");
        {
          // todo
          thermal_boundary_width_factor_in = prm.get_double("Thermal boundary width factor in");
          thermal_boundary_width_factor_out = prm.get_double("Thermal boundary width factor out");
          //semix = ;
          //semiy = ;
        }
        prm.leave_subsection();

        Utilities::AsciiDataBase<dim>::parse_parameters(prm);
      }
      prm.leave_subsection();
    }
    
    template <int dim>
    Point<dim>
    Subduction2T<dim>::cartesian_to_polar(const Point<dim> &cartesian_position) const
    {
      Point<dim> polar_position;
      double R = sqrt(pow(cartesian_position[0], 2.0) + pow(cartesian_position[1], 2.0));
      double theta = atan2(cartesian_position[1], cartesian_position[0]);
      if (theta < 0.0){
        theta += (2 * M_PI);
      }
      polar_position[0] = R;
      polar_position[1] = theta;
      return polar_position;
    } 
  }
}

// explicit instantiations
namespace aspect
{
  namespace InitialTemperature
  {
    ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Subduction2T,
                                              "subduction 2d temperature",
                                              "Implementation of an initial temperature for 2d subduction model")
  }
}

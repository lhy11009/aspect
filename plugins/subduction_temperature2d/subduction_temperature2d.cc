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
      const double RC=4.000e+05;
      const double ST=2.000e+05;
      const double SM=2.700e+05;
      const double VSUB=1.58549e-09;
      const double AGEOP=1.26144e15;
      const double TS=2.730e+02;
      const double TM=1.673e+03; 
      const double K=1.000e-06;

      double temperature;


      // a thickness of thermal boundary on the upper slab boundary
      const double Ht_slab = sqrt(K*RC/VSUB);

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
                &&(Ro*phi<=Ro*PHIC+std::min(RC,sqrt(abs(2*RC*(Ro-r)-pow(Ro-r, 2.0)))))
                &&(phi>PHIC))
                {
        // slab
        // temperature of the plate at the same depth
        const double plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, Ro*phi/VSUB);
        
        // the depth withing the initial slab
        const double depth_in_slab = RC-sqrt(pow(Ro*(phi-PHIC), 2.0)+pow((Ro-r)-RC, 2.0));
        const double slab_temperature = half_space_cooling(TM, TS, depth_in_slab, K, Ro*phi/VSUB);

        // cooling temperature as the smaller of these two

        const double cooling_temperature = std::min(plate_temperature, slab_temperature);

        // a perturbation on the upper surface
        // This is implemented as tranform from internal temperature of the slab to surface temperature.
        // The surface temperature is taken as the average of TS and temperature of overiding plate.
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        // todo method for perturbation
        double perturbation_fucntion_value = sin(M_PI/2.0*depth_in_slab/Ht_slab);
        const double perturbation_temperature = (depth_in_slab < Ht_slab)?
                                          (over_plate_temperature -cooling_temperature)
                                          * (1- perturbation_fucntion_value/ 2.0:
                                          0.0;
        
        // Final temperature is the sum of the two parts.
        temperature = cooling_temperature + perturbation_temperature;
                }
      else if (((Ro-r)>ST)
               &&((Ro-r)<=SM)
               &&(Ro*phi<=Ro*PHIC+std::min(RC,sqrt(abs(2*RC*(Ro-r)-pow(Ro-r, 2.0))))))
               {
        // temperature below the slab tip
        // A half-space cooling temperature is assigned once more but
        // an intermediate temperature is taken for the surface temperature.
        // depth of this interface
        // This have the same expression to depth_in_slab, except that this is actually below the initial slab
        const double depth_interface = (RC-sqrt(pow(Ro*(phi-PHIC), 2.0)+pow(abs(Ro-r)-RC, 2.0)));

        // Intermediate temperature on this interface
        const double T_interface = TS+(TM-TS)*(abs(Ro-r)-ST)/(SM-ST);

        const double cooling_temperature = half_space_cooling(TM, T_interface, depth_interface, K, Ro*phi/VSUB);
        
        // the depth withing the initial slab
        const double depth_in_slab = RC-sqrt(pow(Ro*(phi-PHIC), 2.0)+pow((Ro-r)-RC, 2.0));
        
        // temperature overiding plate
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        
        // a perturbation on the upper surface
        // a thickness of thermal boundary on top
        // This is implemented as tranform from internal temperature of the slab to surface temperature.
        // The surface temperature is taken as the average of the intermediate temperature,
        // and temperature of overiding plate.
        const double perturbation_temperature = (depth_in_slab < Ht_slab)?
                                          (over_plate_temperature -cooling_temperature)
                                          * (1-sin(M_PI/2.0*depth_in_slab/Ht_slab)) / 2.0:
                                          0.0;

        // Final temperature is the sum of the two parts.
        temperature = cooling_temperature + perturbation_temperature;
               }
      else if (PHIM-phi>0)
      {
        // overiding plate
        const double over_plate_temperature = half_space_cooling(TM, TS, abs(Ro-r), K, AGEOP);
        
        // the depth withing the initial slab
        // This formula includes absolute value just to make sure it is consistent with the previous 'depth_in_slab'
        const double depth_out_slab = abs(RC-sqrt(pow(Ro*(phi-PHIC), 2.0)+pow((Ro-r)-RC, 2.0)));
        
        // a perturbation on the upper surface of the slab
        // This is implemented as tranform from temperature of the overiding plate to surface temperature.
        // Out of the upper boundary layer of the slab, this is 0.0
        double perturbation_temperature = 0.0;
        if (depth_out_slab<Ht_slab)
        { 
          if ((Ro-r)<ST)
          {
            // Surface temperature is taken as the average of TS and temperature of overiding plate.
            perturbation_temperature = (TS - over_plate_temperature)
                                        * (1-sin(M_PI/2.0*depth_out_slab/Ht_slab)) / 2.0;
          }
          else if ((Ro-r)<SM)
          {
            // Surface temperature is taken as the average of an intermediate temperature
            // and temperature of overiding plate.
            // Intermediate temperature on this interface
            const double T_interface = TS+(TM-TS)*(abs(Ro-r)-ST)/(SM-ST);

            perturbation_temperature = (T_interface - over_plate_temperature)
                                        * (1-sin(M_PI/2.0*depth_out_slab/Ht_slab)) / 2.0;
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
        //todo
        /*
        prm.declare_entry ("Remove temperature heterogeneity down to specified depth",
                           boost::lexical_cast<std::string>(-std::numeric_limits<double>::max()),
                           Patterns::Double (),
                           "This will remove temperature variations prescribed by the input model "
                           "down to the specified depth (in meters). Note that your resolution has "
                           "to be adequate to capture this cutoff. For example if you specify a depth "
                           "of 660km, but your closest spherical depth layers are only at 500km and "
                           "750km (due to a coarse resolution) it will only remove heterogeneities "
                           "down to 500km. Similar caution has to be taken when using adaptive meshing.");
        prm.declare_entry ("Set reference temperature down to specified depth", "1600",
                           Patterns::Double (),
                           "This parameter sets the a constant value of temperature down to the specified depth.");
        prm.declare_entry ("Use Yamauchi and Takei parameterization", "true",
                           Patterns::Bool(),
                           "This parameter determines whether to use the anelasticity model of "
                           "Yamauchi & Takei (2016) to convert absolute Vs into temperature");
        prm.declare_entry ("Use original density and frequency model of Yamauchi and Takei", "true",
                           Patterns::Bool(),
                           "Use original density and frequency model of Yamauchi & Takei (2016) where density"
                           "has simple pressure-dependence and shear wave is period set at 100s");
        */

        Utilities::AsciiDataBase<dim>::declare_parameters(prm,
                                                          "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
                                                          "box_2d_Vs_YT16.txt");
      }
      prm.leave_subsection();
    }


    template <int dim>
    void
    Subduction2T<dim>::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection ("Initial temperature model");
      {
        // todo
        /*
        no_perturbation_depth   = prm.get_double ("Remove temperature heterogeneity down to specified depth");
        reference_temperature   = prm.get_double ("Set reference temperature down to specified depth");
        use_yamauchi_takei = prm.get_bool ("Use Yamauchi and Takei parameterization");
        use_original_model = prm.get_bool ("Use original density and frequency model of Yamauchi and Takei");
        */

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

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


#ifndef _aspect_initial_temperature_anelasticity_temperature_h
#define _aspect_initial_temperature_anelasticity_temperature_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/utilities.h>
#include <aspect/simulator_access.h>
#include <cmath>
#include <algorithm>
#include <functional>

#include <deal.II/base/parsed_function.h>
namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that describes an initial temperature field for a 2D or 3D shear wave velocity (Vs) model.
     * Vs values are converted to temperature using the anelasticity parameterization of Yamauchi & Takei (2016).
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Subduction2T : public Utilities::AsciiDataInitial<dim>, public Interface<dim>
    {
      public:
        /**
         * Constructor. Initialize variables.
         */
        Subduction2T ();

        /**
         * Initialization function. This function is called once at the
         * beginning of the program. Checks preconditions.
         */
        void
        initialize ();

        // Avoid -Woverloaded-virtual:
        using Utilities::AsciiDataInitial<dim>::initialize;

        /**
         * Return the boundary temperature as a function of position. For the
         * current class, this function returns value from the text files.
         */
        double
        initial_temperature (const Point<dim> &position) const;

        /**
         * consist model temperature with mantle phase transition
        */
        double
        below_lith_model (const double r, 
                          const double mantle_temperature, 
                          const double lith_depth,
                          const std::vector<double>& phases_below_lith_depth,
                          const std::vector<double>& phases_below_lith_temperature) const;
        
        /**
         * temperature from a slab temperature model 
         */ 
        double slab_model (const double r, const double phi) const;

        double katrina_mantle_temperature(const double r) const;
 
        /**
         * substract adiabatic temperature if needed
         */ 
        double
        temperature_substract_adiabatic(const double r, const double temperature) const;
        
        /**
         * Declare the parameters for the input files.
         */
        static
        void
        declare_parameters (ParameterHandler  &prm);

        /**
         * Read the parameters from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm);

      private:


        /**
         * todo
         * method of slab boudnary
         */
        // method_for_slab_boundary;
        // todo
        double thermal_boundary_width_factor_in;
        double thermal_boundary_width_factor_out;
        double adiabatic_surface_temperature;
        double adiabatic_surface_temperature_gradient;

        // if using extended boussinesq approximation
        bool extended_boussinesq;
    };
    
    /**
     * A class that implements adiabatic initial conditions for the
     * temperature field and, optional, upper and lower thermal boundary
     * layers calculated using the half-space cooling model. The age of the
     * boundary layers are input parameters.
     *
     * @ingroup InitialTemperatures
     */
    template <int dim>
    class Adiabatic1 : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial temperature as a function of position.
         */
        double initial_temperature (const Point<dim> &position) const override;

        /**
         * Declare the parameters this class takes through input files.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
        /**
         * Age of the upper thermal boundary layer at the surface of the
         * model. If set to zero, no boundary layer will be present in the
         * model.
         */
        double age_top_boundary_layer;
        /* Age of the lower thermal boundary layer. */
        double age_bottom_boundary_layer;

        /**
         * Radius (in m) of the initial temperature perturbation at the bottom
         * of the model domain.
         */
        double radius;
        /**
         * Amplitude (in K) of the initial temperature perturbation at the
         * bottom of the model domain.
         */
        double amplitude;
        /*
         * Position of the initial temperature perturbation (in the
         * center or at the boundary of the model domain).
         */
        std::string perturbation_position;

        /*
         * Deviation from adiabaticity in a subadiabatic mantle
         * temperature profile. 0 for an adiabatic temperature
         * profile.
         */
        double subadiabaticity;

        /**
         * A function object representing the compositional fields that will
         * be used as a reference profile for calculating the thermal
         * diffusivity. The function depends only on depth.
         */
        std::unique_ptr<Functions::ParsedFunction<1> > function;
        
        /**
         * adiabatic surface temperature
         */ 
        double adiabatic_surface_temperature;
    };
  }

  namespace InitialComposition
  {
    using namespace dealii;

    /**
     * A class that implements initial conditions for the compositional fields
     * based on a functional description provided in the input file.
     *
     * @ingroup InitialCompositionModels
     */
    template <int dim>
    class Subduction2C : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the initial composition as a function of position and number
         * of compositional field.
         */
        double initial_composition (const Point<dim> &position, const unsigned int n_comp) const override;

        /**
         * Declare the parameters this class takes through input files. The
         * default implementation of this function does not describe any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        static
        void
        declare_parameters (ParameterHandler &prm);

        /**
         * Read the parameters this class declares from the parameter file.
         * The default implementation of this function does not read any
         * parameters. Consequently, derived classes do not have to overload
         * this function if they do not take any runtime parameters.
         */
        void
        parse_parameters (ParameterHandler &prm) override;

      private:
    };
  }
}
        
/**
 * Calculate temperature from a half-space cooling model
 * We use the internal temperature and the surface temperature so that it is easier to plug in.
 */ 
double 
half_space_cooling (const double internal_temperature, 
                    const double surface_temperature,
                    const double depth, 
                    const double thermal_diffusivity,
                    const double age);
        

std::pair<double, int> 
ellipse_distance_sqr_shortest1(const double AC, const double BC,
                               const double xp, const double yp, const double tolerance);
        
/**
 * compute the distance derivative to theta
 */ 
double
ellipse_distance_sqr_div(const double AC, const double BC, const double xp, 
                         const double yp, const double theta);
        
        
/**
 * compute the distance between a point in the ellipse to a point on the surface
 */
double
ellipse_distance_sqr(const double AC, const double BC, const double xp, 
                     const double yp, const double theta);
        
/**
 * get the right side of a ellipse
 */
double 
ellipse_equation_right(const double AC, const double BC, const double x, const double y);

/**
 * derive the distance between two points allinged with the ellipse focus
 */
double
ellipse_spherical_dr(const double AC, const double BC, const double x, const double y);
        

#endif

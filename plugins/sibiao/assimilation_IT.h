/*
  Copyright (C) 2011 - 2018 by the authors of the ASPECT code.

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

#ifndef _aspect_initial_temperature_assimilation_h
#define _aspect_initial_temperature_assimilation_h

#include <aspect/initial_temperature/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/utilities.h>
#include <aspect/compat.h>

#include <array>
#include <deal.II/base/function_lib.h>

namespace aspect
{
  namespace InitialTemperature
  {
    using namespace dealii;

    /**
     * A class that implements a prescribed temperature field determined from
     * a Ascii age data file.
     *
     * @ingroup InitialTemperatures
    */
    namespace internal
	{
    /**
          *AssimilationLookup handles all kinds of tasks around looking up a certain
          * top thermal boundary condition from a Gplates surface age file.
          */
    template <int dim>
          class AssimilationLookup
          {
            public:

              /**
               * This lookup is currently same as part of the AsciiDataLookup.
               * Initialize all members.
               */
    	       AssimilationLookup(const unsigned int components,
    	                          const double scale_factor);

              /**
               * Outputs the Assimilation module information at model start.
               */
              std::string
              screen_output() const;

              /**
               * Loads an ascii age data file. Throws an exception if the
               * file does not exist.
               */
              void
			  load_file(const std::string &filename, const MPI_Comm &comm);

              /**
               * Return the ascii age file in cartesian coordinates.
               * Takes as input the position. Actual velocity interpolation is
               * performed in spherical coordinates.
               * @param position The current position to compute age
              */
              double
              get_age(const Point<dim> &position, const unsigned int component) const;


            private:
              /**
               * The number of data components read in (=columns in the data file).
                */
                   unsigned int components;

                   /**
                    * Interpolation functions to access the data.
                    * Either InterpolatedUniformGridData or InterpolatedTensorProductGridData;
                    * the type is determined from the grid specified in the data file.
                    */
                    std::vector<std::unique_ptr<Function<dim>>> data;

                    /**
                      * The names of the data components in the columns of the read file.
                      * Does not contain any strings if none are provided in the first
                      * uncommented line of the file.
                      */
                     std::vector<std::string> data_component_names;

                     /**
                       * The coordinate values in each direction as specified in the data file.
                       */
                     std::array<std::vector<double>,dim> coordinate_values;

                     /**
                       * The maximum value of each component
                       */
                     std::vector<double> maximum_component_value;

                     /**
                       * The min and max of the coordinates in the data file.
                       */
                     std::array<std::pair<double,double>,dim> grid_extent;

                     /**
                       * Number of points in the data grid as specified in the data file.
                       */
                     TableIndices<dim> table_points;

                      /**
                        * Scales the data boundary condition by a scalar factor. Can be used
                        * to transform the unit of the data.
                        */
                     const double scale_factor;

                      /**
                        * Stores whether the coordinate values are equidistant or not,
                        * this determines the type of data function stored.
                        */
                      bool coordinate_values_are_equidistant;

                      /**
                        * Computes the table indices of each entry in the input data file.
                        * The index depends on dim, grid_dim and the number of components.
                        */
                      TableIndices<dim>
                      compute_table_indices(const unsigned int i) const;

          };
	}
      template <int dim>
      class Assimilation : public Interface<dim>, public SimulatorAccess<dim>
      {
        public:
          /**
           * Empty Constructor.
           */
          Assimilation ();

          /**
           * Initialization function. This function is called once at the
           * beginning of the program. Checks preconditions.
           */
          virtual
          void
          initialize ();

          /**
           * A function that is called at the beginning of each time step. For
           * the current plugin, this function loads the next age files if
           * necessary and outputs a warning if the end of the set of age
           * files is reached.
           */
          virtual
          void
          update ();

          /**
            * Returns the real age component at the given position.
          */
          double
		  ascii_age (const Point<dim> &position) const;

          /**
           * Return the initial temperature as a function of position. For the
           * current class, this function returns value from the text files.
           */
          double
          initial_temperature (const Point<dim> &position) const;

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
          parse_parameters (ParameterHandler &prm);

        private:
          /**
            *The parameters needed for the age data assimilation model
            */
          /**
            * A variable that stores the currently used data file of a series. It
            * gets updated if necessary by update().
            */
           int current_file_number;

           /**
            * Time from which on the data file with number 'First data file
            * number' is used as boundary condition. Previous to this time, 0 is
            * returned for every field. Depending on the setting of the global
            * 'Use years in output instead of seconds' flag in the input file,
            * this number is either interpreted as seconds or as years."
            */
           double first_data_file_model_time;

           /**
            * Number of the first data file to be loaded when the model time is
            * larger than 'First data file model time'.
            */
           int first_data_file_number;

           /**
            * In some cases the boundary files are not numbered in increasing but
            * in decreasing order (e.g. 'Ma BP'). If this flag is set to 'True'
            * the plugin will first load the file with the number 'First data
            * file number' and decrease the file number during the model run.
            */
           bool decreasing_file_order;

          /**
           * Directory in which the ascii age files are present.
           */
          std::string data_directory;

          /**
           * First part of filename of ascii age files. The files have to have
           * the pattern age_file_name.n.dat where n is the number of the
           * current timestep (starts from 0).
           */
          std::string data_file_name;

          /**
           * Time in model units (depends on other model inputs) between two
           * velocity files.
           */
          double data_file_time_step;

          /**
           * Weight between ge file n and n+1 while the current time is
           * between the two values t(n) and t(n+1).
           */
          double time_weight;

          /**
           * State whether we have time_dependent boundary conditions. Switched
           * off after finding no more velocity files to suppress attempts to read
           * in new files.
           */
          bool time_dependent;

          /**
           * Scale the age-dependent thermal boundary condition by a scalar factor.
           */
          double scale_factor;

          /**
           * Pointer to an object that reads and processes data we get from
           * gplates files.
           */
          std::unique_ptr<internal::AssimilationLookup<dim> > lookup;

          /**
           * Pointer to an object that reads and processes data we get from
           * gplates files. This saves the previous data time step.
           */
          std::unique_ptr<internal::AssimilationLookup<dim> > old_lookup;

          /**
                  * Handles the update of the age data in lookup. The input
                  * parameter makes sure that both age files (n and n+1) can be
                  * reloaded if the model time step is larger than the velocity file
                  * time step.
          */
          void
          update_data (const bool load_both_files);

          /**
                  * Handles settings and user notification in case the time-dependent
                  * part of the top thermal boundary condition is over.
           */
           void
           end_time_dependence ();

           /**
                  * Create a filename out of the name template.
            */
            std::string
            create_filename (const int timestep) const;

          /**
           * The maximum thickness and temperature of an oceanic plate
           * and a continental plate when time goes to infinity
           */
          double d_cmax;
          double d_omax;
          double Tom;
          double Tcm;

    };
  }
}


#endif

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
#include "assimilation_IT.h"

#include <aspect/global.h>
#include <aspect/utilities.h>

#include <aspect/geometry_model/spherical_shell.h>
#include <aspect/geometry_model/chunk.h>
#include <aspect/geometry_model/sphere.h>

#include <aspect/boundary_temperature/spherical_constant.h>
#include <aspect/boundary_temperature/initial_temperature.h>
#include <aspect/boundary_temperature/interface.h>

#include <boost/lexical_cast.hpp>


namespace aspect
{
  namespace InitialTemperature
  {
    namespace internal
	{
       template <int dim>
       AssimilationLookup<dim>::AssimilationLookup(const unsigned int components,
                                             const double scale_factor)
         :
         components(components),
         data(components),
         maximum_component_value(components),
         scale_factor(scale_factor),
         coordinate_values_are_equidistant(false)
       {}

       template <int dim>
       std::string
       AssimilationLookup<dim>::screen_output() const
       {

         std::ostringstream output;

         output << std::setprecision (3) << std::setw(3) << std::fixed << std::endl
                << " Setting up Age Assimilation - Initial temperature plugin."
				<< std::endl
                << std::endl;

         return output.str();
       }

       template <int dim>
       void
       AssimilationLookup<dim>::load_file(const std::string &filename, const MPI_Comm &comm)
  	   {

  	        // Read data from disk and distribute among processes
  	        std::stringstream in(Utilities::read_and_distribute_file_content(filename,comm));

  	        // Read header lines and table size
  	        while (in.peek() == '#')
  	        {
  	          std::string line;
  	          std::getline(in,line);
  	          std::stringstream linestream(line);
  	          std::string word;
  	          while (linestream >> word)
  	             if (word == "POINTS:")
  	                for (unsigned int i = 0; i < dim; i++)
  	                {
  	                  unsigned int temp_index;
  	                  linestream >> temp_index;

  	                  if (table_points[i] == 0)
  	                    table_points[i] = temp_index;
  	                  else
  	                    AssertThrow (table_points[i] == temp_index,
  	                                 ExcMessage("The file grid must not change over model runtime. "
  	                                            "Either you prescribed a conflicting number of points in "
  	                                            "the input file, or the POINTS comment in your data files "
  	                                            "is changing between following files."));
  	                  }
  	        }
  	         for (unsigned int i = 0; i < dim; i++)
  	         {
  	           AssertThrow(table_points[i] != 0,
  	           ExcMessage("Could not successfully read in the file header of the "
  	                      "ascii data file <" + filename + ">. One header line has to "
  	                      "be of the format: '#POINTS: N1 [N2] [N3]', where N1 and "
  	                      "potentially N2 and N3 have to be the number of data points "
  	                      "in their respective dimension. Check for typos in this line "
  	                      "(e.g. a missing space character)."));
  	         }
  	         // Read column lines if present
  	          unsigned int name_column_index = 0;
  	          double temp_data;

  	          while (true)
  	            {
  	              AssertThrow (name_column_index < 100,
  	                           ExcMessage("The program found more than 100 columns in the first line of the data file. "
  	                                      "This is unlikely intentional. Check your data file and make sure the data can be "
  	                                      "interpreted as floating point numbers. If you do want to read a data file with more "
  	                                      "than 100 columns, please remove this assertion."));

  	              std::string column_name_or_data;
  	              in >> column_name_or_data;
  	              try
  	                {
  	                  // If the data field contains a name this will throw an exception
  	                  temp_data = boost::lexical_cast<double>(column_name_or_data);

  	                  // If there was no exception we have left the line containing names
  	                  // and have read the first data field. Save number of components, and
  	                  // make sure there is no contradiction if the components were already given to
  	                  // the constructor of this class.
  	                  if (components == numbers::invalid_unsigned_int)
  	                    components = name_column_index - dim;
  	                  else if (name_column_index != 0)
  	                    AssertThrow (components == name_column_index,
  	                                 ExcMessage("The number of expected data columns and the "
  	                                            "list of column names at the beginning of the data file "
  	                                            + filename + " do not match. The file should contain "
  	                                            "one column name per column (one for each dimension "
  	                                            "and one per data column)."));

  	                break;
  	                }
  	                catch (const boost::bad_lexical_cast &e)
  	                        {
  	                          // The first dim columns are coordinates and contain no data
  	                          if (name_column_index >= dim)
  	                            {
  	                              // Transform name to lower case to prevent confusion with capital letters
  	                              // Note: only ASCII characters allowed
  	                              std::transform(column_name_or_data.begin(), column_name_or_data.end(), column_name_or_data.begin(), ::tolower);

  	                              AssertThrow(std::find(data_component_names.begin(),data_component_names.end(),column_name_or_data)
  	                                          == data_component_names.end(),
  	                                          ExcMessage("There are multiple fields named " + column_name_or_data +
  	                                                     " in the data file " + filename + ". Please remove duplication to "
  	                                                     "allow for unique association between column and name."));

  	                              data_component_names.push_back(column_name_or_data);
  	                            }
  	                        ++name_column_index;
  	                        }
  	            }
  	        /**
  	         * Create table for the data. This peculiar reinit is necessary, because
  	         * there is no constructor for Table, which takes TableIndices as
  	         * argument.
  	         */
  	        data.resize(components);
  	        maximum_component_value.resize(components,-std::numeric_limits<double>::max());
  	        Table<dim,double> data_table;
  	        data_table.TableBase<dim,double>::reinit(table_points);
  	        std::vector<Table<dim,double> > data_tables(components+dim,data_table);

  	      // Read data lines
  	      unsigned int read_data_entries = 0;
  	      do
  	        {
  	          const unsigned int column_num = read_data_entries%(components+dim);

  	          if (column_num >= dim)
  	            {
  	              temp_data *= scale_factor;
  	              maximum_component_value[column_num-dim] = std::max(maximum_component_value[column_num-dim], temp_data);
  	            }

  	          data_tables[column_num](compute_table_indices(read_data_entries)) = temp_data;

  	          ++read_data_entries;
  	        }
  	      while (in >> temp_data);

  	      AssertThrow(in.eof(),
  	                  ExcMessage ("While reading the data file '" + filename + "' the ascii data "
  	                              "plugin has encountered an error before the end of the file. "
  	                              "Please check for malformed data values (e.g. NaN) or superfluous "
  	                              "lines at the end of the data file."));

  	      const unsigned int n_expected_data_entries = (components + dim) * data_table.n_elements();
  	      AssertThrow(read_data_entries == n_expected_data_entries,
  	                  ExcMessage ("While reading the data file '" + filename + "' the ascii data "
  	                              "plugin has reached the end of the file, but has not found the "
  	                              "expected number of data values considering the spatial dimension, "
  	                              "data columns, and number of lines prescribed by the POINTS header "
  	                              "of the file. Please check the number of data "
  	                              "lines against the POINTS header in the file."));

  	      // In case the data is specified on a grid that is equidistant
  	      // in each coordinate direction, we only need to store
  	      // (besides the data) the number of intervals in each direction and
  	      // the begin- and endpoints of the coordinates.
  	      // In case the grid is not equidistant, we need to keep
  	      // all the coordinates in each direction, which is more costly.
  	      // Here we fill the data structures needed for both cases,
  	      // and check whether the coordinates are equidistant or not.
  	      // We also check the requirement that the coordinates are
  	      // strictly ascending.

  	      // The number of intervals in each direction
  	      std::array<unsigned int,dim> table_intervals;

  	      // Whether or not the grid is equidistant
  	      coordinate_values_are_equidistant = true;

  	      for (unsigned int i = 0; i < dim; i++)
  	        {

  	          table_intervals[i] = table_points[i] - 1;

  	          TableIndices<dim> idx;
  	          double temp_coord = data_tables[i](idx);
  	          double new_temp_coord = 0;

  	          // The minimum coordinates
  	          grid_extent[i].first = temp_coord;

  	          // The first coordinate value
  	          coordinate_values[i].clear();
  	          coordinate_values[i].push_back(temp_coord);

  	          // The grid spacing
  	          double grid_spacing = numbers::signaling_nan<double>();

  	          // Loop over the rest of the coordinate points
  	          for (unsigned int n = 1; n < table_points[i]; n++)
  	            {
  	              idx[i] = n;
  	              new_temp_coord = data_tables[i](idx);
  	              AssertThrow(new_temp_coord > temp_coord,
  	                          ExcMessage ("Coordinates in dimension "
  	                                      + Utilities::int_to_string(i)
  	                                      + " are not strictly ascending. "));

  	              // Test whether grid is equidistant
  	              if (n == 1)
  	                grid_spacing = new_temp_coord - temp_coord;
  	              else
  	                {
  	                  const double current_grid_spacing = new_temp_coord - temp_coord;
  	                  // Compare current grid spacing with first grid spacing,
  	                  // taking into account roundoff of the read-in coordinates
  	                  if (std::abs(current_grid_spacing - grid_spacing) > 0.005*(current_grid_spacing+grid_spacing))
  	                    coordinate_values_are_equidistant = false;
  	                }

  	              // Set the coordinate value
  	              coordinate_values[i].push_back(new_temp_coord);

  	              temp_coord = new_temp_coord;
  	            }

  	          // The maximum coordinate
  	          grid_extent[i].second = temp_coord;
  	        }

  	      // For each data component, set up a GridData,
  	      // its type depending on the read-in grid.
  	      for (unsigned int i = 0; i < components; i++)
  	        {
  	          if (coordinate_values_are_equidistant)
  	            data[i]
  	              = std_cxx14::make_unique<Functions::InterpolatedUniformGridData<dim>> (grid_extent,
  	                                                                                     table_intervals,
  	                                                                                     data_tables[dim+i]);
  	          else
  	            data[i]
  	              = std_cxx14::make_unique<Functions::InterpolatedTensorProductGridData<dim>> (coordinate_values,
  	                                                                                           data_tables[dim+i]);
  	        }
  	   }

       template <int dim>
       double
       AssimilationLookup<dim>::get_age(const Point<dim> &position,
                                     const unsigned int component) const
      {
        Assert(component<components, ExcMessage("Invalid component index"));
        return data[component]->value(position);
      }

      template <int dim>
      TableIndices<dim>
      AssimilationLookup<dim>::compute_table_indices(const unsigned int i) const
      {
        TableIndices<dim> idx;
        idx[0] = (i / (components+dim)) % table_points[0];
        if (dim >= 2)
          idx[1] = ((i / (components+dim)) / table_points[0]) % table_points[1];
        if (dim == 3)
          idx[2] = (i / (components+dim)) / (table_points[0] * table_points[1]);

        return idx;
      }
    }

    template <int dim>
    Assimilation<dim>::Assimilation ()
    :
    current_file_number(0),
    first_data_file_model_time(0.0),
    first_data_file_number(0),
	data_file_time_step(0.0),
    time_weight(0.0),
    time_dependent(true),
	old_lookup(),
    lookup()
  {}

    template <int dim>
    void
    Assimilation<dim>::initialize ()
	{

    	AssertThrow ((Plugins::plugin_type_matches<GeometryModel::SphericalShell<dim> >(this->get_geometry_model()) ||
	    	          Plugins::plugin_type_matches<GeometryModel::Chunk<dim> >(this->get_geometry_model()) ||
                      Plugins::plugin_type_matches<GeometryModel::Sphere<dim> >(this->get_geometry_model())),
	    ExcMessage ("This age assimilation plugin currenlty is only used in "
    	            "a spherical shell, chunk, sphere geometry."));

	    lookup = std_cxx14::make_unique<internal::AssimilationLookup<dim>>(1, this->scale_factor);
	    old_lookup = std_cxx14::make_unique<internal::AssimilationLookup<dim>>(1, this->scale_factor);

	    this->get_pcout() << lookup->screen_output();

	    // Set the first age-file number and load the first files
	    current_file_number = first_data_file_number;

	    this->get_pcout() << std::endl << " Loading first age data for the initial temperature calculation: "
	     	    		  << create_filename (current_file_number) << "." << std::endl << std::endl;

        const int next_file_number = (decreasing_file_order) ? current_file_number - 1 : current_file_number + 1;

	    const std::string filename (create_filename (current_file_number));
	    if (Utilities::fexists(filename))
	    	lookup->load_file(filename,this->get_mpi_communicator());
	    else
	    	AssertThrow(false, ExcMessage (std::string("Age data file <")+ filename + "> not found!"));
	    // If the top thermal boundary condition is constant, switch off time_dependence
	    // immediately. If not, also load the second file for interpolation.
	    // This catches the case that many files are present, but the
	    // parameter file requests a single file.
	    if (create_filename (current_file_number) == create_filename (current_file_number+1))
	    {
	    	end_time_dependence ();
	    }
	    else
	    {
	    	const std::string filename (create_filename (next_file_number));
	        this->get_pcout() << std::endl << " Update age data for prescribing temperature "
	                          << filename << "." << std::endl << std::endl;
	        if (Utilities::fexists(filename))
	        {
	        	lookup.swap(old_lookup);
	            lookup->load_file(filename,this->get_mpi_communicator());
	        }
	        else
	        	end_time_dependence ();
	     }
	  }

    template <int dim>
    std::string
    Assimilation<dim>::create_filename (const int timestep) const
    {
    	std::string templ = data_directory+data_file_name;
        const int size = templ.length();
        std::vector<char> buffer(size+10);
        snprintf(buffer.data(), size + 10, templ.c_str(), timestep);
        std::string str_filename (buffer.data());
        return str_filename;
     }

    /* Update the age file */
    template <int dim>
    void
    Assimilation<dim>::update ()
    {
    	const double time_since_start = this->get_time() - first_data_file_model_time;

  	    if (time_dependent && (time_since_start >= 0.0))
  	    {
  	    	const int time_steps_since_start = static_cast<int> (time_since_start / data_file_time_step);
            this->get_pcout() << std::endl << " time steps since start "
  	            	     	  << time_steps_since_start << "." << std::endl << std::endl;

  	        // whether we need to update our data files. This looks so complicated
  	        // because we need to catch increasing and decreasing file orders and all
            // possible first_data_file_model_times and first_data_file_numbers.
  	        const bool need_update = time_steps_since_start > std::abs(current_file_number - first_data_file_number);

  	        if (need_update)
  	        {
  	        	// The last file, which was tried to be loaded was
  	            // number current_file_number +/- 1, because current_file_number
  	            // is the file older than the current model time
  	            const int old_file_number =
  	                  (decreasing_file_order)
  	                  ?
  	                  (current_file_number - 1)
  	                  :
  	                  (current_file_number + 1);

  	            // Calculate new file_number
  	            current_file_number =
  	                  (decreasing_file_order)
  	                  ?
  	                  (first_data_file_number - time_steps_since_start)
  	                  :
  	                  (first_data_file_number + time_steps_since_start);

  	            const bool load_both_files = std::abs(current_file_number - old_file_number) >= 1;   //LS
  	            update_data(load_both_files);

  	           }

  	           time_weight = (time_since_start / data_file_time_step)
  	                         - std::abs(current_file_number - first_data_file_number);

  	           Assert ((0 <= time_weight) && (time_weight <= 1),
  	           ExcMessage ("Error in set_current_time. Time_weight has to be in [0,1]"));
  	    }
        }

    template <int dim>
    void
    Assimilation<dim>::update_data (const bool load_both_files)
    {
    	// If the time step was large enough to move forward more
        // then one data file we need to load both current files
        // to stay accurate in interpolation
        if (load_both_files)
        {
        	const std::string filename (create_filename (current_file_number));
            this->get_pcout() << std::endl << "   Loading Age data assimilation file "
                              << filename << "." << std::endl << std::endl;
            if (Utilities::fexists(filename))
            {
            	lookup.swap(old_lookup);
                lookup->load_file(filename,this->get_mpi_communicator());
            }

         // If loading current_time_step failed, end time dependent part with old_file_number.
         else
            end_time_dependence ();
        }

      // Now load the next data file. This part is the main purpose of this function.
      const int next_file_number =
                (decreasing_file_order) ?
                current_file_number - 1
                :
                current_file_number + 1;

      const std::string filename (create_filename (next_file_number));
      this->get_pcout() << std::endl << "   Loading next age data assimilation file "
                        << filename << "." << std::endl << std::endl;
      if (Utilities::fexists(filename))
        {
          lookup.swap(old_lookup);
          lookup->load_file(filename,this->get_mpi_communicator());
        }

      // If next file does not exist, end time dependent part with current_time_step.
      else
        end_time_dependence ();
    }

    template <int dim>
    void
    Assimilation<dim>::end_time_dependence ()
    {
      // no longer consider the problem time dependent from here on out
      // this cancels all attempts to read files at the next time steps
      time_dependent = false;
      // Give warning if first processor
      this->get_pcout() << std::endl
              << "   Loading new age data file did not succeed." << std::endl
              << "   Assuming constant age-dependent conditions for the rest of the simulation "
              << std::endl << std::endl;
    }

    // Get the current plate age at each surface point
    template <int dim>
    double
    Assimilation<dim>::ascii_age (const Point<dim> &position) const
    {
      const std::array<double,dim> spherical_position =
	        Utilities::Coordinates::cartesian_to_spherical_coordinates(position);

	  Point<dim> internal_position = position;

      //Since we only consider the chunk/sphere,
	  for (unsigned int i = 0; i < dim; i++)
		  internal_position[i] = spherical_position[i];

	  if(!time_dependent)
	  {
		  return lookup->get_age(internal_position, 0);
	  }
	  else
	  {
		  // Get the old plate age at each surface point first
		  // then calculate the real age
		  const double current_age = lookup->get_age(internal_position, 0);
		  const double old_age = old_lookup->get_age(internal_position,0);
		  const double real_age = time_weight * current_age + (1 - time_weight) * old_age;
		  return real_age;
	  }

    }

    template <int dim>
    double
    Assimilation<dim>::initial_temperature (const Point<dim> &position) const
    {

	    /* Note that this object is time-dependent once the prescribed_temperature is used
	     * Otherwise, it exactly means the initial temperature.
	     */

    	/* Get the depth of the point with respect to the reference surface */
  	    const double depth = this->get_geometry_model().depth(position);
  	    const double maxdepth = this->get_geometry_model().maximal_depth();

        /* Get the temperature at the top and bottom boundary of the model */
        const double Ts = this->get_boundary_temperature_manager().minimal_temperature(this->get_fixed_temperature_boundary_indicators());
        const double Tb = this->get_boundary_temperature_manager().maximal_temperature(this->get_fixed_temperature_boundary_indicators());

//  	  const double dt = (Tb-Ts-1300)*1000/maxdepth;

          double kapa = 1e-6;
          double age_ma = 1e6*year_in_seconds;
          double temp = 0.0;
          double temp1 = 0.0;
          double dt = 0.4;
          double magic_number = 1e-7 * this->get_geometry_model().maximal_depth();

        // Get the plate ages at each point
    	const double real_age = ascii_age(position);

        // Here we use plate cooling model of the lithosphere; from Geodynamic book.
        int n_sum=100;
        double sum_OP = 0.0;
        double sum_CP = 0.0;
        for (int i=1; i<=n_sum; i++)
        {
        	sum_OP += (1/i)*(exp((-kapa*i*i*numbers::PI*numbers::PI*real_age*age_ma)/(d_omax*d_omax)))*(sin(i*numbers::PI*depth/d_omax));
            sum_CP += (1/i)*(exp((-kapa*i*i*numbers::PI*numbers::PI*real_age*age_ma)/(d_cmax*d_cmax)))*(sin(i*numbers::PI*depth/d_cmax));
        }
        if (real_age>0 && real_age<300 && depth <= d_omax+magic_number)
        	{
        		//temp = Ts + (Tom-Ts)*(1.0-erfc(depth/(2*sqrt(kapa*current_age))));
        	    //temp = Ts + depth*(Tom-Ts)/d_omax;
        	    temp1 = Ts + (Tom-Ts)*((depth/d_omax) + (2/numbers::PI)*sum_OP);
        	    temp = std::min(Tom, temp1);
        	}
        // Assume the age of the continent is 300 Ma
        else if (real_age>=300 && depth <= d_cmax+magic_number)
        {
        		temp1 = Ts + (Tcm-Ts)*((depth/d_cmax) + (2/numbers::PI)*sum_CP);
        	    temp = std::min(Tcm, temp1);
        }
        else
        	    temp = Ts+1300+dt*depth/1000;

        /*		if(temp < Ts || temp > Tb)
        	            std::cout << "T " << temp << " " << position << std::endl;
        */
        return std::max(Ts,std::min(temp,Tb));

    }

  template <int dim>
  void
  Assimilation<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.enter_subsection ("Initial temperature model");
    {
      prm.enter_subsection("Age assimilation model");
      {
      	prm.declare_entry ("Data directory",
      	                   "$ASPECT_SOURCE_DIR/data/initial-temperature/ascii-data/test/",
      	                   Patterns::DirectoryName (),
      	                   "The name of a directory that contains the age data. This path "
      	                   "may either be absolute (if starting with a '/') or relative to "
      	                   "the current directory. The path may also include the special "
      	                   "text '$ASPECT_SOURCE_DIR' which will be interpreted as the path "
      	                   "in which the ASPECT source files were located when ASPECT was "
      	                   "compiled. This interpretation allows, for example, to reference "
      	                   "files located in the `data/' subdirectory of ASPECT. ");
      	prm.declare_entry ("Data file name", "phi.%d",
      	                   Patterns::Anything (),
      	                   "The file name of the material data. Provide file in format: "
      	                   "assim_age.\\%d.txt where \\%d is any sprintf integer "
      	                   "qualifier, specifying the format of the current file number.");
      	prm.declare_entry ("Data file time step", "1e6",
      	      	           Patterns::Double (0),
      	      	           "Time step between following age files. "
      	      	           "Depending on the setting of the global 'Use years in output instead of seconds' flag "
      	      	           "in the input file, this number is either interpreted as seconds or as years. "
      	      	           "The default is one million, i.e., either one million seconds or one million years.");
      	prm.declare_entry ("Decreasing file order", "false",
      	                   Patterns::Bool (),
      	                   "In some cases the boundary files are not numbered in increasing "
      	                   "but in decreasing order (e.g. 'Ma BP'). If this flag is set to "
      	                   "'True' the plugin will first load the file with the number "
      	                   "'First velocity file number' and decrease the file number during "
      	                   "the model run.");
      	prm.declare_entry ("First data file model time", "0",
      	                   Patterns::Double (0),
      	                   "Time from which on the age file with number 'First data "
      	                   "file number' is used as boundary condition. Previous to this "
      	                   "time, a no-slip boundary condition is assumed. Depending on the setting "
      	                   "of the global 'Use years in output instead of seconds' flag "
      	                   "in the input file, this number is either interpreted as seconds or as years.");
      	prm.declare_entry ("First data file number", "0",
      	                   Patterns::Integer (),
      	                   "Number of the first age file to be loaded when the model time "
      	                   "is larger than 'First data file model time'.");
      	prm.declare_entry ("Scale factor", "1",
      	                   Patterns::Double (),
      	                   "Scalar factor, which is applied to the boundary velocity. "
      	                   "You might want to use this to scale the velocities to a "
      	                   "reference model (e.g. with free-slip boundary) or another "
      	                   "plate reconstruction.");
      	prm.declare_entry ("Maximum oceanic plate thickness", "70000.0",
      	                   Patterns::Double (0),
      	                   "The maximum thickness of an oceanic plate in the half-space cooling model "
      	                   "for when time goes to infinity. Units: m. " );
        prm.declare_entry ("Maximum oceanic plate temperature", "1601.00",
      	                   Patterns::Double (0),
      	                   "The maximum temperature of an oceanic plate in the plate cooling model "
      	                   "for when time goes to infinity. Units: K. " );
        prm.declare_entry ("Maximum continental plate thickness", "150000.0",
      	               	   Patterns::Double (0),
      	               	   "The maximum thickness of an oceanic plate in the half-space cooling model "
      	               	   "for when time goes to infinity. Units: m. " );
      	prm.declare_entry ("Maximum continental plate temperature", "1633.00",
      	               	   Patterns::Double (0),
      	               	   "The maximum temperature of an oceanic plate in the plate cooling model "
      	               	   "for when time goes to infinity. Units: K. " );

      }
      prm.leave_subsection();
    }
    prm.leave_subsection();
  }

  template <int dim>
  void
  Assimilation<dim>::parse_parameters (ParameterHandler &prm)
  {
	  AssertThrow (dim == 3,
  	  ExcMessage ("The ' Age data assimilation' model for the initial "
  	              "temperature is only available for 3d sphere/chunk computations currently."));

  	  prm.enter_subsection ("Initial temperature model");
  	  {
  		  prm.enter_subsection("Age assimilation model");
  	  	  {
  			  data_directory = Utilities::expand_ASPECT_SOURCE_DIR(prm.get ("Data directory"));
  		      data_file_name = prm.get ("Data file name");
  		      data_file_time_step        = prm.get_double ("Data file time step");
  		      decreasing_file_order      = prm.get_bool   ("Decreasing file order");
  		      first_data_file_model_time = prm.get_double ("First data file model time");
  		      first_data_file_number     = prm.get_integer("First data file number");
  		      scale_factor               = prm.get_double ("Scale factor");
              d_omax = prm.get_double ("Maximum oceanic plate thickness");
              d_cmax = prm.get_double ("Maximum continental plate thickness");
              Tom = prm.get_double ("Maximum oceanic plate temperature");
              Tcm = prm.get_double ("Maximum continental plate temperature");
              if (this->convert_output_to_years()==true)
              {
            	  data_file_time_step        *= year_in_seconds;
                  first_data_file_model_time *= year_in_seconds;
              }

  		  prm.leave_subsection ();
  	      }
  	  prm.leave_subsection ();
	  }

  }
 }
}

  // explicit instantiations
  namespace aspect
  {
    namespace InitialTemperature
    {
      ASPECT_REGISTER_INITIAL_TEMPERATURE_MODEL(Assimilation,
                                              "age_assimilation_T",
                                              "Implementation of a model in which the initial temperature "
                                              "is derived from age files the GPlates program")
    }
  }

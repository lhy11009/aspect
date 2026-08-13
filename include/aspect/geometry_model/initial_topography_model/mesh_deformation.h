#ifndef _aspect_geometry_model_initial_topography_model_mesh_deformation_h
#define _aspect_geometry_model_initial_topography_model_mesh_deformation_h

#include <aspect/geometry_model/initial_topography_model/interface.h>
#include <aspect/simulator_access.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    /**
     * An initial topography model that derives the topography from the
     * initial deformation prescribed by the mesh deformation model.
     *
     * During early initialization, the mesh deformation handler may not yet
     * be initialized. In this case, this model returns zero topography.
     */
    template <int dim>
    class MeshDeformation : public Interface<dim>,
      public SimulatorAccess<dim>
    {
      public:
        /**
         * Return the value of the initial topography from assigned initial mesh deformation
         */
        double
        value (const Point<dim-1> &surface_point) const override;

        /**
         * Return the maximum value of the elevation.
         */
        double
        max_topography () const override;
    };
  }
}

#endif

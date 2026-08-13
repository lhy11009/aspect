#include <aspect/geometry_model/initial_topography_model/mesh_deformation.h>
#include <aspect/mesh_deformation/interface.h>
#include <aspect/geometry_model/box.h>


namespace aspect
{
  namespace InitialTopographyModel
  {
    template <int dim>
    double
    MeshDeformation<dim>::value (const Point<dim-1> &surface_point) const
    {
      /*
       * The initial topography model can be queried before the mesh
       * deformation handler has been initialized. Return zero topography
       * until the handler is ready to provide the initial deformation.
       */
      // todo_df
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Getting initial topography from the mesh deformation model currently "
                             "only supports the box geometry model."));

      double topography = 0.0;

      if (this->has_mesh_deformation_handler())
        {
          const auto &mesh_deformation = this->get_mesh_deformation_handler();
          if (mesh_deformation.is_post_initialization())
            topography = mesh_deformation.initial_surface_deformation(surface_point);
        }

      return topography;
    }



    template <int dim>
    double
    MeshDeformation<dim>::max_topography () const
    {
      // todo_df
      AssertThrow(Plugins::plugin_type_matches<GeometryModel::Box<dim>>(this->get_geometry_model()),
                  ExcMessage("Getting max topography from the mesh deformation model currently "
                             "only supports the box geometry model."));

      double maximal_topography = 0.0;

      if (this->has_mesh_deformation_handler())
        {
          const auto &mesh_deformation = this->get_mesh_deformation_handler();
          if (mesh_deformation.is_post_initialization())
            maximal_topography = mesh_deformation.maximal_initial_surface_deformation();
        }

      return maximal_topography;
    }
  }
}

namespace aspect
{
  namespace InitialTopographyModel
  {
    ASPECT_REGISTER_INITIAL_TOPOGRAPHY_MODEL(MeshDeformation,
                                             "mesh deformation",
                                             "An initial topography model that derives "
                                             "the initial surface topography from the "
                                             "initial deformation prescribed by the "
                                             "mesh deformation model.")
  }
}

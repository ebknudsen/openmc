#include "openmc/hex_mesh.h"
#include <algorithm> // for copy, equal, min, min_element
#include <cmath>     // for ceil, sqrt
#include <cstddef>   // for size_t
#include <gsl/gsl-lite.hpp>
#include <string>

#ifdef OPENMC_MPI
#include "mpi.h"
#endif

#include "xtensor/xbuilder.hpp"
#include "xtensor/xeval.hpp"
#include "xtensor/xmath.hpp"
#include "xtensor/xsort.hpp"
#include "xtensor/xtensor.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h> // for fmt

#include "openmc/capi.h"
#include "openmc/constants.h"
#include "openmc/container_util.h"
#include "openmc/error.h"
#include "openmc/file_utils.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/memory.h"
#include "openmc/message_passing.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle_data.h"
#include "openmc/plot.h"
#include "openmc/random_dist.h"
#include "openmc/search.h"
#include "openmc/settings.h"
#include "openmc/tallies/filter.h"
#include "openmc/tallies/tally.h"
#include "openmc/volume_calc.h"
#include "openmc/xml_interface.h"

/*#ifdef LIBMESH*/
/*#include "libmesh/mesh_modification.h"*/
/*#include "libmesh/mesh_tools.h"*/
/*#include "libmesh/numeric_vector.h"*/
/*#endif*/

#ifdef DAGMC
#include "moab/FileOptions.hpp"
#endif

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

#ifdef LIBMESH
const bool LIBMESH_ENABLED = true;
#else
const bool LIBMESH_ENABLED = false;
#endif

namespace model {

std::unordered_map<int32_t, int32_t> mesh_map;
vector<unique_ptr<Mesh>> meshes;

} // namespace model

HexagonalMesh::HexagonalMesh(pugi::xml_node node) : StructuredMesh {node}
{
  // Determine number of dimensions for mesh
  if (!check_for_node(node, "dimension")) {
    fatal_error("Must specify <dimension> on a regular mesh.");
  }

  xt::xtensor<int, 1> shape = get_node_xarray<int>(node, "dimension");
  int n = n_dimension_ = shape.size();
  if (n != 1 && n != 2) {
    fatal_error("Mesh must be one or two dimensions.");
  }
  std::copy(shape.begin(), shape.end(), shape_.begin());

  // Check that dimensions are all greater than zero
  if (xt::any(shape <= 0)) {
    fatal_error("All entries on the <dimension> element for a tally "
                "mesh must be positive.");
  }

  // Check for lower-left coordinates
  if (check_for_node(node, "lower_left")) {
    // Read mesh lower-left corner location
    lower_left_ = get_node_xarray<double>(node, "lower_left");
  } else {
    fatal_error("Must specify <lower_left> on a hexagonal mesh.");
  }

  // Make sure lower_left and dimension match
  if (shape.size() != lower_left_.size()) {
    fatal_error("Number of entries on <lower_left> must be the same "
                "as the number of entries on <dimension>.");
  }

  if (check_for_node(node, "width")) {
    // Make sure one of upper-right or width were specified
    if (check_for_node(node, "upper_right")) {
      fatal_error("Cannot specify both <upper_right> and <width> on a hexgonal mesh.");
    }

    width_ = get_node_xarray<double>(node, "width");

    // Check to ensure width has same dimensions
    auto n = width_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <width> must be the same as "
                  "the number of entries on <lower_left>.");
    }

    // Check for negative widths
    if (xt::any(width_ < 0.0)) {
      fatal_error("Cannot have a negative <width> on a tally hexgonal mesh.");
    }

    // Set width and upper right coordinate
    upper_right_ = xt::eval(lower_left_ + shape * width_);

  } else if (check_for_node(node, "upper_right")) {
    upper_right_ = get_node_xarray<double>(node, "upper_right");

    // Check to ensure width has same dimensions
    auto n = upper_right_.size();
    if (n != lower_left_.size()) {
      fatal_error("Number of entries on <upper_right> must be the "
                  "same as the number of entries on <lower_left>.");
    }

    // Check that upper-right is above lower-left
    if (xt::any(upper_right_ < lower_left_)) {
      fatal_error("The <upper_right> coordinates must be greater than "
                  "the <lower_left> coordinates on a tally hexgonal mesh.");
    }

    // Set width
    width_ = xt::eval((upper_right_ - lower_left_) / shape);
  } else {
    fatal_error("Must specify either <upper_right> or <width> on a hexagonal mesh.");
  }

  // n.b. must check that the number of hexes is odd
  //   or allow a parameter crown/ring
  int max_a_ = (shape[0] - 1) / 2;
  hex_count_ = 1+ 6/2*(a_max+1)*(a_max);

  element_volume_ = width_[1] * width_[0] * width_[0] * sqrt(3);
  // Set material volumes
  volume_frac_ = 1.0 / hex_count_;

  //size of hex is defined as the radius of the circumscribed circle
  size_ = (width_[0]/shape[0])/sqrt(3.0);

}

const std::string Hegagonal::mesh_type = "hexagonal";

std::string HexgonalMesh::get_mesh_type() const
{
  return mesh_type;
}

int HexagonalMesh::get_index_in_direction(const Position& r, int i) const
{
  switch (i){
    case 0:
      return std::ceil( (0.5 * r.x - 1.0/(2*sqrt(3)) * r.y) / this->abs_a );
    case 1:
      return std::ceil( (1.0/sqrt(3) * r.y) / this->abs_a );
    case 2:
      return std::ceil( ( (0.5 * r.x - 1.0/(2*sqrt(3)) * r.y) - (1.0/sqrt(3) * r.y) )/ this->abs_a );
    case 3:
      // z
      return std::ceil((r - lower_left_[i]) / width_[i]);
  }
}

Position HexagonalMesh::get_positition_of_hexindex(HexMeshIndex ijkl) const
{
  //return the cartesian position of a hex indexed by abcz
  Position r;
  r.x = ijkl[0]*basis_[0][0]*width_[0] + ijkl[1]*basis_[1][0]*width_[0];
  r.y = ijkl[0]*basis_[0][1]*width_[0] + ijkl[1]*basis_[1][1]*width_[0];;
  r.z = ijkl[3]*width_[1];
}

StructuredMesh::MeshIndex HexgonalMesh::get_indices(
  Position r, bool& in_mesh) const
{
  //return index of mesh element in hexes
  local_coords(r);

  HexMeshIndex idx;
  idx[0] = get_index_in_direction(r,0);
  idx[1] = get_index_in_direction(r,1);
  idx[2] = -idx[0] -idx[1];
  idx[3] = get_index_in_direction(r, 3);

  //check if either index is out of bounds
  in_mesh = in_bounds();
  return idx;
}

bool HexagonalMesh::in_bounds(HexMeshIndex& ijkl) const
{
  for (auto it=ijkl.begin(); it<ijkl.end(); it++){
    auto elem = *it;
    if (abs(elem)>max_a_) return false;
  }
  return true;
}


Position HexgonalMesh::sample_element(
  const HexMeshIndex& ijk, uint64_t* seed) const
{
  //return random position inside mesh element
  double a_r = uniform_distribution();
  double b_r = uniform_distribution();

  x = this->basis_[0][0] * (a_r + (ijk[0]-1) ) + this->basis_[1][0] * (b_r + (ijk[1]-1) );
  y = this->basis_[0][1] * (a_r + (ijk[0]-1) ) + this->basis_[1][1] * (b_r + (ijk[1]-1) );

  double z = uniform_distribution(-width_[1]/2.0, width_[1]/2.0, seed);

  return origin_ + Position(x, y, z);
}


double HexgonalMesh::find_r_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  //finds the cylinder radius where r+u crosses a plane
  //maybe relevant to hexagonal since it wil tell us approximately which
  //hex-ring we're in
  if ((shell < 0) || (shell > shape_[0]))
    return INFTY;

  // solve r.x^2 + r.y^2 == r0^2
  // x^2 + 2*s*u*x + s^2*u^2 + s^2*v^2+2*s*v*y + y^2 -r0^2 = 0
  // s^2 * (u^2 + v^2) + 2*s*(u*x+v*y) + x^2+y^2-r0^2 = 0

  const double r0 = grid_[0][shell];
  if (r0 == 0.0)
    return INFTY;

  const double denominator = u.x * u.x + u.y * u.y;

  // Direction of flight is in z-direction. Will never intersect r.
  if (std::abs(denominator) < FP_PRECISION)
    return INFTY;

  // inverse of dominator to help the compiler to speed things up
  const double inv_denominator = 1.0 / denominator;

  const double p = (u.x * r.x + u.y * r.y) * inv_denominator;
  double c = r.x * r.x + r.y * r.y - r0 * r0;
  double D = p * p - c * inv_denominator;

  if (D < 0.0)
    return INFTY;

  D = std::sqrt(D);

  // the solution -p - D is always smaller as -p + D : Check this one first
  if (std::abs(c) <= RADIAL_MESH_TOL)
    return INFTY;

  if (-p - D > l)
    return -p - D;
  if (-p + D > l)
    return -p + D;

  return INFTY;
}

StructuredMesh::MeshDistance HexgonalMesh::find_z_crossing(
  const Position& r, const Direction& u, double l, int shell) const
{
  //finds plane-crossing in z - should be the same as for Cylindrical.

  MeshDistance d;
  d.next_index = shell;

  // Direction of flight is within xy-plane. Will never intersect z.
  if (std::abs(u.z) < FP_PRECISION)
    return d;

  d.max_surface = (u.z > 0.0);
  if (d.max_surface && (shell <= shape_[2])) {
    d.next_index += 1;
    d.distance = (grid_[2][shell] - r.z) / u.z;
  } else if (!d.max_surface && (shell > 0)) {
    d.next_index -= 1;
    d.distance = (grid_[2][shell - 1] - r.z) / u.z;
  }
  return d;
}

StructuredMesh::MeshDistance HexgonalMesh::distance_to_grid_boundary(
  const MeshIndex& ijk, int i, const Position& r0, const Direction& u,
  double l) const
{
  // Compute the distance to the element boundary of index i
  // Since MeshIndex is defined to have 3 and only 3 numbers (can we override?) we need to do a bit of gymnastics to
  // convert to a HexMeshIndex
  //
  Position r = r0 - origin_;

  if (i < 4) {
    //get_index in direction i
    int abc_i = get_index_in_direction(r0,i);
    //get_position_of_hex_at_index

    //get_distance to hex boundary of difference
    return std::min(
      MeshDistance(ijk[i] + 1, true, find_hex_crossing(r, u, l, ijk[i])),
      MeshDistance(ijk[i] - 1, false, find_hex_crossing(r, u, l, ijk[i] - 1)));

  } else {
    return find_z_crossing(r, u, l, ijk[i]);
  }
}




int HexagonalMesh::set_grid()
{
  // not necessary to do this - since it is a regular mesh we can do it on the fly
}

MeshDistance find_hex_crossing(const Position& r0, const direction& u,double l, int i) const
{
  // return the distance to the hex-plane in the desired direction
  // i:  {0:a, 1:b, 2:c}
  ijk =

}


std::pair<vector<double>, vector<double>> HexagonalMesh::plot(
  Position plot_ll, Position plot_ur) const
{
  fatal_error("Plot of hexagonal Mesh not implemented");

  // Figure out which axes lie in the plane of the plot.
  array<vector<double>, 2> axis_lines;
  return {axis_lines[0], axis_lines[1]};
}

void HexgonalMesh::to_hdf5(hid_t group) const
{
  hid_t mesh_group = create_group(group, "mesh " + std::to_string(id_));

  write_dataset(mesh_group, "type", "hexgonal");
  write_dataset(mesh_group, "a_grid", grid_[0]);
  write_dataset(mesh_group, "b_grid", grid_[1]);
  write_dataset(mesh_group, "z_grid", grid_[2]);
  write_dataset(mesh_group, "origin", origin_);

  close_group(mesh_group);
}

double HexgonalMesh::volume(const MeshIndex& ijk) const
{
  double a_i = grid_[0][ijk[0] - 1];
  double a_o = grid_[0][ijk[0]];

  double b_i = grid_[1][ijk[1] - 1];
  double b_o = grid_[1][ijk[1]];

  double z_i = grid_[2][ijk[2] - 1];
  double z_o = grid_[2][ijk[2]];

  return 6 * sqrt(3.0) * (a_o-a_i)*(b_o-b_i)*(z_o - z_i);
}

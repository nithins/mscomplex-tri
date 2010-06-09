#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <cpputils.h>
#include <glutils.h>
#include <tri_edge.h>

namespace trimesh
{
  const uint gc_max_cell_dim = tri_cell_complex_t::cc_dim;
  typedef tri_cell_complex_t::cellid_t  cellid_t;
  typedef double                        cell_fn_t;

  typedef std::vector<cellid_t>   cellid_list_t;
  typedef std::vector<cell_fn_t>  cell_fn_list_t;

  typedef glutils::tri_idx_list_t tri_idx_list_t;
  typedef glutils::vertex_list_t  vertex_list_t;

  const cellid_t invalid_cellid = -1;

  enum eGradientDirection
  {
    GRADDIR_DESCENDING,
    GRADDIR_ASCENDING,
    GRADDIR_COUNT,
  };
}

#endif

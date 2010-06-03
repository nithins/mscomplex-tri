#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <cpputils.h>
#include <glutils.h>

namespace grid
{
  const uint gc_grid_dim = 2;

  typedef double                 cell_fn_t;
  typedef uint                   cellid_t;
  typedef glutils::tri_idx_t     tri_idx_t;
  typedef glutils::vertex_t      vertex_t;

  typedef std::vector<cellid_t>  cellid_list_t;
  typedef std::vector<tri_idx_t> tri_idx_list_t;
  typedef std::vector<vertex_t>  vertex_list_t;
  typedef std::vector<cell_fn_t> cell_fn_list_t;

  const cellid_t invalid_cellid = -1;

  enum eGradientDirection
  {
    GRADDIR_DESCENDING,
    GRADDIR_ASCENDING,
    GRADDIR_COUNT,
  };

}

#endif

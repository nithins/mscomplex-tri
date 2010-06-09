#ifndef GRID_DATASET_ENSURE_H_INCLUDED
#define GRID_DATASET_ENSURE_H_INCLUDED

#include <grid_dataset.h>

// bunch of predicates that throw when I suspect something could be logically
// wrong with the state of the dataset.. to be disabled in release builds

//#ifndef NDEBUG
#define USE_ENSURE_PREDICATES
//#endif

namespace trimesh
{
  inline void ensure_cell_paired(const dataset_t *ds,cellid_t c )
  {
    if (ds->m_cell_flags[c] &dataset_t::CELLFLAG_PAIRED == 0)
      throw std::logic_error ("invalid pair requested");
  }

  inline void ensure_valid_trilist_indexes(const tri_idx_list_t & trilist,uint vert_ct)
  {
    for(uint i = 0 ;i < trilist.size(); ++i)
      for(uint j = 0 ; j < trilist[i].size(); ++j)
        if(!(trilist[i][j] < vert_ct))
          throw std::logic_error("invalid triangle index found");
  }

}
#endif

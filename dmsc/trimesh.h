#ifndef TRIMESH_H_INCLUDED
#define TRIMESH_H_INCLUDED

#include <vector>
#include <cpputils.h>
#include <glutils.h>
#include <tri_edge.h>

namespace trimesh
{
  const uint gc_max_cell_dim =    tri_cc_t::cc_dim;
  typedef tri_cc_t::cellid_t      cellid_t;
  typedef double                  cell_fn_t;

  typedef std::vector<cellid_t>   cellid_list_t;
  typedef std::vector<int>        int_list_t;
  typedef std::vector<char>       char_list_t;
  typedef std::vector<cell_fn_t>  cell_fn_list_t;
  typedef std::vector<char>       bool_list_t;

  typedef glutils::tri_idx_list_t tri_idx_list_t;
  typedef glutils::vertex_list_t  vertex_list_t;

  const cellid_t invalid_cellid = -1;

  enum eGDIR  {GDIR_DES,GDIR_ASC,GDIR_CT};

  class dataset_t;
  class mscomplex_t;
  class data_manager_t;

  typedef boost::shared_ptr<dataset_t>            dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t>          mscomplex_ptr_t;
  typedef boost::shared_ptr<data_manager_t>       data_manager_ptr_t;

  typedef boost::shared_ptr<const dataset_t>      dataset_const_ptr_t;
  typedef boost::shared_ptr<const mscomplex_t>    mscomplex_const_ptr_t;
  typedef boost::shared_ptr<const data_manager_t> data_manager_const_ptr_t;

}

#define _FFL            (std::string("\n")+FILEFUNCLINE)

#endif

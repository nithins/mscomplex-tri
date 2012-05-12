#ifndef TRIMESH_H_INCLUDED
#define TRIMESH_H_INCLUDED

#include <tri_edge.h>

namespace trimesh
{
  typedef tri_cc_t::tri_idx_t      tri_idx_t;
  typedef tri_cc_t::tri_idx_list_t tri_idx_list_t;

  const uint gc_max_cell_dim =    tri_cc_t::cc_dim;
  typedef tri_cc_t::cellid_t      cellid_t;
  typedef double                  fn_t;
  typedef std::vector<fn_t>       fn_list_t;

  typedef boost::numeric::ublas::bounded_vector<int,2>  int_pair_t;
  typedef std::vector<int_pair_t>                       int_pair_list_t;


  typedef std::vector<cellid_t>   cellid_list_t;
  typedef std::vector<int>        int_list_t;
  typedef std::vector<char>       char_list_t;
  typedef std::vector<fn_t>       fn_list_t;
  typedef std::vector<char>       bool_list_t;
  typedef std::vector<cellid_t>   mfold_t;

  const cellid_t invalid_cellid = tri_cc_t::INVALID_VALUE;

  enum eGDIR  {DES,ASC,GDIR_CT};

  class dataset_t;
  class mscomplex_t;

  typedef boost::shared_ptr<dataset_t>            dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t>          mscomplex_ptr_t;

  typedef boost::shared_ptr<const dataset_t>      dataset_const_ptr_t;
  typedef boost::shared_ptr<const mscomplex_t>    mscomplex_const_ptr_t;
}

#define _FFL            (std::string("\n")+FILEFUNCLINE)

#endif

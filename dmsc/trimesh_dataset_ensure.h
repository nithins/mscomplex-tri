#ifndef TRIMESH_DATASET_ENSURE_H_INCLUDED
#define TRIMESH_DATASET_ENSURE_H_INCLUDED

#include <fstream>
#include <trimesh_dataset.h>

namespace trimesh
{
  inline int dataset_t::cell_dim(cellid_t c) const
  {return m_tcc.get_cell_dim(c);}

  inline bool dataset_t::is_boundry(cellid_t c) const
  {return m_tcc.is_cell_boundry(c);}

  template<>
  inline uint dataset_t::get_cets<GDIR_DES>(cellid_t c,cellid_t *cets) const
  {return m_tcc.get_cell_facets(c,cets);}

  template<>
  inline uint dataset_t::get_cets<GDIR_ASC>(cellid_t c,cellid_t *cets) const
  {return m_tcc.get_cell_co_facets(c,cets);}

  template<eGDIR dir>
  inline uint dataset_t::get_co_cets(cellid_t c,cellid_t *cets) const
  {return get_cets<(dir == GDIR_DES)?(GDIR_ASC):(GDIR_DES)>(c,cets);}





  inline const cellid_t& dataset_t::max_fct(cellid_t c) const
  {ASSERT(m_cell_mxfct[c] != invalid_cellid);return m_cell_mxfct[c];}

  inline cellid_t& dataset_t::max_fct(cellid_t c)
  {ASSERT(m_cell_mxfct[c] == invalid_cellid);return m_cell_mxfct[c];}

  template <int dim> inline cellid_t dataset_t::max_vert(cellid_t c) const
  {return max_vert<dim-1>(max_fct(c));}

  template <> inline cellid_t dataset_t::max_vert<0>(cellid_t c) const
  {return c;}

  template <> inline cellid_t dataset_t::max_vert<-1>(cellid_t c) const
  {
    switch(cell_dim(c))
    {
      case 2: c = max_fct(c);//break;
      case 1: c = max_fct(c);//break;
      case 0: return c;
    }
    ASSERT(false&&"invalid celldim");
    return c;
  }

  inline cell_fn_t dataset_t::cell_fn(cellid_t c) const
  {return m_vert_fns[max_vert<-1>(c)];}



  template <int dim>
  inline bool dataset_t::compare_cells(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t f1 = max_fct(c1);
    cellid_t f2 = max_fct(c2);

    if(f1 != f2)
      return compare_cells<dim-1>(f1,f2);

    f1 = m_tcc.get_opp_cell(f1,c1);
    f2 = m_tcc.get_opp_cell(f2,c2);

    bool is_bnd_f1 = m_tcc.is_cell_boundry(f1);
    bool is_bnd_f2 = m_tcc.is_cell_boundry(f2);

    if(is_bnd_f1 != is_bnd_f2)
      return (is_bnd_f1);

    return compare_cells<0>(f1,f2);
  }

  template <>
  inline bool dataset_t::compare_cells<0>(const cellid_t & c1, const cellid_t &c2) const
  {
    ASSERT(m_tcc.get_cell_dim(c1) == 0);
    ASSERT(m_tcc.get_cell_dim(c2) == 0);

    cell_fn_t f1 = m_vert_fns[c1];
    cell_fn_t f2 = m_vert_fns[c2];

    if (f1 != f2)
      return f1 < f2;

    return c1 < c2;
  }


  inline const cellid_t& dataset_t::pair(cellid_t c) const
  {ASSERT(m_cell_pairs[c] != invalid_cellid);return m_cell_pairs[c];}

  inline void dataset_t::pair(cellid_t c,cellid_t p)
  {ASSERT(m_cell_pairs[c] == invalid_cellid);m_cell_pairs[c] = p;
   ASSERT(m_cell_pairs[p] == invalid_cellid);m_cell_pairs[p] = c;}

  inline bool dataset_t::is_paired(cellid_t c) const
  {return m_cell_pairs[c] != invalid_cellid;}

  inline bool dataset_t::is_critical(cellid_t c) const
  {return m_cell_pairs[c] == invalid_cellid;}

  inline const cellid_t& dataset_t::owner(cellid_t c) const
  {ASSERT(m_cell_own[c] != invalid_cellid);return m_cell_own[c];}

  inline cellid_t& dataset_t::owner(cellid_t c)
  {ASSERT(m_cell_own[c] == invalid_cellid);return m_cell_own[c];}


  inline void  dataset_t::save_manifolds(std::string &s,mscomplex_ptr_t msc)
  {
    std::fstream fs(s.c_str(),std::ios::in|std::ios::binary);
    ensure(fs.is_open(),"Unable to open file ");
    save_manifolds(fs,msc);
  }
}
#endif

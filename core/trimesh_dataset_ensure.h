#ifndef TRIMESH_DATASET_ENSURE_H_INCLUDED
#define TRIMESH_DATASET_ENSURE_H_INCLUDED

#include <fstream>
#include <stack>

#include <boost/range.hpp>
#include <boost/typeof/typeof.hpp>

#include <trimesh_dataset.h>

namespace trimesh
{
  inline int dataset_t::cell_dim(cellid_t c) const
  {return m_tcc->get_cell_dim(c);}

  inline bool dataset_t::is_boundry(cellid_t c) const
  {return m_tcc->is_cell_boundry(c);}

  template<>
  inline uint dataset_t::get_cets<DES>(cellid_t c,cellid_t *cets) const
  {return m_tcc->get_cell_facets(c,cets);}

  template<>
  inline uint dataset_t::get_cets<ASC>(cellid_t c,cellid_t *cets) const
  {return m_tcc->get_cell_co_facets(c,cets);}

  template<eGDIR dir>
  inline uint dataset_t::get_co_cets(cellid_t c,cellid_t *cets) const
  {return get_cets<(dir == DES)?(ASC):(DES)>(c,cets);}

  template<>
  inline uint dataset_t::get_points<DES>(cellid_t c,cellid_t *pts) const
  {return m_tcc->get_cell_points(c,pts);}

  template<>
  inline uint dataset_t::get_points<ASC>(cellid_t c,cellid_t *tris) const
  {return m_tcc->get_cell_tris(c,tris);}



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

  template<>
  inline fn_t dataset_t::fn<dataset_t::CFI_MAX>(cellid_t c) const
  {return m_vert_fns[max_vert<-1>(c)];}

  template<>
  inline fn_t dataset_t::fn<dataset_t::CFI_AVE>(cellid_t c) const
  {
    fn_t fn = 0;
    cellid_t pts[20];
    int n  = m_tcc->get_cell_points(c,pts);

    for( int i = 0 ; i < n; ++i)
      fn += m_vert_fns[pts[i]];

    return fn/fn_t(n);
  }


  template <int dim>
  inline bool dataset_t::compare_cells(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t f1 = max_fct(c1);
    cellid_t f2 = max_fct(c2);

    if(f1 != f2)
      return compare_cells<dim-1>(f1,f2);

    f1 = m_tcc->get_opp_cell(f1,c1);
    f2 = m_tcc->get_opp_cell(f2,c2);

    bool is_bnd_f1 = m_tcc->is_cell_boundry(f1);
    bool is_bnd_f2 = m_tcc->is_cell_boundry(f2);

    if(is_bnd_f1 != is_bnd_f2)
      return (is_bnd_f1);

    return compare_cells<0>(f1,f2);
  }

  template <>
  inline bool dataset_t::compare_cells<0>(const cellid_t & c1, const cellid_t &c2) const
  {
    ASSERT(m_tcc->get_cell_dim(c1) == 0);
    ASSERT(m_tcc->get_cell_dim(c2) == 0);

    fn_t f1 = m_vert_fns[c1];
    fn_t f2 = m_vert_fns[c2];

    if (f1 != f2)
      return f1 < f2;

    return c1 < c2;
  }

  template <int di,int dj>
  inline bool dataset_t::compare_cells(const cellid_t & c1, const cellid_t &c2) const
  {
    if ( di < dj)
      return ! compare_cells<dj,di>(c2,c1);

    ASSERT(cell_dim(c1) == di && cell_dim(c2) == dj);

    if( di == dj)
      return compare_cells<di>(c1,c2);

    return compare_cells<di-1,dj>(max_fct(c1),c2);
  }


  template <int dim>
  inline bool dataset_t::compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t oc1 = c1,oc2 = c2;

    if(is_paired(c1) && cell_dim(pair(c1)) == dim+1)
      oc1 = max_fct(pair(c1));

    if(is_paired(c2) && cell_dim(pair(c2)) == dim+1)
      oc2 = max_fct(pair(c2));

    return compare_cells<dim>(oc1,oc2);
  }

  template <int di,int dj>
  inline bool dataset_t::compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t oc1 = c1,oc2 = c2;

    if(is_paired(c1) && cell_dim(pair(c1)) == di+1)
      oc1 = max_fct(pair(c1));

    if(is_paired(c2) && cell_dim(pair(c2)) == dj+1)
      oc2 = max_fct(pair(c2));

    return compare_cells<di,dj>(oc1,oc2);
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

  template <eGDIR dir,typename rng_t>
  inline void dataset_t::get_mfold(mfold_t &mfold,rng_t rng)
  {
    std::stack<cellid_t> stk;

    BOOST_AUTO(b,boost::begin(rng));
    BOOST_AUTO(e,boost::end(rng));

    int dim = cell_dim(*b);

    while ( b!= e) stk.push(*b++);

    cellid_t f[20],*fe,*fb;

    while(!stk.empty())
    {
      cellid_t c = stk.top(); stk.pop();

      ASSERT(cell_dim(c) == dim);

      mfold.push_back(c);

      fb = f; fe = f + get_cets<dir>(c,f);

      for (; fb != fe; ++fb )
        if( is_paired(*fb))
        {
          cellid_t p = pair(*fb);
          if(p != c && cell_dim(p) == dim)
            stk.push(p);
        }
    }
  }


}
#endif

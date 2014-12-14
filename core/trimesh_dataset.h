/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef TRIMESH_DATASET_H_INCLUDED
#define TRIMESH_DATASET_H_INCLUDED

#include <trimesh.h>

#include <boost/iterator/counting_iterator.hpp>

namespace trimesh
{
  class dataset_t
  {
  public:
    const fn_list_t    &m_vert_fns;
    cellid_list_t       m_cell_own;
    cellid_list_t       m_cell_pairs;
    cellid_list_t       m_cell_mxfct;

    boost::shared_ptr<tri_cc_t> m_tcc;

  public:

    dataset_t (const fn_list_t &vert_fns,const tri_idx_list_t & trilist);
    dataset_t (const fn_list_t &vert_fns,tri_cc_ptr_t tcc);

    ~dataset_t ();

  public:
    void  work(mscomplex_ptr_t);

  public:
    inline int cell_dim(cellid_t) const;
    inline bool is_boundry(cellid_t) const;

    template<eGDIR dir>
    inline uint get_cets(cellid_t c,cellid_t *cets) const;

    template<eGDIR dir>
    inline uint get_co_cets(cellid_t c,cellid_t *cets) const;

    template<eGDIR dir>
    inline uint get_points(cellid_t c,cellid_t *pts) const;


    inline const cellid_t& max_fct(cellid_t ) const;
    inline cellid_t& max_fct(cellid_t );
    template <int dim> inline cellid_t max_vert(cellid_t c) const;

    enum eCellFnInterpolant {CFI_MAX,CFI_AVE};
    template<eCellFnInterpolant CFI> inline fn_t fn(cellid_t c) const;

    template <int dim>
    inline bool compare_cells(const cellid_t & c1, const cellid_t &c2) const;

    template <int di,int dj>
    inline bool compare_cells(const cellid_t & c1, const cellid_t &c2) const;

    template <int dim>
    inline bool compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const;

    template <int di,int dj>
    inline bool compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const;



    inline const cellid_t& pair(cellid_t ) const;
    inline void pair(cellid_t,cellid_t );
    inline bool is_paired(cellid_t) const;
    inline bool is_critical(cellid_t) const;

    inline const cellid_t& owner(cellid_t ) const;
    inline cellid_t& owner(cellid_t );

    template <eGDIR dir,typename rng_t>
    inline void get_mfold(mfold_t &,rng_t rng);


  };
}

#include <trimesh_dataset_ensure.h>
#endif

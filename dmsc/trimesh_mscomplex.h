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

#ifndef __GRID_MSCOMPLEX_H_INCLUDED_
#define __GRID_MSCOMPLEX_H_INCLUDED_

#include <set>
#include <map>

#include <trimesh.h>

namespace trimesh
{

  typedef std::vector<uint>         critpt_idx_list_t;
  typedef std::vector<cell_fn_t>    cp_fn_list_t;
  typedef n_vector_t<uint,2>        uint_pair_t;
  typedef std::vector<uint_pair_t>  uint_pair_list_t;

  struct critpt_t
  {
    typedef std::multiset<uint>     conn_set_t;
    typedef std::vector<cellid_t>   disc_t;
    typedef std::vector<uint>       disc_contrib_t;

    cellid_t     cellid;
    uint         pair_idx;
    cell_fn_t    fn;
    uchar        index;

    bool isCancelled;
    bool is_paired;
    bool is_boundry;

    critpt_t()
    {
      isCancelled           = false;
      is_paired             = false;
      pair_idx              = -1;
    }

    // list of idx's of cancelled cps that contribute their disc to this cp

    disc_contrib_t contrib[GRADDIR_COUNT];
    disc_t         disc[GRADDIR_COUNT] ;
    conn_set_t     conn[GRADDIR_COUNT];
  };


  class mscomplex_t
  {
  public:

    typedef std::map<cellid_t,uint>  id_cp_map_t;
    typedef std::vector<critpt_t *>  critpt_list_t;

    critpt_list_t m_cps;
    id_cp_map_t   m_id_cp_map;

    void connect_cps(cellid_t c1,cellid_t c2);

    void connect_cps(uint_pair_t p);

    void add_critpt(cellid_t c,uchar i,cell_fn_t f,bool bflg);

    void simplify(uint_pair_list_t &,double simplification_treshold);

    void un_simplify(const uint_pair_list_t &);

    void simplify_un_simplify(double simplification_treshold );

    void add_disc_tracking_seed_cps();

    void clear();

    mscomplex_t(){}

    ~mscomplex_t();

    void write_discs(const std::string &fn_prefix);

    void print_connections(std::ostream & os);
  };

  typedef mscomplex_t::critpt_list_t           critpt_list_t;
  typedef critpt_t::conn_set_t                 conn_set_t;
  typedef critpt_t::disc_t                     critpt_disc_t;
  typedef critpt_t::conn_set_t::iterator       conn_iter_t;
  typedef critpt_t::conn_set_t::const_iterator const_conn_iter_t;
  typedef critpt_t::disc_contrib_t             disc_contrib_t;

  inline void order_pr_by_cp_index(mscomplex_t *msc,uint_pair_t &e)
  {
    if(msc->m_cps[e[0]]->index < msc->m_cps[e[1]]->index)
      std::swap(e[0],e[1]);
  }
}
#endif

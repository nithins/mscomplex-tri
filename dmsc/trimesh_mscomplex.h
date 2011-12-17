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

#ifndef TRIMESH_MSCOMPLEX_H_INCLUDED
#define TRIMESH_MSCOMPLEX_H_INCLUDED

#include <set>
#include <map>

#include <iostream>
#include <fstream>

#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/function.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/functional.hpp>

#include <trimesh.h>


namespace trimesh
{
  typedef std::vector<cell_fn_t>     cp_fn_list_t;
  typedef n_vector_t<int,2>          int_pair_t;
  typedef std::vector<int_pair_t>    int_pair_list_t;

  typedef std::multiset<uint>                 conn_t;
  typedef std::multiset<uint>::iterator       conn_iter_t;
  typedef std::multiset<uint>::const_iterator const_conn_iter_t;
  typedef std::vector<conn_t>                 conn_list_t;

  class mscomplex_t:public boost::enable_shared_from_this<mscomplex_t>
  {
  public:

    cellid_list_t   m_cp_cellid;
    cellid_list_t   m_cp_vertid;
    int_list_t      m_cp_pair_idx;
    char_list_t     m_cp_index;
    bool_list_t     m_cp_is_cancelled;
    bool_list_t     m_cp_is_boundry;
    cell_fn_list_t  m_cp_fn;

    int_pair_list_t m_canc_list;

    conn_list_t   m_conn[GDIR_CT];
    conn_list_t  &m_des_conn;
    conn_list_t  &m_asc_conn;

    mscomplex_t();
    ~mscomplex_t();

    inline int  get_num_critpts() const;

    void resize(int i);
    void set_critpt(int i,cellid_t c,char idx,cell_fn_t f,cellid_t vert_cell,bool is_bndry);

    void connect_cps(int p, int q);
    void dir_connect_cps(int p , int q);
    void pair_cps(int p , int q);

    inline char& index(int i);
    inline const char& index(int i) const;

    inline int& pair_idx(int i);
    inline const int& pair_idx(int i) const;
    inline bool is_paired(int i) const;

    inline char&       is_canceled(int i);
    inline const char& is_canceled(int i) const;

    inline char&       is_boundry(int i);
    inline const char& is_boundry(int i) const;

    inline cellid_t& cellid(int i);
    inline const cellid_t& cellid(int i) const;

    inline cellid_t& vertid(int i);
    inline const cellid_t& vertid(int i) const;

    inline cell_fn_t& fn(int i);
    inline const cell_fn_t& fn(int i) const;

    inline bool is_extrema(int i) const;
    inline bool is_saddle(int i) const;

  public:

    void simplify(double simplification_treshold);
    void un_simplify();

    void invert_for_collection();

    void cancel_pair(int p, int q);
    void uncancel_pair( int p, int q);

    void clear();

    void write_graph(std::ostream & os) const;
    void write_graph(const std::string & fn) const;

    void stow(std::ostream &os,bool purge_data=true);
    void load(std::istream &is);

    inline void stow(const std::string &f,bool purge_data=true)
    {std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);stow(fs,purge_data);}
    void load(const std::string &f)
    {std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);load(fs);}

    inline std::string cp_info (int cp_no) const;
    inline std::string cp_conn (int cp_no) const;

    typedef boost::counting_iterator<int> iterator_t;
    typedef boost::function<bool (int)>   filter_t;
    typedef bool (mscomplex_t::*memb_filter_t)(int) const;
    typedef boost::filter_iterator<filter_t,iterator_t> fiterator_t;
    template <typename it_t>            class cp_id_iterator;
    typedef cp_id_iterator<fiterator_t> cp_id_fiterator;

    inline iterator_t begin() const;
    inline iterator_t end() const;

    inline fiterator_t fbegin(filter_t f) const;
    inline fiterator_t fend(filter_t f) const;

    inline fiterator_t fbegin(memb_filter_t f) const;
    inline fiterator_t fend(memb_filter_t f) const;

    inline cp_id_fiterator cp_id_fbegin(filter_t f) const;
    inline cp_id_fiterator cp_id_fend(filter_t f) const;
  };

  template <typename it_t>
  class mscomplex_t::cp_id_iterator:public std::iterator
      <std::bidirectional_iterator_tag,cellid_t,int,cellid_t,cellid_t>
  {
  public:
    cp_id_iterator(){}

    cp_id_iterator(mscomplex_const_ptr_t msc,it_t i):m_msc(msc),m_i(i){};
    mscomplex_const_ptr_t m_msc;
    it_t m_i;

    inline cp_id_iterator& operator++(){++m_i; return *this;}
    inline cp_id_iterator& operator--(){--m_i; return *this;}
    inline reference operator*() const {return m_msc->cellid(*m_i);}

    inline bool operator== (const cp_id_iterator &rhs) const
    {return (m_i == rhs.m_i);}

    inline bool operator!= (const cp_id_iterator &rhs) const
    {return !(*this == rhs);}
  };
}

#include <trimesh_mscomplex_ensure.h>
#endif

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
#include <boost/enable_shared_from_this.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <trimesh.h>


namespace trimesh
{
  typedef std::multiset<uint>                 conn_t;
  typedef std::multiset<uint>::iterator       conn_iter_t;
  typedef std::multiset<uint>::const_iterator const_conn_iter_t;
  typedef std::vector<conn_t>                 conn_list_t;

  typedef cellid_list_t            mfold_t;
  typedef std::vector<mfold_t>     mfold_list_t;


  class mscomplex_t:public boost::enable_shared_from_this<mscomplex_t>
  {
  public:

    cellid_list_t   m_cp_cellid;
    cellid_list_t   m_cp_vertid;
    int_list_t      m_cp_pair_idx;
    char_list_t     m_cp_index;
    bool_list_t     m_cp_is_boundry;
    fn_list_t       m_cp_fn;

    int_pair_list_t m_canc_list;
    fn_list_t       m_canc_pers;

    conn_list_t   m_conn[GDIR_CT];
    conn_list_t  &m_des_conn;
    conn_list_t  &m_asc_conn;

    mfold_list_t  m_mfolds[GDIR_CT];
    mfold_list_t &m_des_mfolds;
    mfold_list_t &m_asc_mfolds;

    mscomplex_t();
    ~mscomplex_t();

    inline int  get_num_critpts() const;

    void resize(int i);
    void set_critpt(int i,cellid_t c,char idx,fn_t f,cellid_t vert_cell,bool is_bndry);

    void connect_cps(int p, int q);
    void dir_connect_cps(int p , int q);

    inline void pair_cps(int p , int q);
    inline void pair_cps(const int_pair_t &);

    inline const char& index(int i) const;
    inline const int& pair_idx(int i) const;
    inline const char& is_boundry(int i) const;
    inline const cellid_t& cellid(int i) const;
    inline const cellid_t& vertid(int i) const;
    inline const fn_t& fn(int i) const;

    inline bool is_paired(int i) const;
    inline bool is_not_paired(int i) const;

    inline bool is_extrema(int i) const;
    inline bool is_saddle(int i) const;

    template<eGDIR dir>
    inline mfold_t& mfold(int i){return m_mfolds[dir][i];}

  public:

    void simplify(double simplification_treshold);
    void un_simplify();

    void get_mfolds(dataset_ptr_t ds);
    void clear_mfolds();

    void cancel_pair(int p, int q);
    void uncancel_pair( int p, int q);

    void clear();

    void write_graph(std::ostream & os) const;
    void write_graph(const std::string & fn) const;

    void save(std::ostream &os);
    void load(std::istream &is);

    void save_ascii(const std::string &f);

    void save_mfolds(std::ostream &os,dataset_ptr_t ds);
    inline void save_mfolds(const std::string &f,dataset_ptr_t ds)
    {std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);save_mfolds(fs,ds);}

    inline void save(const std::string &f)
    {std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);save(fs);}
    void load(const std::string &f)
    {std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);load(fs);}

    inline std::string cp_info (int cp_no) const;
    inline std::string cp_conn (int cp_no) const;

    typedef boost::counting_iterator<int> iterator_t;
    typedef boost::iterator_range<iterator_t> range_t;

    inline range_t cp_range()
    {return boost::make_iterator_range
          (iterator_t(0),iterator_t(get_num_critpts()));}

  };

}

#include <trimesh_mscomplex_ensure.h>
#endif

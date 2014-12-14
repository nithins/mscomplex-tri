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

  class merge_dag_t;

  class mscomplex_t:public boost::enable_shared_from_this<mscomplex_t>
  {
  public:

    cellid_list_t   m_cp_cellid;
    cellid_list_t   m_cp_vertid;
    int_list_t      m_cp_pair_idx;
    char_list_t     m_cp_index;
    int_list_t      m_cp_cancno;
    bool_list_t     m_cp_is_boundry;
    fn_list_t       m_cp_fn;
    fn_t            m_fmax;
    fn_t            m_fmin;


    int_pair_list_t m_canc_list;
    fn_list_t       m_canc_pers;

    conn_list_t   m_conn[GDIR_CT];
    conn_list_t  &m_des_conn;
    conn_list_t  &m_asc_conn;

    mfold_list_t  m_mfolds[GDIR_CT];
    mfold_list_t &m_des_mfolds;
    mfold_list_t &m_asc_mfolds;

    int m_multires_version;

    boost::shared_ptr<merge_dag_t>
                  m_merge_dag;



    mscomplex_t(std::string fname);
    mscomplex_t();
    ~mscomplex_t();

    inline int  get_num_critpts() const;

    void resize(int i);
    void set_critpt(int i,cellid_t c,char idx,fn_t f,cellid_t vert_cell,bool is_bndry);

    void connect_cps(int p, int q);
    void disconnect_cps(int p, int q);

    inline int index(int i)   const;
    inline int pair_idx(int i) const;
    inline bool is_boundry(int i) const;
    inline cellid_t cellid(int i) const;
    inline cellid_t vertid(int i) const;
    inline fn_t fn(int i) const;

    inline bool is_paired(int i) const;
    inline bool is_not_paired(int i) const;

    inline bool is_extrema(int i) const;
    inline bool is_saddle(int i) const;
    inline bool is_maxima(int i) const;
    inline bool is_minima(int i) const;

    template <int i>
    inline bool is_index_i_cp(int cp) const;


    template<eGDIR dir>
    inline mfold_t& mfold(int i){return m_mfolds[dir][i];}

  public:

    void simplify(double f_tresh, bool is_nrm=false, int req_nmin=0, int req_nmax=0);
    void set_multires_version(int version);
    inline int  get_multires_version() const {return m_multires_version;}
    int  get_multires_version_for_thresh(double t,bool is_nrm=false) const;

    void collect_mfolds(eGDIR dir, int dim, dataset_ptr_t ds);
    void collect_mfolds(dataset_ptr_t ds);
    void get_mfold(eGDIR dir, int cp,cellid_list_t &mfold,int ver=-1);
    void get_contrib(eGDIR dir, int cp,int_list_t &contrib,int ver=-1);



    void cancel_pair(int p,int q);
    void cancel_pair();
    void anticancel_pair();

    void clear();

    inline void save(const std::string &f) const;
    inline void load(const std::string &f);

    void save_bin(std::ostream &os) const;
    void load_bin(std::istream &is);

    template<class Archive>
    void serialize(Archive & ar, const unsigned int /* file_version */);

    inline std::string cp_info (int cp_no) const;
    inline std::string cp_conn (int cp_no) const;

    typedef boost::counting_iterator<int> iterator_t;
    typedef boost::iterator_range<iterator_t> range_t;

    inline range_t cp_range() const
    {return boost::make_iterator_range
          (iterator_t(0),iterator_t(get_num_critpts()));}

    void save_conn_graph_vtk(const std::string & fn, tri_cc_geom_ptr_t tcc) const;

  };

}

#include <trimesh_mscomplex_ensure.h>
#endif

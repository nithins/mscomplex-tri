#include <cmath>
#include <queue>

#include <logutil.h>

#include <trimesh_mscomplex.h>
#include <trimesh_mscomplex_ensure.h>
#include <limits>

using namespace std;

namespace trimesh
{
  void mark_cancel_pair(mscomplex_t *msc,uint_pair_t e)
  {
    critpt_t * cp1 = msc->m_cps[e[0]];
    critpt_t * cp2 = msc->m_cps[e[1]];

    cp1->is_paired = true;
    cp2->is_paired = true;

    cp1->pair_idx  = e[1];
    cp2->pair_idx  = e[0];

  }

  void mscomplex_t::add_critpt(cellid_t c,uchar i,cell_fn_t f,bool bflg,uint vert_idx)
  {
    critpt_t * cp  = new critpt_t;
    cp->cellid     = c;
    cp->index      = i;
    cp->fn         = f;
    cp->is_boundry = bflg;
    cp->vert_idx   = vert_idx;
    m_id_cp_map.insert(std::make_pair(c,m_cps.size()));
    m_cps.push_back(cp);
  }

  void mscomplex_t::connect_cps(cellid_t c0,cellid_t c1)
  {
    ensure_cellid_critical(this,c0);
    ensure_cellid_critical(this,c1);

    connect_cps(uint_pair_t(m_id_cp_map[c0],m_id_cp_map[c1]));
  }

  bool is_saddle_extremum_pair(mscomplex_t * msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    ensure_index_one_separation(msc,e);

    if(msc->m_cps[e[0]]->index == gc_max_cell_dim ||
       msc->m_cps[e[1]]->index == 0)
      return true;

    return false;

  }

  void mscomplex_t::connect_cps(uint_pair_t e)
  {
    order_pr_by_cp_index(this,e);

    ensure_ordered_index_one_separation(this,e);

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      if(m_cps[e[dir]]->conn[dir].count(e[dir^1]) == 2)
        return;

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      m_cps[e[dir]]->conn[dir].insert(e[dir^1]);

  }

  void cancelPairs ( mscomplex_t *msc,uint_pair_t e,
                     uint_pair_list_t * new_edges = NULL)
  {

    order_pr_by_cp_index(msc,e);

    ensure_single_connectivity(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    conn_iter_t it[GRADDIR_COUNT];

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir].erase(e[dir^1]);

    // cps in lower of u except l
    for(it[0] = cp[0]->conn[0].begin();it[0] != cp[0]->conn[0].end();++it[0])
      for(it[1] = cp[1]->conn[1].begin();it[1] != cp[1]->conn[1].end();++it[1])
      {
        ensure_not_cancelled(msc,*it[0]);
        ensure_not_cancelled(msc,*it[1]);

        msc->connect_cps(uint_pair_t(*it[0],*it[1]));
      }

    if(new_edges)
      for(it[0] = cp[0]->conn[0].begin();it[0] != cp[0]->conn[0].end();++it[0])
        for(it[1] = cp[1]->conn[1].begin();it[1] != cp[1]->conn[1].end();++it[1])
            new_edges->push_back(uint_pair_t(*it[1],*it[0]));

    for(uint dir = 0 ; dir<2;++dir)
     for(conn_iter_t it = cp[dir]->conn[dir].begin();it != cp[dir]->conn[dir].end();++it)
      {
        msc->m_cps[*it]->conn[dir^1].erase(e[dir]);
      }

    for(uint dir = 0 ; dir<2;++dir)
      for(conn_iter_t it = cp[dir]->conn[dir^1].begin();it != cp[dir]->conn[dir^1].end();++it)
      {
        msc->m_cps[*it]->conn[dir].erase(e[dir]);
      }

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir].insert(e[dir^1]);

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->isCancelled = true;

    for(uint dir = 0 ; dir < 2;++dir)
      cp[dir]->conn[dir^1].clear();
  }

  void uncancel_pair( mscomplex_t  *msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    ensure_ordered_connectivity(msc,e);

    ensure_cp_is_cancelled(msc,e[0]);
    ensure_cp_is_cancelled(msc,e[1]);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    conn_iter_t i_it,j_it;

    for(uint dir = 0 ; dir <2 ; ++dir)
    {
      conn_set_t new_conn;

      for(i_it = cp[dir]->conn[dir].begin();i_it != cp[dir]->conn[dir].end() ; ++i_it )
      {
        if(*i_it == e[dir^1]) continue;

        critpt_t *conn_cp = msc->m_cps[*i_it];

        if(!conn_cp->is_paired)
        {
          new_conn.insert(*i_it);
          continue;
        }

        ensure_cp_is_not_cancelled(msc,*i_it);

        ensure_cp_is_paired(msc,*i_it);

        critpt_t *cp_pr = msc->m_cps[conn_cp->pair_idx];

        for(j_it = cp_pr->conn[dir].begin(); j_it != cp_pr->conn[dir].end() ; ++j_it )
        {
          ensure_cp_is_not_paired(msc,*j_it);

          ensure_index_one_separation(msc,uint_pair_t(*j_it,e[dir]));

          if(new_conn.count(*j_it) == 0)
            new_conn.insert(*j_it);
        }
      }

      msc->m_cps[e[dir]]->isCancelled = false;
      cp[dir]->conn[dir].clear();
      cp[dir]->conn[dir].insert(new_conn.begin(),new_conn.end());
    }
  }

  void mscomplex_t::save(std::ostream & os)
  {
    os<<"# Num Cps"<<std::endl;

    os<<m_cps.size()<<std::endl;

    os<<"# SL.No  cpIdx isPaired pair_idx vertNo fn "<<std::endl;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      os<<i<<" ";
      os<<(int)cp->index<<" ";
      os<<(bool)cp->is_paired<<" ";
      os<<(int)cp->pair_idx<<" ";
      os<<cp->vert_idx<<" ";
      os<<(cell_fn_t)cp->fn<<" ";

      os<<std::endl;
    }

    os<<"#slno numDes numAsc connList"<<std::endl;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      os<<(int)i<<" ";

      if(cp->is_paired == true)
      {
        os<<"0 0"<<std::endl;
        continue;
      }
      os<<(int)cp->conn[0].size()<<" ";
      os<<(int)cp->conn[1].size()<<" ";

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        conn_set_t &conn = cp->conn[dir];

        for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
        {
          os<<*it<<" ";
        }
      }
      os<<std::endl;
    }
  }

  void mscomplex_t::save_manifolds(ostream & os,const tri_cc_geom_t &geom)
  {
    os<<"# Ascending and descending manifolds of critical points "<<endl;
    os<<"# ------------------------------------------------------"<<endl;
    os<<"# des manifold of a maximum is a set of tri indices.    "<<endl;
    os<<"#                                                       "<<endl;
    os<<"# des manifold of a saddle is a set of edges ..         "<<endl;
    os<<"# ..edges are written as u v..                          "<<endl;
    os<<"# ..(u,v) are the pair of vertices they separate.       "<<endl;
    os<<"#                                                       "<<endl;
    os<<"# asc manifold of a minimum is a set of vert indices.   "<<endl;
    os<<"# ..verts are written as v n v0 t0 v1 t1 .. vn tn [vn+1]"<<endl;
    os<<"# ..v is the vertex                                     "<<endl;
    os<<"# ..n is number of vi,ti pairs                          "<<endl;
    os<<"# ..vi ti are the verts of the edge (v,vi) and          "<<endl;
    os<<"#   tri ti in the star of v                             "<<endl;
    os<<"# ..[vn+1] is present if v is on boundry.               "<<endl;
    os<<"# ..(v,vn+1) is the final edge of the star              "<<endl;
    os<<"#                                                       "<<endl;
    os<<"# asc manifold of a saddle is a set of edges            "<<endl;
    os<<"# ..edges are written as u v t1 v1 [t2 v2]              "<<endl;
    os<<"# ..(u,v) the pair of vertices they separate            "<<endl;
    os<<"# ..t1,v1 tri adjacent to edge and the vert opposite    "<<endl;
    os<<"# ..[t2,v2] other incident tri and opp vert  .          "<<endl;
    os<<"# ..[t2,v2] will not exist only for edges on boundry.   "<<endl;
    os<<"# ..if [t2,v2] does not exist then cell is critical.    "<<endl;
    os<<"#                                                       "<<endl;
    os<<"# All indices are w.r.t original tri file and 0 based   "<<endl;
    os<<"#-------------------------------------------------------"<<endl;
    os<<"#Num Cps                                                "<<endl;
    os<<"#                                                       "<<endl;
    os<<"#cpIdx isPaired pair_idx vertNo fn desCellCt ascCellCt  "<<endl;
    os<<"#descending cells (newline separated)                   "<<endl;
    os<<"#ascending cells  (newline separated)                   "<<endl;
    os<<"#                                                       "<<endl;
    os<<"#...                                                    "<<endl;

    int tri_offset = geom.get_num_cells_max_dim(1);

    os<<m_cps.size()<<endl;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      critpt_t * cp = m_cps[i];

      cellid_list_t mfold_list[2];

      if(cp->is_paired == false)
      {
        for(int dir = 0 ; dir < GRADDIR_COUNT; ++dir)
        {
          set<cellid_t> cset;

          for(uint j = 0 ; j < cp->contrib[dir].size();++j)
          {
            critpt_t *cp_contrib = m_cps[cp->contrib[dir][j]];

            if(cp_contrib->index != cp->index)
              throw logic_error("contrib and cp must have same idx");

            for(uint k = 0; k < cp_contrib->disc[dir].size(); ++k)
            {
              cellid_t c = cp_contrib->disc[dir][k];

              if(cset.count(c) == 0)
                cset.insert(c);
            }
          }

          mfold_list[dir].insert(mfold_list[dir].begin(),cset.begin(),cset.end());
        }
      }

      switch(cp->index)
      {
      case 0: mfold_list[0].clear();break;
      case 2: mfold_list[1].clear();break;
      };

      os<<i<<" ";
      os<<(int)cp->index<<" ";
      os<<(bool)cp->is_paired<<" ";
      os<<(int)cp->pair_idx<<" ";
      os<<cp->vert_idx<<" ";
      os<<(cell_fn_t)cp->fn<<" ";

      os<<mfold_list[0].size()<<" ";
      os<<mfold_list[1].size()<<" ";
      os<<endl;

      switch(cp->index)
      {
      case 0:
        {
          for(int j = 0 ; j < mfold_list[GRADDIR_ASCENDING].size();++j)
          {
            cellid_t c = mfold_list[GRADDIR_ASCENDING][j];

            cellid_t star[40];

            int star_ct = geom.get_vert_star(c,star);

            os<<c<<" "<<(int)(star_ct/2)<<" ";

            for(int k = 0 ; k < star_ct;++k)
            {
              if(k%2 == 0)
              {
                cellid_t fct[2];

                if(geom.get_cell_facets(star[k],fct) != 2)
                  throw std::runtime_error("unexpected facet ct for edge");

                cellid_t ov =fct[0];

                if(ov == c)
                  ov = fct[1];

                os<<ov<<" ";
              }
              else
              {
                os<<(int)(star[k] - tri_offset)<<" ";
              }
            }
            os<<endl;
          }
          break;
        }
      case 2:
        {
          for(int j = 0 ; j < mfold_list[GRADDIR_DESCENDING].size();++j)
          {
            os<<(mfold_list[GRADDIR_DESCENDING][j]-tri_offset)<<endl;
          }
          break;
        }
      case 1:
        {
          for(int j = 0 ; j < mfold_list[GRADDIR_DESCENDING].size();++j)
          {
            cellid_t fct[10];

            cellid_t c = mfold_list[GRADDIR_DESCENDING][j];

            if(geom.get_cell_facets(c,fct) != 2)
              throw std::runtime_error("incorrect cofacet count");

            os<<fct[0]<<" "<<fct[1]<<endl;
          }

          for(int j = 0 ; j < mfold_list[GRADDIR_ASCENDING].size();++j)
          {
            cellid_t fct[10],cfct[10],pt[10],op_pt[10];

            cellid_t c = mfold_list[GRADDIR_ASCENDING][j];

            if(geom.get_cell_facets(c,fct) != 2)
              throw std::runtime_error("incorrect facet count");

            uint cfct_ct = geom.get_cell_co_facets(c,cfct);

            if(((cfct_ct == 2) || ((cfct_ct == 1) && (c == cp->cellid)))==false)
              throw std::runtime_error("incorrect cofacet count");

            for (int k = 0 ; k < cfct_ct; ++k)
            {
              if(geom.get_cell_points(cfct[k],pt) != 3 )
                throw std::runtime_error("invalid pt count");

              op_pt[k] = pt[0];

              if(op_pt[k] == fct[0] ||op_pt[k] == fct[1])
                op_pt[k] = pt[1];

              if(op_pt[k] == fct[0]||op_pt[k] == fct[1])
                op_pt[k] = pt[2];

              cfct[k] -= tri_offset;
            }

          os<<fct[0]<<" "<<fct[1]<<" "<<cfct[0]<<" "<<op_pt[0]<<" ";

          if(cfct_ct == 2)
            os<<cfct[1]<<" "<<op_pt[1]<<" ";

          os<<endl;
          }
        }
      }

    }
  }

  mscomplex_t::~mscomplex_t()
  {
    clear();
  }

  void shallow_replicate_cp(mscomplex_t &msc, const critpt_t &cp)
  {
    if(msc.m_id_cp_map.count(cp.cellid) != 0)
      throw std::logic_error("this cp is present in msc");

    critpt_t * dest_cp = new critpt_t;

    dest_cp->isCancelled              = cp.isCancelled;
    dest_cp->is_paired                = cp.is_paired;
    dest_cp->cellid                   = cp.cellid;
    dest_cp->index                    = cp.index;
    dest_cp->fn                       = cp.fn;

    msc.m_id_cp_map[dest_cp->cellid]  = msc.m_cps.size();
    msc.m_cps.push_back(dest_cp);
  }

  void mscomplex_t::clear()
  {
    std::for_each(m_cps.begin(),m_cps.end(),&delete_ftor<critpt_t>);
    m_cps.clear();
    m_id_cp_map.clear();
  }

  struct persistence_comparator_t
  {
    mscomplex_t *m_msc;

    persistence_comparator_t(mscomplex_t *m):m_msc(m){}

    bool operator()(const uint_pair_t & p0, const uint_pair_t &p1)
    {
      return cmp_lt(p1,p0);
    }

    bool cmp_lt(uint_pair_t p0, uint_pair_t p1)
    {
      order_pr_by_cp_index(m_msc,p0);
      order_pr_by_cp_index(m_msc,p1);

      int v00 = m_msc->m_cps[p0[0]]->vert_idx;
      int v01 = m_msc->m_cps[p0[1]]->vert_idx;
      int v10 = m_msc->m_cps[p1[0]]->vert_idx;
      int v11 = m_msc->m_cps[p1[1]]->vert_idx;

      cellid_t c00 = m_msc->m_cps[p0[0]]->cellid;
      cellid_t c01 = m_msc->m_cps[p0[1]]->cellid;
      cellid_t c10 = m_msc->m_cps[p1[0]]->cellid;
      cellid_t c11 = m_msc->m_cps[p1[1]]->cellid;

      if( (v00 == v01 ) != (v10 == v11))
        return (v00 == v01 );

      if( (v00 == v01 ) &&(v10 == v11))
      {
        if(v00 == v10)
        {
          if(c00 != c10)
            return c00 < c10;
          else
            return c01 < c11;
        }
        else
        {
          return (v00 < v10);
        }
      }

      cell_fn_t f00 = m_msc->m_cps[p0[0]]->fn;
      cell_fn_t f01 = m_msc->m_cps[p0[1]]->fn;
      cell_fn_t f10 = m_msc->m_cps[p1[0]]->fn;
      cell_fn_t f11 = m_msc->m_cps[p1[1]]->fn;

      cell_fn_t d1 = std::abs(f01-f00);
      cell_fn_t d2 = std::abs(f11-f10);

      if(d1 != d2)
        return d1 < d2;

      if(c00 != c10)
        return c00 < c10;

      return c01 < c11;
    }
  };

  bool is_valid_canc_edge(mscomplex_t *msc,uint_pair_t e )
  {
    order_pr_by_cp_index(msc,e);

    critpt_t * cp[] = {msc->m_cps[e[0]],msc->m_cps[e[1]]};

    for(uint dir = 0 ; dir < 2; ++dir)
      if(cp[dir]->isCancelled)
        return false;

    for(uint dir = 0 ; dir < 2; ++dir)
      if(!cp[dir]->is_boundry && cp[dir^1]->is_boundry)
        return false;

    for(uint dir = 0 ; dir < 2; ++dir)
      if(cp[dir]->conn[dir].count(e[dir^1]) != 1)
        return false;

    return true;
  }

  bool is_epsilon_persistent(mscomplex_t *msc,uint_pair_t e )
  {
    return (msc->m_cps[e[0]]->vert_idx == msc->m_cps[e[1]]->vert_idx);
  }

  void mscomplex_t::simplify(uint_pair_list_t & canc_pairs_list,
                               double simplification_treshold)
  {
    typedef std::priority_queue
        <uint_pair_t,uint_pair_list_t,persistence_comparator_t>
        canc_pair_priq_t;

    persistence_comparator_t comp(this);

    canc_pair_priq_t  canc_pair_priq(comp);



    cell_fn_t max_val = std::numeric_limits<cell_fn_t>::min();
    cell_fn_t min_val = std::numeric_limits<cell_fn_t>::max();

    for(uint i = 0 ;i < m_cps.size();++i)
    {
      critpt_t *cp = m_cps[i];

      max_val = std::max(max_val,m_cps[i]->fn);

      min_val = std::min(min_val,m_cps[i]->fn);

      for(const_conn_iter_t it = cp->conn[0].begin();it != cp->conn[0].end() ;++it)
      {
        if(is_valid_canc_edge(this,uint_pair_t(i,*it)))
          canc_pair_priq.push(uint_pair_t(i,*it));
      }
    }

    double max_persistence = max_val - min_val;

    uint num_cancellations = 0;

    while (canc_pair_priq.size() !=0)
    {
      uint_pair_t pr = canc_pair_priq.top();

      canc_pair_priq.pop();

      if(is_valid_canc_edge(this,pr) == false)
        continue;
      if(is_epsilon_persistent(this,pr) == false)
      {
        double persistence = std::abs(m_cps[pr[0]]->fn-m_cps[pr[1]]->fn)/max_persistence;

        if(persistence >= simplification_treshold)
          break;
      }

      uint_pair_list_t new_edges;

      cancelPairs ( this,pr,&new_edges);
      num_cancellations++;

      for(uint dir = 0 ; dir < 2 ; ++dir)
        m_cps[pr[dir]]->is_paired = true;

      for(uint dir = 0 ; dir < 2 ; ++dir)
        m_cps[pr[dir]]->pair_idx = pr[dir^1];

      canc_pairs_list.push_back(pr);

      for(uint i = 0 ; i < new_edges.size(); i++)
      {
        canc_pair_priq.push(new_edges[i]);
      }
    }
    _LOG_VAR(num_cancellations);
  }

  void mscomplex_t::un_simplify(const uint_pair_list_t &canc_pairs_list)
  {
    for(uint_pair_list_t::const_reverse_iterator it = canc_pairs_list.rbegin();
    it != canc_pairs_list.rend() ; ++it)
      uncancel_pair(this,*it);
  }

  void mscomplex_t::simplify_un_simplify(double simplification_treshold)
  {
    uint_pair_list_t canc_pairs_list;

    simplify(canc_pairs_list,simplification_treshold);

    un_simplify(canc_pairs_list);
  }

  void mscomplex_t::add_disc_tracking_seed_cps()
  {
    for(uint i = 0 ; i < m_cps.size(); ++i)
    {
      if(!m_cps[i]->is_paired)
      {
        for(uint dir = 0 ; dir < 2 ;++dir)
        {
          m_cps[i]->disc[dir].push_back(m_cps[i]->cellid);
          m_cps[i]->contrib[dir].push_back(i);
        }
        continue;
      }

      uint_pair_t e(i,m_cps[i]->pair_idx);

      order_pr_by_cp_index(this,e);

      if(e[0] != i) continue;

      critpt_t * cp[] = {m_cps[e[0]],m_cps[e[1]]};

      ensure_ordered_index_one_separation(this,e);

      ensure_pairing(this,e);

      for(uint dir = 0 ; dir < 2 ;++dir)
      {
        bool need_disc = false;

        for(conn_iter_t it  = cp[dir]->conn[dir].begin();
                        it != cp[dir]->conn[dir].end(); ++it)
        {
          ensure_cp_is_not_paired(this,*it);

          m_cps[*it]->contrib[dir^1].push_back(e[dir^1]);

          need_disc = true;
        }

        if(need_disc)
          cp[dir^1]->disc[dir^1].push_back(cp[dir^1]->cellid);

      }
    }
  }

  void write_disc(const critpt_disc_t *disc,
                  const std::string &prefix,
                  const cellid_t &c )
  {

    std::stringstream ss;
    ss<<prefix;
    ((std::ostream&)ss)<<c;

    std::ofstream os;

    os.open(ss.str().c_str(),std::ios::out|std::ios::ate|std::ios::app|std::ios::binary);


    if(os.is_open() == false)
      throw std::runtime_error("asc/des disc file not writeable");
    os.write((const char*)(const void*)disc->data(),disc->size()*sizeof(cellid_t));

    os.close();
  }

  void mscomplex_t::write_discs(const std::string &fn_prefix)
  {
    critpt_t * cp ;

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      cp = m_cps[i];

      if(cp->is_paired == true)
        continue;

      if(cp->disc[0].size() != 0 )
        write_disc(&cp->disc[0],fn_prefix+"des_",cp->cellid);

      if(cp->disc[1].size() != 0 )
        write_disc(&cp->disc[1],fn_prefix+"asc_",cp->cellid);

    }
  }
}



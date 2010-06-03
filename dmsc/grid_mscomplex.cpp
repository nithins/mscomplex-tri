#include <cmath>
#include <queue>

#include <logutil.h>

#include <grid_mscomplex.h>
#include <grid_mscomplex_ensure.h>
#include <limits>

namespace grid
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

  void mscomplex_t::add_critpt(cellid_t c,uchar i,cell_fn_t f)
  {
    critpt_t * cp = new critpt_t;
    cp->cellid    = c;
    cp->index     = i;
    cp->fn        = f;
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

    if(msc->m_cps[e[0]]->index == gc_grid_dim ||
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

  void print_cp_connections(std::ostream & os,const mscomplex_t &msc,
                            const conn_set_t &conn)
  {

    os<<"{ ";
    for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
    {
      if(msc.m_cps[*it]->is_paired)
        os<<"*";
      os<<msc.m_cps[*it]->cellid;
      os<<", ";
    }
    os<<"}";
  }

  void print_cp_contrib(std::ostream & os,const mscomplex_t &msc,
                            const disc_contrib_t &contrib)
  {
    os<<"{ ";
    for(disc_contrib_t::const_iterator it = contrib.begin(); it != contrib.end(); ++it)
    {
      if(msc.m_cps[*it]->is_paired)
        os<<"*";
      os<<msc.m_cps[*it]->cellid;
      os<<", ";
    }
    os<<"}";
  }



  void mscomplex_t::print_connections(std::ostream & os)
  {
    const char *conn_dir_txt[] = {"des","asc"};

    const char *contrib_dir_txt[] = {"des_contrib","asc_contrib"};

    for(uint i = 0 ; i < m_cps.size();++i)
    {
      if(m_cps[i]->is_paired)
        os<<"*";

      os<<"cellid = "<<m_cps[i]->cellid<<"\n";

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        os<<conn_dir_txt[dir]<<"=";
        print_cp_connections(os,*this,m_cps[i]->conn[dir]);
        os<<'\n';
      }
      os<<'\n';
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

  mscomplex_t * mscomplex_t::merge_up(const mscomplex_t& msc1,const mscomplex_t& msc2)
  {

    // form the intersection rect
    rect_t ixn;

    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
      throw std::logic_error ("rects should intersect for merge");

    if ( ixn.eff_dim() != gc_grid_dim -1)
      throw std::logic_error ("rects must merge along a 1 manifold");

    // TODO: ensure that the union of  rects is not including anything extra

    rect_t r = msc1.m_rect.bounding_box(msc2.m_rect);

    rect_t e = msc1.m_ext_rect.bounding_box(msc2.m_ext_rect);

    mscomplex_t * out_msc = new mscomplex_t(r,e);

    const mscomplex_t* msc_arr[] = {&msc1,&msc2};

    // make a union of the critical points in this
    for (uint i = 0 ; i <2;++i)
    {
      const mscomplex_t * msc = msc_arr[i];

      for (uint j = 0 ; j <msc->m_cps.size();++j)
      {
        const critpt_t *src_cp = msc->m_cps[j];

        // if it is contained or not
        if (i == 1 && (out_msc->m_id_cp_map.count(src_cp->cellid) == 1))
          continue;

        if(src_cp->isCancelled)
          continue;

        shallow_replicate_cp(*out_msc,*src_cp);

      }
    }

    for (uint i = 0 ; i <2;++i)
    {
      const mscomplex_t * msc = msc_arr[i];

      // copy over connectivity information
      for (uint j = 0 ; j <msc->m_cps.size();++j)
      {
        const critpt_t *src_cp = msc->m_cps[j];

        if(src_cp->isCancelled)
          continue;

        critpt_t *dest_cp = out_msc->m_cps[out_msc->m_id_cp_map[src_cp->cellid]];

        if(src_cp->is_paired)
        {
          critpt_t *src_pair_cp = msc->m_cps[src_cp->pair_idx];

          dest_cp->pair_idx = out_msc->m_id_cp_map[src_pair_cp->cellid];
        }

        bool is_src_cmn_bndry = (ixn.contains(src_cp->cellid) && i == 1);

        for (uint j = 0 ; j < 2; ++j)
        {
          for (const_conn_iter_t it = src_cp->conn[j].begin();
          it != src_cp->conn[j].end();++it)
          {
            const critpt_t *conn_cp = msc->m_cps[*it];

            // common boundry connections would have been found along the boundry
            if( is_src_cmn_bndry && ixn.contains(conn_cp->cellid))
              continue;

            if (conn_cp->isCancelled)
              continue;

            dest_cp->conn[j].insert (out_msc->m_id_cp_map[conn_cp->cellid]);
          }
        }
      }
    }

    // carry out the cancellation
    for(cell_coord_t y = ixn[1][0]; y <= ixn[1][1];++y)
    {
      for(cell_coord_t x = ixn[0][0]; x <= ixn[0][1];++x)
      {
        cellid_t c(x,y);

        if(out_msc->m_id_cp_map.count(c) != 1)
          throw std::logic_error("missing common bndry cp");

        u_int src_idx = out_msc->m_id_cp_map[c];

        critpt_t *src_cp = out_msc->m_cps[src_idx];

        if(src_cp->isCancelled || !src_cp->is_paired)
          continue;

        u_int pair_idx = src_cp->pair_idx;

        cellid_t p = out_msc->m_cps[pair_idx]->cellid;

        if(!out_msc->m_rect.isInInterior(c)&& !out_msc->m_ext_rect.isOnBoundry(c))
          continue;

        if(!out_msc->m_rect.isInInterior(p)&& !out_msc->m_ext_rect.isOnBoundry(p))
          continue;

        cancelPairs(out_msc,uint_pair_t(src_idx,pair_idx));
      }
    }

    return out_msc;
  }

  void mscomplex_t::merge_down(mscomplex_t& msc1,mscomplex_t& msc2)
  {
    rect_t ixn;

    if (!msc2.m_rect.intersection (msc1.m_rect,ixn))
      throw std::logic_error ("rects should intersect for merge");

    if ( ixn.eff_dim() != gc_grid_dim -1)
      throw std::logic_error ("rects must merge along a 1 manifold");

    for(cell_coord_t y = ixn[1][1]; y >= ixn[1][0];--y)
    {
      for(cell_coord_t x = ixn[0][1]; x >= ixn[0][0];--x)
      {

        cellid_t c(x,y);

        if(this->m_id_cp_map.count(c) != 1)
          throw std::logic_error("missing common bndry cp");

        u_int src_idx = this->m_id_cp_map[c];

        critpt_t *src_cp = this->m_cps[src_idx];

        if(!src_cp->isCancelled )
          continue;

        u_int pair_idx = src_cp->pair_idx;

        cellid_t p = this->m_cps[pair_idx]->cellid;

        if(!this->m_rect.isInInterior(c)&& !this->m_ext_rect.isOnBoundry(c))
          continue;

        if(!this->m_rect.isInInterior(p)&& !this->m_ext_rect.isOnBoundry(p))
          continue;

        uncancel_pair(this,uint_pair_t(src_idx,pair_idx));
      }
    }

    // identify and copy the results to msc1 and msc2

    mscomplex_t* msc_arr[] = {&msc1,&msc2};

    for (uint i = 0 ; i <2;++i)
    {
      mscomplex_t * msc = msc_arr[i];

      // adjust connections for uncancelled cps in msc
      for(uint j = 0 ; j < m_cps.size();++j)
      {
        critpt_t * src_cp = m_cps[j];

        if(src_cp->isCancelled)
          throw std::logic_error("all cps ought to be uncancelled by now");

        if(!src_cp->is_paired)
          continue;

        critpt_t * src_pair_cp = m_cps[src_cp->pair_idx];

        bool src_in_msc      = (msc->m_id_cp_map.count(src_cp->cellid) != 0);
        bool src_pair_in_msc = (msc->m_id_cp_map.count(src_pair_cp->cellid) != 0);

        if(!src_in_msc && !src_pair_in_msc)
          continue;

        if(!src_in_msc)
        {
          shallow_replicate_cp(*msc,*src_cp);
        }

        if(!src_pair_in_msc)
        {
          shallow_replicate_cp(*msc,*src_pair_cp);
        }

        uint dest_cp_idx = msc->m_id_cp_map[src_cp->cellid];

        critpt_t *dest_cp = msc->m_cps[dest_cp_idx];

        if(!src_in_msc || !src_pair_in_msc || !dest_cp->is_paired)
        {
          uint dest_pair_cp_idx  = msc->m_id_cp_map[src_pair_cp->cellid];
          critpt_t *dest_pair_cp = msc->m_cps[dest_pair_cp_idx];

          dest_cp->is_paired      = true;
          dest_pair_cp->is_paired = true;

          dest_cp->pair_idx      = dest_pair_cp_idx;
          dest_pair_cp->pair_idx = dest_cp_idx;
        }

        for(uint k = 0 ; k < 2;++k)
        {
          dest_cp->conn[k].clear();

          for(conn_iter_t it = src_cp->conn[k].begin(); it!=src_cp->conn[k].end();++it)
          {
            critpt_t *src_conn_cp = m_cps[*it];

            if(src_conn_cp->is_paired == true)
              throw std::logic_error("only non cancellable cps must be remaining");

            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
            {
              shallow_replicate_cp(*msc,*src_conn_cp);
            }

            dest_cp->conn[k].insert(msc->m_id_cp_map[src_conn_cp->cellid]);
          }// end it
        }// end k
      }// end j

      // adjust connections for non uncancelled cps in msc
      for(uint j = 0 ; j < m_cps.size();++j)
      {
        critpt_t * src_cp = m_cps[j];

        if(src_cp->is_paired)
          continue;

        if(msc->m_id_cp_map.count(src_cp->cellid) != 1)
          continue;

        critpt_t *dest_cp = msc->m_cps[msc->m_id_cp_map[src_cp->cellid]];

        for(uint k = 0 ; k < 2;++k)
        {
          dest_cp->conn[k].clear();

          for(conn_iter_t it = src_cp->conn[k].begin(); it!=src_cp->conn[k].end();++it)
          {
            critpt_t *src_conn_cp = m_cps[*it];

            if(src_conn_cp->is_paired == true)
            {
              _LOG("appears that "<<src_conn_cp->cellid<<"is still connected to"<<
                   src_cp->cellid);
              //throw std::logic_error("only non cancellable cps must be remaining 1");
            }

            if(msc->m_id_cp_map.count(src_conn_cp->cellid) == 0)
              continue;

            dest_cp->conn[k].insert(msc->m_id_cp_map[src_conn_cp->cellid]);
          }// end it
        }// end k
      }// end j
    }//end i
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

    bool operator()(const uint_pair_t & p1, const uint_pair_t &p2)
    {
      cell_fn_t f1 = m_msc->m_cps[p1[0]]->fn;
      cell_fn_t f2 = m_msc->m_cps[p1[1]]->fn;
      cell_fn_t f3 = m_msc->m_cps[p2[0]]->fn;
      cell_fn_t f4 = m_msc->m_cps[p2[1]]->fn;

      cell_fn_t d1 = std::abs(f2-f1);
      cell_fn_t d2 = std::abs(f4-f3);

      if(d1 != d2)
        return d1>d2;

      cellid_t c1 = m_msc->m_cps[p1[0]]->cellid;
      cellid_t c2 = m_msc->m_cps[p1[1]]->cellid;

      cellid_t c3 = m_msc->m_cps[p2[0]]->cellid;
      cellid_t c4 = m_msc->m_cps[p2[1]]->cellid;

      d1 = (c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]);
      d2 = (c3[0]-c4[0])*(c3[0]-c4[0]) + (c3[1]-c4[1])*(c3[1]-c4[1]);

      if(d1 != d2)
        return d1>d2;

      if(c1 > c2)
        std::swap(c1,c2);

      if(c3 > c4)
        std::swap(c3,c4);

      if(c1 != c3)
        return c1 > c3;

      return c2 > c4;
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
      if(!msc->m_rect.isOnBoundry(cp[dir]->cellid) &&
         msc->m_rect.isOnBoundry(cp[dir^1]->cellid))
        return false;


    for(uint dir = 0 ; dir < 2; ++dir)
      if(cp[dir]->conn[dir].count(e[dir^1]) != 1)
        return false;

    return true;
  }

  void mscomplex_t::simplify(uint_pair_list_t & canc_pairs_list,
                               double simplification_treshold)
  {
    typedef std::priority_queue
        <uint_pair_t,uint_pair_list_t,persistence_comparator_t>
        canc_pair_priq_t;

    persistence_comparator_t comp(this);

    canc_pair_priq_t  canc_pair_priq(comp);

    cell_fn_t max_persistence = 0.0;

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

    max_persistence = max_val - min_val;

    _LOG_VAR(max_persistence);

    uint num_cancellations = 0;

    while (canc_pair_priq.size() !=0)
    {
      uint_pair_t pr = canc_pair_priq.top();

      canc_pair_priq.pop();

      cell_fn_t persistence = std::abs(m_cps[pr[0]]->fn-m_cps[pr[1]]->fn);

      if((double)persistence/(double)max_persistence > simplification_treshold)
        break;

      if(is_valid_canc_edge(this,pr) == false)
        continue;

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
      if(!m_ext_rect.contains(m_cps[i]->cellid))
        continue;

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

#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/binary_object.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

//#include <grid_dataset.h>

namespace boost
{
  namespace serialization
  {
    template<class Archive>
    void serialize(Archive & ar, grid::rect_point_t & r, const unsigned int )
    {
      typedef boost::array<grid::rect_point_t::value_type,grid::rect_point_t::static_size> base_t;

      ar & boost::serialization::base_object<base_t>(r);
    }

    template<class Archive>
    void serialize(Archive & ar, grid::rect_range_t & r, const unsigned int )
    {
      typedef boost::array<grid::rect_range_t::value_type,grid::rect_range_t::static_size> base_t;

      ar & boost::serialization::base_object<base_t>(r);
    }

    template<class Archive>
    void serialize(Archive & ar, grid::rect_t & r, const unsigned int )
    {
      typedef boost::array<grid::rect_t::value_type,grid::rect_t::static_size> base_t;

      ar & boost::serialization::base_object<base_t>(r);
    }

    template<class Archive>
    void serialize(Archive & ar, grid::critpt_t & c, const unsigned int )
    {
      ar & c.cellid;
      ar & c.conn[0];
      ar & c.conn[1];
      ar & c.is_paired;
      ar & c.isCancelled;
      ar & c.fn;
      ar & c.pair_idx;
    }


    template<class Archive>
    void serialize(Archive & ar, grid::mscomplex_t & g, const unsigned int )
    {
      ar & g.m_rect;
      ar & g.m_ext_rect;
      ar & g.m_id_cp_map;
      ar & g.m_cps;
    }

    //    template<class Archive>
    //    void serialize(Archive & ar, GridDataset & ds, const unsigned int )
    //    {
    //       ar & ds.m_rect;
    //       ar & ds.m_ext_rect;
    //
    //       GridDataset::rect_size_t ext_sz = ds.m_ext_rect.size();
    //       uint num_data_items = (ext_sz[0]+1)*(ext_sz[1]+1);
    //
    //       if(Archive::is_loading::value)
    //         ds.init(NULL);
    //
    //       ar & make_binary_object(ds.(*m_cell_flags).data(),num_data_items*sizeof(GridDataset::cell_flag_t));
    //       ar & make_binary_object(ds.m_cell_own.data(),num_data_items*sizeof(GridDataset::cellid_t));
    //       ar & make_binary_object(ds.m_cell_pairs.data(),num_data_items*sizeof(GridDataset::cellid_t));
    //    }
  }
}

//// without the explicit instantiations below, the program will
//// fail to link for lack of instantiantiation of the above function
//// The impls are visible only in this file to save compilation time..
//
//template void boost::serialization::serialize<boost::archive::text_iarchive>(
//    boost::archive::text_iarchive & ar,
//    GridDataset & g,
//    const unsigned int file_version
//);
//
//template void boost::serialization::serialize<boost::archive::text_oarchive>(
//    boost::archive::text_oarchive & ar,
//    GridDataset & g,
//    const unsigned int file_version
//);


// without the explicit instantiations below, the program will
// fail to link for lack of instantiantiation of the above function
// The impls are visible only in this file to save compilation time..

template void boost::serialization::serialize<boost::archive::binary_iarchive>(
    boost::archive::binary_iarchive & ar,
    grid::mscomplex_t & g,
    const unsigned int file_version
    );
template void boost::serialization::serialize<boost::archive::binary_oarchive>(
    boost::archive::binary_oarchive & ar,
    grid::mscomplex_t & g,
    const unsigned int file_version
    );



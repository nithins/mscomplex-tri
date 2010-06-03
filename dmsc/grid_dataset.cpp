#include <grid_dataset.h>

#include <vector>
#include <queue>

#include <timer.h>

#include <QFile>
#include <logutil.h>

#include <grid_mscomplex.h>

namespace grid
{
  cellid_t get_cp_cellid(mscomplex_t *msgraph,uint idx)
  {
    return msgraph->m_cps[idx]->cellid;
  }

  static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
  {
    &dataset_t::getCellFacets,
    &dataset_t::getCellCofacets
  };

  inline uint   dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
  {
    inc[0] = cellid_t (c[0]  ,c[1]+1);
    inc[1] = cellid_t (c[0]  ,c[1]-1);
    inc[2] = cellid_t (c[0]-1,c[1]);
    inc[3] = cellid_t (c[0]+1,c[1]);
    return 4;
  }

  void compute_disc_bfs
      (dataset_t *dataset,
       critpt_disc_t *disc,
       cellid_t start_cellId,
       eGradientDirection gradient_dir
       )
  {
    typedef cellid_t id_type;

    std::queue<id_type> cell_queue;

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      id_type top_cell = cell_queue.front();

      cell_queue.pop();

      disc->push_back(top_cell);

      id_type cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( !dataset->isCellCritical ( cets[i] ) )
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            id_type next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  int dataset_t::postMergeFillDiscs(mscomplex_t *msgraph)
  {
    msgraph->add_disc_tracking_seed_cps();

    for(uint i = 0 ; i < msgraph->m_cps.size() ; ++i)
    {
      critpt_t * cp = msgraph->m_cps[i];

//      if(cp->index != 1) continue;

      for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      {
        if(cp->disc[dir].size() == 1)
        {
          cp->disc[dir].clear();
          compute_disc_bfs(this,&cp->disc[dir],cp->cellid,(eGradientDirection)dir);
        }
      }
    }

    return 0;
  }

  void connectCps (mscomplex_t *msgraph,cellid_t c1,cellid_t c2)
  {
    if(msgraph->m_rect.contains(c1) && msgraph->m_rect.contains(c2))
      msgraph->connect_cps(c1,c2);
  }

  void connectCps (mscomplex_t *msgraph,uint i1,uint i2)
  {
    connectCps(msgraph,msgraph->m_cps[i1]->cellid,msgraph->m_cps[i2]->cellid);
  }

  inline bool lowestPairableCoFacet
      (dataset_t *dataset,cellid_t cellId,cellid_t& pairid
       )
  {
    cellid_t cofacets[20];
    bool    cofacet_usable[20];

    uint cofacet_count = dataset->getCellCofacets ( cellId,cofacets );

    bool isTrueBoundryCell = dataset->isTrueBoundryCell ( cellId ) ;

    // for each co facet
    for ( uint i = 0 ; i < cofacet_count ; i++ )
    {
      cellid_t facets[20];
      uint facet_count = dataset->getCellFacets ( cofacets[i],facets );

      cofacet_usable[i] = true;

      if ( isTrueBoundryCell &&
           !dataset->isTrueBoundryCell ( cofacets[i] ) )
      {
        cofacet_usable[i] = false;
        continue;
      }

      for ( uint j = 0 ; j < facet_count ; j++ )
      {
        if ( dataset->compareCells ( cellId,facets[j] ))
        {
          cofacet_usable[i] = false;
          break;
        }
      }
    }

    bool pairid_usable = false;

    for ( uint i =0 ; i < cofacet_count;i++ )
    {
      if ( cofacet_usable[i] == false )
        continue;

      if(pairid_usable == false)
      {
        pairid_usable = true;
        pairid = cofacets[i];
        continue;
      }

      if ( dataset->compareCells ( cofacets[i],pairid ) )
        pairid = cofacets[i];

    }
    return pairid_usable;
  }


  void track_gradient_tree_bfs
      (dataset_t *dataset,cellid_t start_cellId,eGradientDirection gradient_dir)
  {
    std::queue<cellid_t> cell_queue;

    // mark here that that cellid has no parent.

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      cellid_t top_cell = cell_queue.front();

      cell_queue.pop();

      (*dataset->m_cell_own)(top_cell) = start_cellId;

      cellid_t      cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( dataset->isCellCritical ( cets[i] ) )
        {
          //        connectCps(msgraph,start_cellId,cets[i]);
        }
        else
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            cellid_t next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              (*dataset->m_cell_own)(cets[i]) = start_cellId;

              // mark here that the parent of next cell is top_cell
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e) :
      m_rect (r),m_ext_rect (e),m_ptcomp(this)
  {

    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
    //
    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_own     = NULL;
  }

  dataset_t::dataset_t () :
      m_ptcomp(this)
  {
    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_own     = NULL;
  }

  dataset_t::~dataset_t ()
  {
    clear();

    clear_fnref();
  }

  void dataset_t::init()
  {
    rect_size_t   s = m_ext_rect.size();

    m_cell_flags = new cellflag_array_t( (boost::extents[1+s[0]][1+s[1]]));
    m_cell_pairs = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));
    m_cell_own   = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));

    for (int y = 0 ; y<=s[1];++y)
      for (int x = 0 ; x<=s[0];++x)
        (*m_cell_flags)[x][y] = CELLFLAG_UNKNOWN;

    rect_point_t bl = m_ext_rect.lower_corner();

    (*m_cell_flags).reindex (bl);
    (*m_cell_pairs).reindex (bl);
    (*m_cell_own).reindex (bl);
  }

  void  dataset_t::clear()
  {
    if(m_cell_flags != NULL)
      delete m_cell_flags;

    if(m_cell_pairs != NULL)
      delete m_cell_pairs;

    if(m_cell_own != NULL)
      delete m_cell_own;

    m_critical_cells.clear();

    m_saddle_incidence_idx.clear();
    m_saddle_incidence_idx_offset.clear();

    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_own     = NULL;
  }

  void dataset_t::init_fnref(cell_fn_t * pData)
  {
    rect_size_t   s = m_ext_rect.size();

    if(pData != NULL)
      m_vert_fns_ref =
          new varray_ref_t(pData,boost::extents[1+s[0]/2][1+s[1]/2],
                           boost::fortran_storage_order());

    rect_point_t bl = m_ext_rect.lower_corner();

    if(pData != NULL)
      (*m_vert_fns_ref).reindex (bl/2);

  }

  void dataset_t::clear_fnref()
  {
    if(m_vert_fns_ref != NULL)
      delete m_vert_fns_ref;

    m_vert_fns_ref = NULL;
  }

  cellid_t   dataset_t::getCellPairId (cellid_t c) const
  {
    if ((*m_cell_flags) (c) &CELLFLAG_PAIRED == 0)
      throw std::logic_error ("invalid pair requested");

    return (*m_cell_pairs) (c);
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    if(getCellDim(c1) == 0)
      return ptLt(c1,c2);

    cellid_t pts1[20];
    cellid_t pts2[20];

    uint pts1_ct = getCellPoints ( c1,pts1);
    uint pts2_ct = getCellPoints ( c2,pts2);

    std::sort ( pts1,pts1+pts1_ct,m_ptcomp );
    std::sort ( pts2,pts2+pts2_ct,m_ptcomp);

    return std::lexicographical_compare
        ( pts1,pts1+pts1_ct,pts2,pts2+pts2_ct,
          m_ptcomp );
  }

  cell_fn_t dataset_t::get_cell_fn (cellid_t c) const
  {

    if(m_vert_fns_ref == NULL)
      return 0.0;

    cell_fn_t  fn = 0.0;

    cellid_t pts[20];

    uint pts_ct = getCellPoints (c,pts);

    for (int j = 0 ; j <pts_ct ;++j)
      fn += (*m_vert_fns_ref) (pts[j]/2);

    fn /= pts_ct;

    return fn;
  }

  void dataset_t::set_cell_fn (cellid_t c,cell_fn_t f)
  {
    if (getCellDim (c) != 0)
      throw std::logic_error ("values only for vertices are specified");

    c[0] /=2;

    c[1] /=2;

    (*m_vert_fns_ref) (c) = f;
  }

  uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
  {
    switch (getCellDim (c))
    {
    case 0:
      p[0] = c;
      return 1;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        p[0] = cellid_t (c[0]+d0,c[1]+d1);
        p[1] = cellid_t (c[0]-d0,c[1]-d1);
      }

      return 2;
    case 2:
      p[0] = cellid_t (c[0]+1,c[1]+1);
      p[1] = cellid_t (c[0]+1,c[1]-1);
      p[2] = cellid_t (c[0]-1,c[1]-1);
      p[3] = cellid_t (c[0]-1,c[1]+1);
      return 4;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }
  }

  uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
  {
    switch (getCellDim (c))
    {
    case 0:
      return 0;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        f[0] = cellid_t (c[0]+d0,c[1]+d1);
        f[1] = cellid_t (c[0]-d0,c[1]-d1);
      }

      return 2;
    case 2:
      f[0] = cellid_t (c[0]  ,c[1]+1);
      f[1] = cellid_t (c[0]  ,c[1]-1);
      f[2] = cellid_t (c[0]-1,c[1]);
      f[3] = cellid_t (c[0]+1,c[1]);
      return 4;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }
  }

  uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = 0;

    switch (getCellDim (c))
    {
    case 0:
      cf[0] = cellid_t (c[0]  ,c[1]+1);
      cf[1] = cellid_t (c[0]  ,c[1]-1);
      cf[2] = cellid_t (c[0]-1,c[1]);
      cf[3] = cellid_t (c[0]+1,c[1]);
      cf_ct =  4;
      break;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        cf[0] = cellid_t (c[0]+d1,c[1]+d0);
        cf[1] = cellid_t (c[0]-d1,c[1]-d0);
        cf_ct =  2;
      }

      break;
    case 2:
      return 0;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }

    // position in cf[] where the next valid cf should be placed
    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;

  }

  uint dataset_t::getMaxCellDim() const
  {
    return 2;
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

  bool dataset_t::isCellMarked (cellid_t c) const
  {
    return ! ((*m_cell_flags) (c) == CELLFLAG_UNKNOWN);
  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_CRITCAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_PAIRED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    (*m_cell_pairs) (c) = p;
    (*m_cell_pairs) (p) = c;

    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_PAIRED;
    (*m_cell_flags) (p) = (*m_cell_flags) (p) |CELLFLAG_PAIRED;
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_CRITCAL;
  }

  bool dataset_t::isTrueBoundryCell (cellid_t c) const
  {
    return (m_ext_rect.isOnBoundry (c));
  }

  bool dataset_t::isFakeBoundryCell (cellid_t c) const
  {
    return (m_rect.isOnBoundry (c) && (!m_ext_rect.isOnBoundry (c)));
  }

  bool dataset_t::isCellExterior (cellid_t c) const
  {
    return (!m_rect.contains (c) && m_ext_rect.contains (c));
  }

  std::string dataset_t::getCellFunctionDescription (cellid_t c) const
  {
    std::stringstream ss;

    ( (std::ostream &) ss) <<c;

    return ss.str();

  }

  std::string dataset_t::getCellDescription (cellid_t c) const
  {

    std::stringstream ss;

    ( (std::ostream &) ss) <<c;

    return ss.str();

  }

  void dataset_t::work()
  {
    assignGradients();

    collateCriticalPoints();

    assignCellOwnerExtrema();
  }

  void  dataset_t::assignGradients()
  {

    // determine all the pairings of all cells in m_rect
    for (cell_coord_t y = m_rect.lower_corner()[1]; y <= m_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y),p;

      if (isCellMarked (c))
        continue;

      if (lowestPairableCoFacet (this,c,p))
        pairCells (c,p);
    }

    for (cell_coord_t y = m_rect.lower_corner()[1]; y <= m_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y);

      if (!isCellMarked (c)) markCellCritical (c);
    }

    // mark artificial boundry as critical

    for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
    {
      cellid_t bcs[] = {cellid_t (x,m_rect.lower_corner()[1]),cellid_t (x,m_rect.upper_corner()[1]) };

      for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
      {
        cellid_t &c = bcs[i];

        if (isCellCritical (c)) continue;

        cellid_t cf[20];

        u_int cf_ct =  getCellCofacets (c,cf);

        for (u_int j = 0 ; j <cf_ct;++j)
        {
          if (isCellExterior (cf[j]))
          {
            markCellCritical (c);
            markCellCritical (getCellPairId (c));
            break;
          }
        }
      }
    }

    for (cell_coord_t y = m_rect.lower_corner()[1] +1; y < m_rect.upper_corner()[1];y += 1)
    {
      cellid_t bcs[] = {cellid_t (m_rect.lower_corner()[0],y),cellid_t (m_rect.upper_corner()[0],y) };

      for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
      {
        cellid_t &c = bcs[i];

        if (isCellCritical (c)) continue;

        cellid_t cf[20];

        u_int cf_ct =  getCellCofacets (c,cf);

        for (u_int j = 0 ; j <cf_ct;++j)
        {
          if (isCellExterior (cf[j]))
          {
            markCellCritical (c);
            markCellCritical (getCellPairId (c));
            break;
          }
        }
      }
    }
  }

  void  dataset_t::collateCriticalPoints()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y);

      if (isCellCritical (c))
        m_critical_cells.push_back(c);
    }
  }


  void  dataset_t::assignCellOwnerExtrema()
  {
    for (cellid_list_t::iterator it = m_critical_cells.begin() ;
    it != m_critical_cells.end();++it)
    {

      (*m_cell_own)(*it) = *it;

      switch (getCellDim (*it))
      {
      case 0:
        track_gradient_tree_bfs(this,*it,GRADDIR_ASCENDING);
        break;
      case 2:
        track_gradient_tree_bfs(this,*it,GRADDIR_DESCENDING);
        break;
      default:
        break;
      }
    }

  }

  void  dataset_t::writeout_connectivity(mscomplex_t *msgraph)
  {

    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c));
    }

    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      if(!isCellPaired(c))  continue;

      uint cp_idx = i;

      msgraph->m_cps[cp_idx]->is_paired = true;

      msgraph->m_cps[cp_idx]->pair_idx =
          msgraph->m_id_cp_map[getCellPairId(c)];
    }

    for (cellid_list_t::iterator it = m_critical_cells.begin() ;
    it != m_critical_cells.end();++it)
    {
      cellid_t c = *it;

      if(getCellDim(c) == 1)
      {
        cellid_t f[4],cf[4];

        uint f_ct = getCellFacets(c,f);
        uint cf_ct = getCellCofacets(c,cf);

        for(uint i = 0 ; i < f_ct;++i)
        {
          cellid_t f_own_cp = (*m_cell_own)(f[i]);

          if(f_own_cp != cellid_t(-1,-1))
            connectCps(msgraph,c,f_own_cp);
        }

        for(uint i = 0 ; i < cf_ct;++i)
        {
          cellid_t cf_own_cp = (*m_cell_own)(cf[i]);

          if(cf_own_cp != cellid_t(-1,-1))
            connectCps(msgraph,c,cf_own_cp);
        }
      }
    }
  }

  void dataset_t::getCellCoord (cellid_t c,double &x,double &y,double &z)
  {
    x = c[0];
    y = 0;
    z = c[1];

    cellid_t pts[20];

    if(m_ext_rect.contains(c))
    {
      y= get_cell_fn(c);

    }
  }

  void dataset_t::log_flags()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
    {
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
        cellid_t c(x,y);

        int val = (*m_cell_flags)(c);

        std::cout<<val<<" ";
      }
      std::cout<<std::endl;
    }
  }

  void dataset_t::log_pairs()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
    {
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
        cellid_t c(x,y);
        std::cout<<(*m_cell_pairs)(c)<<" ";
      }
      std::cout<<std::endl;
    }
  }
}

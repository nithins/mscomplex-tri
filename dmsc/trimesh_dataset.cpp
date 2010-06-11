#include <trimesh_dataset.h>
#include <trimesh_dataset_ensure.h>

#include <vector>
#include <queue>

#include <timer.h>

#include <QFile>
#include <logutil.h>

#include <trimesh_mscomplex.h>

namespace trimesh
{

  class pt_comp_t
  {
    dataset_t *pOwn;
  public:
    pt_comp_t(dataset_t *o):pOwn(o){}

    bool operator()(cellid_t c1,cellid_t c2)
    {
      return pOwn->ptLt(c1,c2);
    }
  };

  inline bool dataset_t::ptLt(cellid_t c1, cellid_t c2) const
  {
    cell_fn_t f1 = m_vert_fns[c1];
    cell_fn_t f2 = m_vert_fns[c2];

    if (f1 != f2)
      return f1 < f2;

    return c1<c2;
  }


  cellid_t get_cp_cellid(mscomplex_t *msgraph,uint idx)
  {
    return msgraph->m_cps[idx]->cellid;
  }

  static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
  {
    &dataset_t::getCellFacets,
    &dataset_t::getCellCofacets
  };

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
//          if ( !dataset->isCellExterior ( cets[i] ) )
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
    msgraph->connect_cps(c1,c2);
  }


  inline bool lowestPairableCoFacet
      (dataset_t *dataset,cellid_t cellId,cellid_t& pairid
       )
  {
    cellid_t cofacets[20];
    bool    cofacet_usable[20];

    uint cofacet_count = dataset->getCellCofacets ( cellId,cofacets );

    bool isTrueBoundryCell = dataset->isBoundryCell ( cellId ) ;

    // for each co facet
    for ( uint i = 0 ; i < cofacet_count ; i++ )
    {
      cellid_t facets[20];
      uint facet_count = dataset->getCellFacets ( cofacets[i],facets );

      cofacet_usable[i] = true;

      if ( isTrueBoundryCell &&
           !dataset->isBoundryCell ( cofacets[i] ) )
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

      dataset->m_cell_own[top_cell] = start_cellId;

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
//          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            cellid_t next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              dataset->m_cell_own[cets[i]] = start_cellId;

              // mark here that the parent of next cell is top_cell
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }


  dataset_t::dataset_t () :
      m_ptcomp(new pt_comp_t(this))
  {
  }

  dataset_t::~dataset_t ()
  {
    delete m_ptcomp;

    clear();
  }

  void dataset_t::init(const cell_fn_list_t &vert_fns,const tri_idx_list_t & trilist)
  {
    ensure_valid_trilist_indexes(trilist,vert_fns.size());

    m_vert_fns.resize(vert_fns.size());

    std::copy(vert_fns.begin(),vert_fns.end(),m_vert_fns.begin());

    m_tri_cc.init(trilist,vert_fns.size());

    m_cell_flags.resize(m_tri_cc.get_num_cells(),0);

    m_cell_pairs.resize(m_tri_cc.get_num_cells(),-1);

    m_cell_own.resize(m_tri_cc.get_num_cells(),-1);
  }

  void  dataset_t::clear()
  {
    m_vert_fns.clear();

    m_cell_flags.clear();

    m_cell_own.clear();

    m_cell_pairs.clear();

    m_tri_cc.clear();

    m_critical_cells.clear();
  }

  cellid_t   dataset_t::getCellPairId (cellid_t c) const
  {
    ensure_cell_paired(this,c);

    return m_cell_pairs[c];
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    if(getCellDim(c1) == 0)
      return ptLt(c1,c2);

    cellid_t pts1[20];
    cellid_t pts2[20];

    uint pts1_ct = getCellPoints ( c1,pts1);
    uint pts2_ct = getCellPoints ( c2,pts2);

    std::sort ( pts1,pts1+pts1_ct,*m_ptcomp );
    std::sort ( pts2,pts2+pts2_ct,*m_ptcomp);

    return std::lexicographical_compare
        ( pts1,pts1+pts1_ct,pts2,pts2+pts2_ct,*m_ptcomp );
  }

  cell_fn_t dataset_t::get_cell_fn (cellid_t c) const
  {

    cell_fn_t  fn = 0.0;

    cellid_t pts[20];

    uint pts_ct = getCellPoints (c,pts);

    for (int j = 0 ; j <pts_ct ;++j)
      fn += m_vert_fns[pts[j]];

    fn /= pts_ct;

    return fn;
  }

  uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
  {
    return m_tri_cc.get_cell_points(c,p);
  }

  uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
  {
    return m_tri_cc.get_cell_facets(c,f);
  }

  uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
  {
    return m_tri_cc.get_cell_co_facets(c,cf);
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

  bool dataset_t::isCellMarked (cellid_t c) const
  {
    return ! (m_cell_flags[c] == CELLFLAG_UNKNOWN);
  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return (m_cell_flags[c]& CELLFLAG_CRITCAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return (m_cell_flags[c] & CELLFLAG_PAIRED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    m_cell_pairs[c] = p;
    m_cell_pairs[p] = c;

    m_cell_flags[c] = m_cell_flags[c] |CELLFLAG_PAIRED;
    m_cell_flags[p] = m_cell_flags[p] |CELLFLAG_PAIRED;
  }

  uint dataset_t::getCellDim ( cellid_t c ) const
  {
    return m_tri_cc.get_cell_dim(c);
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    m_cell_flags[c] = m_cell_flags[c] |CELLFLAG_CRITCAL;
  }

  bool dataset_t::isBoundryCell (cellid_t c) const
  {
    return m_tri_cc.is_cell_boundry(c);
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

    uint num_cells      = m_tri_cc.get_num_cells();

    uint num_cells_ltmd = m_tri_cc.get_num_cells_max_dim(m_tri_cc.get_dim() -1);

    for(uint i = 0 ; i < num_cells_ltmd; ++i)
    {
      cellid_t c = i,p;

      if (isCellMarked (c))
        continue;

      if (lowestPairableCoFacet (this,c,p))
        pairCells (c,p);
    }

    for(uint i = 0 ; i < num_cells;++i)
    {
      cellid_t c = i;

      if (!isCellMarked (c))
        markCellCritical (c);
    }
  }

  void  dataset_t::collateCriticalPoints()
  {
    uint num_cells = m_tri_cc.get_num_cells();

    for(uint i = 0 ; i < num_cells;++i )
    {
      cellid_t c = i;

      if (isCellCritical (c))
        m_critical_cells.push_back(c);
    }
  }


  void  dataset_t::assignCellOwnerExtrema()
  {
    for (cellid_list_t::iterator it = m_critical_cells.begin() ;
    it != m_critical_cells.end();++it)
    {

      m_cell_own[*it] = *it;

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

      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c),isBoundryCell(c));
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
          cellid_t f_own_cp = m_cell_own[f[i]];

          if(f_own_cp != invalid_cellid)
            connectCps(msgraph,c,f_own_cp);
        }

        for(uint i = 0 ; i < cf_ct;++i)
        {
          cellid_t cf_own_cp = m_cell_own[cf[i]];

          if(cf_own_cp != invalid_cellid)
            connectCps(msgraph,c,cf_own_cp);
        }
      }
    }
  }

  void dataset_t::log_flags()
  {
    for(uint i = 0 ; i < m_cell_flags.size(); ++i)
      std::cout<<(int)m_cell_flags[i]<<" ";

    std::cout<<std::endl;

  }

  void dataset_t::log_pairs()
  {
    for(uint i = 0 ; i < m_cell_flags.size(); ++i)
      std::cout<<m_cell_pairs[i]<<" ";

    std::cout<<std::endl;

  }
}

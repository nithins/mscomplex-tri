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


#ifndef __GRID_DATASET_H_INCLUDED_
#define __GRID_DATASET_H_INCLUDED_

#include <vector>
#include <trimesh.h>

namespace trimesh
{

  class mscomplex_t;

  class pt_comp_t;

  class dataset_t
  {

  public:

    enum eCellFlags
    {
      CELLFLAG_UNKNOWN = 0,
      CELLFLAG_PAIRED  = 1,
      CELLFLAG_CRITCAL = 2,
    };

  public:

    cell_fn_list_t      m_vert_fns;
    cellid_list_t       m_cell_own;
    cellid_list_t       m_cell_pairs;
    cellid_list_t       m_critical_cells;
    std::vector<uchar>  m_cell_flags;

    pt_comp_t*          m_ptcomp;
    tri_cell_complex_t  m_tri_cc;

  public:

    // initialization of the dataset

    dataset_t ();

    ~dataset_t ();

    void  init(const cell_fn_list_t &vert_fns,const tri_idx_list_t & trilist);

    void  clear();


    // actual algorithm work
  public:

    void  work();

    void  writeout_connectivity(mscomplex_t *msgraph);

    void  assignGradients();

    void  collateCriticalPoints();

    void  assignCellOwnerExtrema();

    int   postMergeFillDiscs(mscomplex_t *msgraph);

    // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    inline bool ptLt ( cellid_t c1,cellid_t c2) const;

    bool   compareCells( cellid_t ,cellid_t ) const;

    uint   getCellPoints ( cellid_t ,cellid_t  * ) const;

    uint   getCellFacets ( cellid_t ,cellid_t * ) const;

    uint   getCellCofacets ( cellid_t ,cellid_t * ) const;

    bool   isPairOrientationCorrect ( cellid_t c, cellid_t p ) const;

    bool   isCellMarked ( cellid_t c ) const;

    bool   isCellCritical ( cellid_t c ) const;

    bool   isCellPaired ( cellid_t c ) const;

    void   pairCells ( cellid_t c,cellid_t p );

    void   markCellCritical ( cellid_t c );

    uint   getCellDim ( cellid_t c ) const;

    bool   isBoundryCell ( cellid_t c ) const;

    std::string  getCellFunctionDescription ( cellid_t pt ) const;

    std::string getCellDescription ( cellid_t cellid ) const;

    // misc functions
  public:

    void log_flags();

    void log_pairs();

    // return fn at point .. averge of points for higher dims
    cell_fn_t get_cell_fn ( cellid_t c ) const;

  };
}
#endif

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

#ifndef GRID_DATAMANAGER_H_INCLUDED_
#define GRID_DATAMANAGER_H_INCLUDED_

#include <fstream>
#include <vector>

#include <grid_dataset.h>

namespace boost
{
  class thread;
}

namespace grid
{
  class  mscomplex_t;

  struct datapiece_t
  {
    dataset_t   *dataset;
    mscomplex_t *msgraph;

    uint level;

    uint m_pieceno;

    datapiece_t (uint l);

    std::string label();
  };



  class data_manager_t
  {
    typedef std::vector<datapiece_t *> pieces_list_t;

  public:

    pieces_list_t                m_pieces;

    std::string                  m_filename;
    cellid_t                     m_size;

    u_int                        m_num_levels;
    double                       m_simp_tresh;
    bool                         m_single_threaded_mode;
    bool                         m_compute_out_of_core;

    boost::thread **             m_threads;

  public:

    uint num_parallel;


    data_manager_t
        ( std::string  filename,
          cellid_t     size,
          u_int        num_levels,
          bool         threaded_mode,
          double       simp_tresh,
          bool         compute_out_of_core,
          uint         np);

    virtual ~data_manager_t ();

    void createPieces_quadtree(rect_t r,rect_t e,u_int level );

    void createDataPieces();

    void readDataAndInit(std::ifstream &data_stream,cell_fn_t *,uint start_offset);

    uint getMaxDataBufItems();

    void waitForThreadsInRange(uint,uint);

    void computeMsGraph ( datapiece_t  * );

    void computeMsGraphInRange(uint ,uint );

    void finalMergeDownPiecesInRange(uint start,uint end);

    void collectManifold( datapiece_t  * );

    void collectManifoldsInRange(uint start,uint end);

    void writeManifoldsInRange(uint start,uint end);

    void computeSubdomainMsgraphs ();

    void mergePiecesUp( );

    void mergePiecesDown( );

    void collectSubdomainManifolds( );

    void logAllConnections(const std::string &prefix);

    void logAllCancelPairs(const std::string &prefix);

  };
}

#endif

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

namespace trimesh
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

    pieces_list_t  m_pieces;

    std::string    m_tri_filename;
    std::string    m_bin_filename;
    uint           m_bin_comp_no;
    double         m_simp_tresh;

  public:

    data_manager_t (
        std::string tri_file,
        std::string bin_file,
        uint bin_comp,
        double simp_tresh
        );

    ~data_manager_t ();

    void work();

    void createDataPieces();

    void destroyDataPieces();

    void init();

    void logAllConnections(const std::string &prefix);

    void logAllCancelPairs(const std::string &prefix);

  };
}

#endif

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

#ifndef TRIMESH_DATAMANAGER_H_INCLUDED
#define TRIMESH_DATAMANAGER_H_INCLUDED

#include <fstream>
#include <vector>

#include <trimesh.h>

namespace boost
{
  class thread;
}

namespace trimesh
{
//  struct datapiece_t
//  {
//    dataset_t   *dataset;
//    mscomplex_t *msgraph;

//    uint level;

//    uint m_pieceno;

//    datapiece_t (uint l);

//    std::string label();
//  };

  class data_manager_t
  {

  public:

    dataset_ptr_t   m_dataset;
    mscomplex_ptr_t m_msgraph;


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

    void init();

    void work();

    void save();

  };
}

#endif

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

#include <iostream>
#include <fstream>

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/regex.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <timer.h>
#include <logutil.h>

#include <grid_datamanager.h>

#include <grid_mscomplex.h>


namespace trimesh
{

  using namespace std;

  void data_manager_t::init()
  {
    tri_idx_list_t tlist;
    vertex_list_t  vlist;

    glutils::read_tri_file(m_tri_filename.c_str(),vlist,tlist);

    cell_fn_list_t cell_fns;

    cell_fns.resize(vlist.size());


    // read in function values

    const   uint fnname_max_size = 32;

    typedef float bin_data_type_t;

    fstream fnfile ( m_bin_filename.c_str(), fstream::in | fstream::binary );

    int num_bin_values, num_bin_comps;

    fnfile.read ( reinterpret_cast<char *> ( &num_bin_values ), sizeof ( int ) );
    fnfile.read ( reinterpret_cast<char *> ( &num_bin_comps ), sizeof ( int ) );

    if(num_bin_values != vlist.size())
      throw std::runtime_error("tri / bin file num verts mismatch");

    if(!(m_bin_comp_no < num_bin_comps))
      throw std::runtime_error("selected scalar component does not exist");

    char *fnnames = new char[num_bin_comps * fnname_max_size];

    for ( uint i = 0; i < ( uint ) num_bin_comps; i++ )
      fnfile.read ( fnnames + i * fnname_max_size, fnname_max_size );

    fnfile.seekg ( sizeof ( bin_data_type_t ) * m_bin_comp_no, ios::cur );

    for ( uint i = 0; i < ( uint ) num_bin_values; i++ )
    {
      bin_data_type_t data;

      fnfile.read ( reinterpret_cast<char *> ( &data ),sizeof(bin_data_type_t));
      fnfile.seekg ( sizeof ( bin_data_type_t)*( num_bin_comps-1), ios::cur );

      cell_fns[i] = data;
    }

    fnfile.close();

    _LOG("------------------------");
    _LOG(" scalar component names ");
    _LOG("------------------------");
    for(uint i = 0 ; i < num_bin_comps;++i)
    {
      _LOG(fnnames + i * fnname_max_size);
    }
    _LOG("------------------------");
    _LOG("selected comp = "<< fnnames + m_bin_comp_no * fnname_max_size);
    _LOG("------------------------");


    delete fnnames;

    datapiece_t * dp = m_pieces[0];

    dp->dataset->init(cell_fns,tlist);
  }

  void data_manager_t::createDataPieces()
  {
    datapiece_t * dp = new datapiece_t(0);

    dp->dataset = new dataset_t();
    dp->msgraph = new mscomplex_t();

    m_pieces.push_back(dp);
  }

  void data_manager_t::destroyDataPieces()
  {
    for(uint i = 0 ; i < m_pieces.size(); ++i)
    {
      datapiece_t * dp = m_pieces[i];

      delete dp->msgraph;
      delete dp->dataset;
      delete dp;
    }
    m_pieces.clear();
  }

  void data_manager_t::work()
  {
    Timer t;
    t.start();

    datapiece_t * dp = m_pieces[0];

    _LOG ( "===================" );
    _LOG ( "Starting Processing" );
    _LOG ( "-------------------" );

    init();

    _LOG ("timer_time = "<<t.getElapsedTimeInMilliSec());

    dp->dataset->work();

    _LOG ("timer_time = "<<t.getElapsedTimeInMilliSec());

    dp->dataset->writeout_connectivity(dp->msgraph);

    _LOG ("timer_time = "<<t.getElapsedTimeInMilliSec());

    if(m_simp_tresh > 0.0)
      dp->msgraph->simplify_un_simplify(m_simp_tresh);

    _LOG ("timer_time = "<<t.getElapsedTimeInMilliSec());

    dp->dataset->postMergeFillDiscs(dp->msgraph);

    _LOG ("timer_time = "<<t.getElapsedTimeInMilliSec());
    _LOG ( "-------------------" );
    _LOG ( "Finished Processing" );
    _LOG ( "===================" );
  }


  data_manager_t::data_manager_t(
      std::string tri_file,
      std::string bin_file,
      uint bin_comp,
      double simp_tresh
      ):
      m_tri_filename(tri_file),
      m_bin_filename(bin_file),
      m_bin_comp_no(bin_comp),
      m_simp_tresh(simp_tresh)
  {
    createDataPieces();
  }

  data_manager_t::~data_manager_t()
  {
    destroyDataPieces();
  }

  void data_manager_t ::logAllConnections(const std::string &prefix)
  {

    for(uint i = 0 ; i <m_pieces.size();++i)
    {
      datapiece_t *dp = m_pieces[i];

      std::string filename(prefix+dp->label()+string(".txt"));

      ofstream outfile;
      outfile.open(filename.c_str(),  std::ios::out|std::ios::trunc);

      if(outfile.is_open() == false )
      {
        _LOG("failed to open log file");
        break;
      }

      std::stringstream ss;

      dp->msgraph->print_connections( (ostream&)ss);

      outfile<<ss.str();
    }

  }

  void data_manager_t::logAllCancelPairs(const std::string &prefix)
  {
  }

  datapiece_t::datapiece_t (uint pno):
      dataset(NULL),
      msgraph(NULL),
      level(0),
      m_pieceno(pno)
  {
  }

  std::string datapiece_t::label()
  {
    std::stringstream ss;
    ss<<m_pieceno;

    return ss.str();
  }
}

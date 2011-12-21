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

#include <trimesh_datamanager.h>

#include <trimesh_mscomplex.h>
#include <trimesh_dataset.h>


namespace trimesh
{

  using namespace std;

  typedef float bin_data_type_t;

  const   uint bin_fnname_max_size = 32;

  template <typename T>
  void read_bin_file(std::vector<T> &cell_fns, const string & fname,int compno)
  {
    fstream fnfile ( fname.c_str(), fstream::in | fstream::binary );

    int num_bin_values, num_bin_comps;

    fnfile.read ( reinterpret_cast<char *> ( &num_bin_values ), sizeof ( int ) );
    fnfile.read ( reinterpret_cast<char *> ( &num_bin_comps ), sizeof ( int ) );

    cell_fns.resize(num_bin_values,-1);

    fnfile.seekg ( bin_fnname_max_size*num_bin_comps, ios::cur );
    fnfile.seekg ( sizeof ( bin_data_type_t ) * compno, ios::cur );

    for ( uint i = 0; i < ( uint ) num_bin_values; i++ )
    {
      bin_data_type_t data;

      fnfile.read ( reinterpret_cast<char *> ( &data ),sizeof(bin_data_type_t));
      fnfile.seekg ( sizeof ( bin_data_type_t)*( num_bin_comps-1), ios::cur );

      cell_fns[i] = data;
    }

    fnfile.close();
  }

  void print_bin_info(const string & fname)
  {
    fstream fnfile ( fname.c_str(), fstream::in | fstream::binary );

    int num_bin_values, num_bin_comps;

    fnfile.read ( reinterpret_cast<char *> ( &num_bin_values ), sizeof ( int ) );
    fnfile.read ( reinterpret_cast<char *> ( &num_bin_comps ), sizeof ( int ) );

    char compname[bin_fnname_max_size];

    cout<<"        component names             "<<endl;
    cout<<"------------------------------------"<<endl;

    for(uint i = 0 ; i < num_bin_comps;++i)
    {
      fnfile.read ( compname, bin_fnname_max_size );
      cout<<i<<". "<<compname<<endl;
    }
    cout<<"------------------------------------"<<endl;
  }

  void read_tri_file( const char *filename,tri_idx_list_t &tlist)
  {
    uint num_v,num_t;

    std::fstream tri_file ( filename, std::fstream::in );

    if(tri_file.is_open() == false)
      throw std::runtime_error("unable to open tri file");

    tri_file >> num_v >> num_t;

    tlist.resize(num_t);

    double vx,vy,vz;

    for ( uint i = 0; i < num_v; ++i )
      tri_file>>vx>>vy>>vz;

    for ( uint i = 0; i < num_t; i++ )
      for ( uint j = 0; j < 3; ++j )
      {
        tri_file >> tlist[i][j];

        if(tlist[i][j] >= num_v||tlist[i][j] < 0)
        {
          throw std::runtime_error("invalid index");
        }
      }

    tri_file.close();
  }


  void data_manager_t::init()
  {
    tri_idx_list_t tlist;
    cell_fn_list_t cell_fns;

    read_tri_file(m_tri_filename.c_str(),tlist);

    print_bin_info(m_bin_filename);

    read_bin_file(cell_fns,m_bin_filename,m_bin_comp_no);

    cout<<"selected comp = "<<m_bin_comp_no<<endl;
    cout<<"------------------------------------"<<endl;

    m_dataset->init(cell_fns,tlist);
  }

  void data_manager_t::work()
  {
    Timer t;
    t.start();

    cout<<"===================================="<<endl;
    cout<<"         Starting Processing        "<<endl;
    cout<<"------------------------------------"<<endl;

    init();
    cout<<"data read ---------------- "<<t.getElapsedTimeInMilliSec()<<endl;

    m_dataset->work(m_msgraph);
    cout<<"gradient done ------------ "<<t.getElapsedTimeInMilliSec()<<endl;

    m_msgraph->stow(m_tri_filename+".full.graph.bin",false);
    cout<<"write graph done --------- "<<t.getElapsedTimeInMilliSec()<<endl;

    m_msgraph->simplify(m_simp_tresh);
    m_msgraph->un_simplify();
    cout<<"simplification done ------ "<<t.getElapsedTimeInMilliSec()<<endl;

    m_msgraph->stow(m_tri_filename+".graph.bin",false);
    cout<<"write graph done --------- "<<t.getElapsedTimeInMilliSec()<<endl;

    m_msgraph->invert_for_collection();
    m_dataset->save_manifolds(m_tri_filename+".mfold.bin",m_msgraph);
    cout<<"write mfolds done --------- "<<t.getElapsedTimeInMilliSec()<<endl;

    cout<<"------------------------------------"<<endl;
    cout<<"        Finished Processing         "<<endl;
    cout<<"===================================="<<endl;

  }

  void data_manager_t::save()
  {
/*    _LOG ( "===================" );
    _LOG ( " Saving Results    " );
    _LOG ( "-------------------" );

    _LOG ( " msgraph.txt       " );
    _LOG ( "-------------------" );
    {
      std::fstream f("msgraph.txt",fstream::out);
      m_msgraph->save(f);
      f.close();
    }

    _LOG ( " msmfolds.txt      " );
    _LOG ( "-------------------" );
    {
      vertex_list_t  vlist;
      tri_cc_geom_t  geom;

      glutils::read_tri_file(m_tri_filename.c_str(),vlist);
      geom.init(m_dataset->m_tri_cc,vlist);

      std::fstream f("msmfolds.txt",fstream::out);
      m_msgraph->save_manifolds(f,geom);
      f.close();
    }

    _LOG ( " Finished Saving   " );
    _LOG ( "===================" );*/
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
      m_simp_tresh(simp_tresh),
      m_dataset(new dataset_t),
      m_msgraph(new mscomplex_t)
  {
  }

  data_manager_t::~data_manager_t()
  {
  }

//  datapiece_t::datapiece_t (uint pno):
//      dataset(NULL),
//      msgraph(NULL),
//      level(0),
//      m_pieceno(pno)
//  {
//  }

//  std::string datapiece_t::label()
//  {
//    std::stringstream ss;
//    ss<<m_pieceno;

//    return ss.str();
//  }
}

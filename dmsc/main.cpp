#include <stdexcept>
#include <iostream>
#include <exception>
#include <string>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include <timer.h>

#include <trimesh_dataset.h>
#include <trimesh_mscomplex.h>
#include <trimesh_mscomplex_simp.h>

using namespace std;
namespace bpo = boost::program_options;

typedef float bin_data_type_t;

const   uint bin_fnname_max_size = 32;

template <typename T>
void read_bin_file(std::vector<T> &fns, const string & fname,int compno)
{
  fstream fnfile ( fname.c_str(), fstream::in | fstream::binary );

  ensure(fnfile.is_open(),"unable to open bin file");

  int num_bin_values, num_bin_comps;

  fnfile.read ( reinterpret_cast<char *> ( &num_bin_values ), sizeof ( int ) );
  fnfile.read ( reinterpret_cast<char *> ( &num_bin_comps ), sizeof ( int ) );

  fns.resize(num_bin_values,-1);

  fnfile.seekg ( bin_fnname_max_size*num_bin_comps, ios::cur );
  fnfile.seekg ( sizeof ( bin_data_type_t ) * compno, ios::cur );

  for ( uint i = 0; i < ( uint ) num_bin_values; i++ )
  {
    bin_data_type_t data;

    fnfile.read ( reinterpret_cast<char *> ( &data ),sizeof(bin_data_type_t));
    fnfile.seekg ( sizeof ( bin_data_type_t)*( num_bin_comps-1), ios::cur );

    fns[i] = data;
  }

  fnfile.close();
}

void print_bin_info(const string & fname)
{
  fstream fnfile ( fname.c_str(), fstream::in | fstream::binary );

  ensure(fnfile.is_open(),"unable to open bin file");

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

void read_tri_tlist( const char *filename,trimesh::tri_idx_list_t &tlist)
{
  uint num_v,num_t;

  std::fstream tri_file ( filename, std::fstream::in );

  ensure(tri_file.is_open(),"unable to open tri file");

  tri_file >> num_v >> num_t;

  tlist.resize(num_t);

  double vx,vy,vz;

  for ( uint i = 0; i < num_v; ++i )
    tri_file>>vx>>vy>>vz;

  for ( uint i = 0; i < num_t; i++ )
    for ( uint j = 0; j < 3; ++j )
    {
      tri_file >> tlist[i][j];

      ensure(is_in_range(tlist[i][j],0,num_v),"invalid index in file");
    }

  tri_file.close();
}

void read_tri_vlist( const char *filename,tri_cc_geom_t::vertex_list_t &vlist)
{
  uint num_v,num_t;

  std::fstream tri_file ( filename, std::fstream::in );

  ensure(tri_file.is_open(),"unable to open tri file");

  tri_file >> num_v >> num_t;

  vlist.resize(num_v);

  for ( uint i = 0; i < num_v; ++i )
    tri_file>>vlist[i][0]
            >>vlist[i][1]
            >>vlist[i][2];


  tri_file.close();
}

namespace ba = boost::algorithm;

template <typename T>
void read_off_file(const string & fname, std::vector<T> &fns,
                   trimesh::tri_idx_list_t &tlist ,int compno)
{
  fstream off_file ( fname.c_str(), fstream::in);

  ensure(off_file.is_open(),"unable to open off file");

  string ln;
  std::vector<string> strs;

  getline(off_file,ln);
  ensure(ln=="OFF","Doesn't seem to be an OFF FILE");

  getline(off_file,ln);
  ba::split(strs,ln,ba::is_any_of("\t \n"));

  int num_v = atoi(strs[0].c_str());
  int num_t = atoi(strs[1].c_str());

  fns.resize(num_v);
  tlist.resize(num_t);

  for ( uint i = 0; i < num_v; ++i )
  {
    getline(off_file,ln);
    ba::split(strs,ln,ba::is_any_of("\t \n"));
    fns[i] = atof(strs[compno].c_str());

    ensure(is_in_range(compno,0,strs.size()),
           "too few components in vinfo line");
  }

  for ( uint i = 0; i < num_t; i++ )
  {
    getline(off_file,ln);
    ba::split(strs,ln,ba::is_any_of("\t \n"));

    int ntv     = atoi(strs[0].c_str());
    tlist[i][0] = atoi(strs[1].c_str());
    tlist[i][1] = atoi(strs[2].c_str());
    tlist[i][2] = atoi(strs[3].c_str());

    ensure(ntv == 3,"Mesh contains non-triangle polys");
    ensure(is_in_range(tlist[i][0],0,num_v),"invalid index in file");
    ensure(is_in_range(tlist[i][1],0,num_v),"invalid index in file");
    ensure(is_in_range(tlist[i][2],0,num_v),"invalid index in file");
  }

  off_file.close();
}

void read_off_vlist(const string & fname, tri_cc_geom_t::vertex_list_t &vlist)
{
  fstream off_file ( fname.c_str(), fstream::in);

  ensure(off_file.is_open(),"unable to open off file");

  string ln;
  std::vector<string> strs;

  getline(off_file,ln);
  ensure(ln=="OFF","Doesn't seem to be an OFF FILE");

  getline(off_file,ln);
  ba::split(strs,ln,ba::is_any_of("\t \n"));

  int num_v = atoi(strs[0].c_str());
  vlist.resize(num_v);

  for ( uint i = 0; i < num_v; ++i )
  {
    getline(off_file,ln);
    ba::split(strs,ln,ba::is_any_of("\t \n"));

    vlist[i][0] = atof(strs[0].c_str());
    vlist[i][1] = atof(strs[1].c_str());
    vlist[i][2] = atof(strs[2].c_str());
  }

  off_file.close();
}


int main(int ac , char **av)
{
  string tri_filename;
  string bin_filename;
  string off_filename;
  string simp_method;

  int    comp_no = 0;
  double simp_tresh  = 0.0;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("tri-file,t",bpo::value(&tri_filename)->default_value(""),
       "tri file name")
      ("bin-file,b",bpo::value(&bin_filename)->default_value(""),
       "bin file name (function file)")
      ("off-file,o",bpo::value(&off_filename)->default_value(""),
       "off file name")
      ("comp-no,c",bpo::value(&comp_no)->default_value(0),
       "scalar component number to use for the MS compelex")
      ("simp-tresh,s",bpo::value(&simp_tresh)->default_value(0.0),
       "simplification treshold\n"\
       "\n"\
       "====For simp-method = \"P\"=========\n"\
       "if s in [0,1] then simplify all features having pers < s\n"\
       "if s < 0 then simplify till there are int(-s) minima\n"\
       "if s > 1 then simplify till there are int(s)  maxima\n"\
       "\n")
      ("simp-method",bpo::value(&simp_method)->default_value("P"),
       "simplification method to use\n"\
       "P  ----> Persistence\n"\
       "AWP ---> Area weighted persistence\n"\
       "ABP ---> Area before persistence");

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);

  if (vm.count("help"))
  {
    cout << desc << endl;
    return 0;
  }
  try
  {
    bpo::notify(vm);
  }
  catch(bpo::required_option e)
  {
    cout<<e.what()<<endl;
    cout<<desc<<endl;
    return 1;
  }

  if((tri_filename.empty() || bin_filename.empty()) && off_filename.empty())
  {
    cout<<"Must specify either tri-bin or off file"<<endl;
    cout<<desc<<endl;
    return 1;
  }

  Timer t;
  t.start();

  cout<<"===================================="<<endl;
  cout<<"         Starting Processing        "<<endl;
  cout<<"------------------------------------"<<endl;

  trimesh::tri_idx_list_t tlist;
  trimesh::fn_list_t      fns;
  cout<<"selected comp = "<<comp_no<<endl;
  cout<<"------------------------------------"<<endl;

  string fn_pfx;

  if(off_filename.empty())
  {
    print_bin_info(bin_filename);
    read_tri_tlist(tri_filename.c_str(),tlist);
    read_bin_file(fns,bin_filename,comp_no);
    fn_pfx = tri_filename;
  }
  else
  {
    read_off_file(off_filename,fns,tlist,comp_no);
    fn_pfx = off_filename;
  }
  cout<<"data read ---------------- "<<t.getElapsedTimeInMilliSec()<<endl;

  trimesh::dataset_ptr_t   ds(new trimesh::dataset_t(fns,tlist));
  trimesh::mscomplex_ptr_t msc(new trimesh::mscomplex_t);
  ds->work(msc);
  cout<<"gradient done ------------ "<<t.getElapsedTimeInMilliSec()<<endl;

  msc->simplify(0.0);
  msc->un_simplify();
  msc->get_mfolds(ds);
  msc->save(fn_pfx+".mscomplex.full.bin");
  msc->save_ascii(fn_pfx+".mscomplex.full.txt");
  cout<<"write unsimplified done -- "<<t.getElapsedTimeInMilliSec()<<endl;


  if( simp_method == "P")
  {
    msc->simplify(simp_tresh);
    msc->un_simplify();
  }
  else if (simp_method == "AWP" || simp_method == "ABP")
  {
    tri_cc_geom_t::vertex_list_t vlist;

    if(off_filename.empty())
      read_tri_vlist(tri_filename.c_str(),vlist);
    else
      read_off_vlist(off_filename.c_str(),vlist);

    tri_cc_geom_ptr_t tcc(new tri_cc_geom_t);
    tcc->init(ds->m_tcc,vlist);

    if(simp_method == "AWP")
      trimesh::simplify_awp(msc,tcc,simp_tresh);
    else if(simp_method == "ABP")
      trimesh::simplify_abp(msc,tcc,simp_tresh);

    msc->un_simplify();
  }

  cout<<"simplification done ------ "<<t.getElapsedTimeInMilliSec()<<endl;

  msc->clear_mfolds();
  msc->get_mfolds(ds);
  msc->save(fn_pfx+".mscomplex.bin");
  msc->save_ascii(fn_pfx+".mscomplex.txt");
  cout<<"write simplified done ---- "<<t.getElapsedTimeInMilliSec()<<endl;

  cout<<"------------------------------------"<<endl;
  cout<<"        Finished Processing         "<<endl;
  cout<<"===================================="<<endl;
}

#include <stdexcept>
#include <iostream>
#include <exception>
#include <string>

#include <boost/program_options.hpp>

#include <timer.h>

#include <trimesh_dataset.h>
#include <trimesh_mscomplex.h>

using namespace std;
namespace bpo = boost::program_options;

typedef float bin_data_type_t;

const   uint bin_fnname_max_size = 32;

template <typename T>
void read_bin_file(std::vector<T> &fns, const string & fname,int compno)
{
  fstream fnfile ( fname.c_str(), fstream::in | fstream::binary );

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

void read_tri_file( const char *filename,trimesh::tri_idx_list_t &tlist)
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

int main(int ac , char **av)
{
  string tri_filename;
  string bin_filename;

  int    bin_comp_no = 0;
  double simp_tresh= 0.0;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("tri-file,t",bpo::value(&tri_filename)->required(), "tri file name")
      ("bin-file,b",bpo::value(&bin_filename)->required(), "bin file name (function file)")
      ("simp-tresh,s",bpo::value(&simp_tresh)->default_value(0.0),"simplification treshold")
      ("bin-file-comp,c",bpo::value(&bin_comp_no)->default_value(0),"scalar component number");


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

  Timer t;
  t.start();

  cout<<"===================================="<<endl;
  cout<<"         Starting Processing        "<<endl;
  cout<<"------------------------------------"<<endl;

  trimesh::tri_idx_list_t tlist;
  trimesh::fn_list_t      fns;
  print_bin_info(bin_filename);
  cout<<"selected comp = "<<bin_comp_no<<endl;
  cout<<"------------------------------------"<<endl;

  read_tri_file(tri_filename.c_str(),tlist);
  read_bin_file(fns,bin_filename,bin_comp_no);
  cout<<"data read ---------------- "<<t.getElapsedTimeInMilliSec()<<endl;

  trimesh::dataset_ptr_t   ds(new trimesh::dataset_t(fns,tlist));
  trimesh::mscomplex_ptr_t msc(new trimesh::mscomplex_t);
  ds->work(msc);
  cout<<"gradient done ------------ "<<t.getElapsedTimeInMilliSec()<<endl;

  msc->simplify(0.0);
  msc->un_simplify();
  msc->get_mfolds(ds);
  msc->save(tri_filename+".mscomplex.full.bin");
  msc->save_ascii(tri_filename+".mscomplex.full.txt");
  cout<<"write unsimplified done -- "<<t.getElapsedTimeInMilliSec()<<endl;

  msc->simplify(simp_tresh);
  msc->un_simplify();
  cout<<"simplification done ------ "<<t.getElapsedTimeInMilliSec()<<endl;

  msc->clear_mfolds();
  msc->get_mfolds(ds);
  msc->save(tri_filename+".mscomplex.bin");
  msc->save_ascii(tri_filename+".mscomplex.txt");
  cout<<"write simplified done ---- "<<t.getElapsedTimeInMilliSec()<<endl;

  cout<<"------------------------------------"<<endl;
  cout<<"        Finished Processing         "<<endl;
  cout<<"===================================="<<endl;
}

#include <stdexcept>
#include <iostream>
#include <exception>
#include <string>

#include <trimesh_datamanager.h>

#include <boost/program_options.hpp>


using namespace std;
namespace bpo = boost::program_options;

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

  trimesh::data_manager_t gdm(tri_filename,bin_filename,bin_comp_no,simp_tresh);

  gdm.work();
}

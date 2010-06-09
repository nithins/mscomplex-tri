#include <exception>
#include <string>

#include <trimesh_viewer_mainwindow.h>
#include <trimesh_datamanager.h>

#include <cpputils.h>

#include <boost/program_options.hpp>

#include <stdexcept>
#include <iostream>

using namespace std;

namespace bpo = boost::program_options ;

int main(int ac , char **av)
{
  string tri_filename;

  string bin_filename;

  int    bin_comp_no = 0;

  double simp_tresh= 0.0;

  bool   gui = false;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("tri-file,t",bpo::value<std::string >(), "tri file name")
      ("bin-file,b",bpo::value<std::string >(), "bin file name (function file)")
      ("simp-tresh,s",bpo::value<double>(),"simplification treshold")
      ("bin-file-comp,c",bpo::value<int>(),"scalar component number")
      ("gui,g","show gui")
      ;


  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);
  bpo::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 1;
  }

  if (vm.count("tri-file"))
    tri_filename = vm["tri-file"].as<std::string>();
  else
    throw std::invalid_argument("no tri filename specified");

  if (vm.count("bin-file"))
    bin_filename = vm["bin-file"].as<std::string>();
  else
    throw std::invalid_argument("no bin filename specified");

  if (vm.count("simp-tresh"))
    simp_tresh = vm["simp-tresh"].as<double>();

  if (vm.count("bin-file-comp"))
    bin_comp_no = vm["bin-file-comp"].as<int>();

  if (vm.count("gui"))
    gui = true;

  trimesh::data_manager_t * gdm = new trimesh::data_manager_t
      (tri_filename,bin_filename,bin_comp_no,simp_tresh);

  gdm->work();

  if(gui)
  {
      QApplication application(ac,av);

      trimesh::viewer_mainwindow gvmw(gdm);

      gvmw.setWindowTitle("ms complex vis");

      gvmw.show();

      application.exec();
  }
  else
  {
    delete gdm;
  }
}

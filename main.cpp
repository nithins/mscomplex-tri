#include <exception>
#include <string>

#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

#include <cpputils.h>

#include <boost/program_options.hpp>

#include <stdexcept>
#include <iostream>

using namespace std;

namespace bpo = boost::program_options ;

int main(int ac , char **av)
{
  string filename;

  n_vector_t<int,2> dim;

  bool   single_thread = false;

  bool   out_of_core_flag = false;

  uint   num_levels  = 1;

  double   simp_tresh= 0.0;

  uint   num_parallel  = 1;

  bool   gui = false;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<std::string >(), "grid file name")
      ("dim,d", bpo::value<n_vector_t<int,2> >(), "dim of grid entered as (x,y)")
      ("single-thread-mode,s", "single threaded mode")
      ("out-of-core-mode,o", "Compute out of Core")
      ("num-levels,n",bpo::value<int>(),"num levels to partition into")
      ("simp-tresh,t",bpo::value<double>(),"simplification treshold")
      ("num-parallel,p",bpo::value<int>(),"num subdomains to process in parallel")
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

  if (vm.count("dim"))
    dim = vm["dim"].as<n_vector_t<int,2> >();
  else
    throw std::invalid_argument("no dim specified");

  if (vm.count("file"))
    filename = vm["file"].as<std::string>();
  else
    throw std::invalid_argument("no filename specified");

  if (vm.count("single-thread-mode"))
    single_thread = true;

  if (vm.count("out-of-core-mode"))
    out_of_core_flag = true;

  if (vm.count("num-levels"))
    num_levels = vm["num-levels"].as<int>();

  if (vm.count("simp-tresh"))
    simp_tresh = vm["simp-tresh"].as<double>();

  if (vm.count("num-parallel"))
    num_parallel = vm["num-parallel"].as<int>();

  if (vm.count("gui"))
    gui = true;

  grid::data_manager_t * gdm = new grid::data_manager_t
      (filename,dim,
       num_levels,
       single_thread,
       simp_tresh,
       out_of_core_flag,
       num_parallel);

  if(gui)
  {
      QApplication application(ac,av);

      grid::viewer_mainwindow gvmw(gdm);

      gvmw.setWindowTitle("ms complex vis");

      gvmw.show();

      application.exec();
  }
  else
  {
    delete gdm;
  }
}

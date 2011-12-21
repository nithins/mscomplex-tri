#include <stdexcept>
#include <iostream>
#include <exception>
#include <string>

#include <trimesh_viewer_mainwindow.h>

#include <boost/program_options.hpp>


using namespace std;
namespace bpo = boost::program_options;

int main(int ac , char **av)
{
  string tri_filename;
  string graph_file;
  string mfold_file;

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("tri-file,t",bpo::value(&tri_filename)->required(), "tri file name")
      ("graph-bin-file,g",bpo::value(&graph_file)->required(), "ms graph file")
      ("mfold-bin-file,m",bpo::value(&mfold_file)->required(), "ms mfold file")
      ;


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

  QApplication application(ac,av);

  trimesh::viewer_mainwindow gvmw(tri_filename,graph_file,mfold_file);

  gvmw.setWindowTitle("ms complex vis");

  gvmw.show();

  application.exec();

}

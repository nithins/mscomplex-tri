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
  QApplication application(ac,av);

  trimesh::viewer_mainwindow gvmw;

  gvmw.setWindowTitle("ms complex vis");

  gvmw.show();

  application.exec();

}

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
  bpo::options_description desc("Allowed options");

  std::string script_file;

  desc.add_options()
      ("help,h", "produce help message")
      ("eval-script-file,e",bpo::value(&script_file)->default_value(""),
       "a script to evaluate")
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

  trimesh::viewer_mainwindow gvmw;

  gvmw.setWindowTitle("ms complex vis");

  gvmw.show();

  if(!script_file.empty())
    gvmw.eval_script(script_file.c_str());

  application.exec();
}

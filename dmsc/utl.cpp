#include "utl.h"

#include <sstream>
#include <fstream>

#include <vector>

#include <tr1/functional>


/*===========================================================================*/

using namespace std;

namespace utl {
/*---------------------------------------------------------------------------*/

void trim(std::string &s)
{
  s.erase(s.begin(), std::find_if
          (s.begin(), s.end(),
           std::not1(std::ptr_fun<int, int>(std::isspace))));
  s.erase(std::find_if
          (s.rbegin(), s.rend(),
           std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

/*---------------------------------------------------------------------------*/

file_line_iterator::file_line_iterator(const char * f,char c_char):
  is(new std::ifstream(f)),c_char(c_char)
{
  ENSUREV(is->is_open(),"cannot read the file",f);
  increment();
}

/*---------------------------------------------------------------------------*/

void file_line_iterator::increment()
{
  value.clear();

  while(is && value.size() == 0)
  {
    if(is->eof())
    {
      is.reset();
      value.clear();
      break;
    }

    std::getline(*is,value);

    int p = value.find(c_char);

    if ( p < value.size() )
      value = value.substr(0,p);

    trim(value);
  }
}

/*---------------------------------------------------------------------------*/

bool file_line_iterator::equal(file_line_iterator const& other) const
{
  if(!is && !other.is)
    return true;

  if(!is || !other.is)
    return false;

  ENSURE(is == other.is, "cannot compare distinct istream iters");

  return is->tellg() == other.is->tellg();
}

/*---------------------------------------------------------------------------*/

const std::string & file_line_iterator::dereference() const
{
  ENSURE(is,"dereferencing past end of line stream");
  return value;
}

/*---------------------------------------------------------------------------*/

boost::mutex logger::s_mutex;
logger       logger::s_logger;

/*---------------------------------------------------------------------------*/

} // namespace utl

/*===========================================================================*/


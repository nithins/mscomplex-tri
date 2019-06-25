#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <iostream>


#include <tri_edge.h>
#include <trimesh_dataset.h>
#include <trimesh_mscomplex.h>

using namespace std;
using namespace boost::python;
using namespace trimesh;
using namespace utl;

namespace bp = boost::python;


/*****************************************************************************/
namespace pymstri {

void read_off(const std::string &fname,
	      tri_cc_geom_t::vertex_list_t &verts,
	      tri_cc_geom_t::tri_idx_list_t &tris)
{
  file_line_iterator lgen(fname.c_str()),lend;

  ENSURES(*lgen++ == "OFF");

  int i = 0, j = 0 , nv = 0 , nt = 0 , np = 0;

  from_string(*lgen++,nv,nt);

  verts.resize(nv);
  tris.resize(nt);

  for(i = 0 ; i < nv && lgen != lend ; ++i)
    from_string(*lgen++,verts[i][0],verts[i][1],verts[i][2]);

  ENSURES(i == nv);

  for(j = 0 ; j < nt && lgen != lend ; ++j)
  {
    from_string(*lgen++,np,tris[j][0],tris[j][1],tris[j][2]);

    ENSURES(np == 3);
  }

  ENSURES(j == nt);
}

void read_off(const std::string &fname,
	      tri_cc_geom_t::vertex_list_t &verts,
	      tri_cc_geom_t::tri_idx_list_t &tris,
	      std::vector<double> &func)
{
  file_line_iterator lgen(fname.c_str()),lend;

  ENSURES(*lgen++ == "OFF");

  int i = 0, j = 0 , nv = 0 , nt = 0 , np = 0;

  from_string(*lgen++,nv,nt);

  verts.resize(nv);
  tris.resize(nt);
  func.resize(nv);

  for(i = 0 ; i < nv && lgen != lend ; ++i)
    from_string(*lgen++,verts[i][0],verts[i][1],verts[i][2],func[i]);

  ENSURES(i == nv);

  for(j = 0 ; j < nt && lgen != lend ; ++j)
  {
    from_string(*lgen++,np,tris[j][0],tris[j][1],tris[j][2]);

    ENSURES(np == 3);
  }

  ENSURES(j == nt);
}


/*****************************************************************************/
/******** Python wrapped Morse-Smale complex class                    ********/
/*****************************************************************************/

// Wrapper to hold on to references to some other objects
class mscomplex_pymstri_t: public mscomplex_t
{
public:
  tri_cc_geom_ptr_t tcc;
  dataset_ptr_t     ds;
};

typedef  boost::shared_ptr<mscomplex_pymstri_t> mscomplex_pymstri_ptr_t;

mscomplex_pymstri_ptr_t new_msc()
{
  mscomplex_pymstri_ptr_t msc(new mscomplex_pymstri_t);
  return msc;
}

tri_cc_geom_ptr_t mscomplex_get_tri_cc(mscomplex_pymstri_ptr_t msc)
{
  ENSURES(msc->tcc != 0);
  return msc->tcc;
}

void __mscomplex_compute_internel__(
    mscomplex_pymstri_ptr_t msc,
    tri_idx_list_t &tris,
    tri_cc_geom_t::vertex_list_t  &verts,
    fn_list_t &func)
{
  msc->tcc.reset(new tri_cc_geom_t());
  msc->tcc->init(tris,verts);

  verts.clear();
  tris.clear();

  msc->ds.reset(new dataset_t(func,msc->tcc->get_tri_cc()));
  msc->ds->work(msc);
}

void mscomplex_compute_off(mscomplex_pymstri_ptr_t msc,std::string off_file)
{
  tri_idx_list_t tris;
  tri_cc_geom_t::vertex_list_t  verts;
  fn_list_t func;

  read_off(off_file,verts,tris,func);
  __mscomplex_compute_internel__(msc,tris,verts,func);
}

template<typename T> void read_bin_arr
(std::string fn,std::vector<fn_t> & funcs,int n)
{
  fstream is(fn.c_str(),std::ios::binary|std::ios::in);

  funcs.resize(n);
  is.read((char*)(void*)(funcs.data()),n*sizeof(T));

  ENSURES(is) << "Ran out of data when reading file=" << fn <<"\n";
}


void mscomplex_compute_off_bin
(mscomplex_pymstri_ptr_t msc,std::string off_file,
 std::string bin_file,std::string bin_fmt="float64")
{
  tri_idx_list_t tris;
  tri_cc_geom_t::vertex_list_t  verts;
  fn_list_t funcs;

  read_off(off_file,verts,tris);

  if(bin_fmt == "float32")
    read_bin_arr<float>(bin_file,funcs,verts.size());
  else if(bin_fmt == "float64")
    read_bin_arr<double>(bin_file,funcs,verts.size());
  else
    throw std::runtime_error("Unknown bin format");

  __mscomplex_compute_internel__(msc,tris,verts,funcs);
}


int mscomplex_num_canc(mscomplex_pymstri_ptr_t msc)
{
  return msc->m_canc_list.size();
}

bp::tuple mscomplex_canc(mscomplex_pymstri_ptr_t msc,int i)
{
  ASSERT(is_in_range(i,0,msc->m_canc_list.size()));
  return bp::make_tuple(msc->m_canc_list[i].first,msc->m_canc_list[i].second);
}

bp::tuple mscomplex_frange(mscomplex_pymstri_ptr_t msc)
{
  return bp::make_tuple(msc->m_fmin,msc->m_fmax);
}

template <eGDIR dir>
bp::list mscomplex_conn(mscomplex_pymstri_ptr_t msc, int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(int c,msc->m_conn[dir][i])
  {
    r.append(c);
  }
  return r;
}

bp::list mscomplex_cps(mscomplex_pymstri_ptr_t msc)
{
  bp::list r;

  for(int i = 0 ; i < msc->get_num_critpts(); ++i)
    if(msc->is_not_paired(i))
      r.append(i);

  return r;
}

void mscomplex_gen_pers_pairs(mscomplex_pymstri_ptr_t msc)
{
  int lastVer = msc->get_multires_version();
  msc->simplify(1.0,true);
  msc->set_multires_version(lastVer);
}

void mscomplex_collect_mfolds(mscomplex_pymstri_ptr_t msc)
{
  ENSURES(msc->ds !=0)
      << "Gradient information unavailable" <<endl
      << "Did you load the mscomplex from a file!!!"<<endl;

  msc->collect_mfolds(msc->ds);
}

template <eGDIR dir>
bp::list mscomplex_geom(mscomplex_pymstri_ptr_t msc, int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(int c,msc->m_mfolds[dir][i])
  {
    r.append(c);
  }
  return r;
}

//bp::list mscomplex_arc_geom(mscomplex_ptr_t msc, int a, int b)
//{
//  bp::list r;

//  if(msc->m_arc_geom.count(int_pair_t(a,b)) != 0)
//  {
//      BOOST_FOREACH(int c,msc->m_arc_geom[int_pair_t(a,b)])
//      {
//        r.append(c);
//      }
//    }

//  return r;
//}


void wrap_mscomplex_t()
{

  docstring_options local_docstring_options(true, false, false);


  class_<mscomplex_pymstri_t,mscomplex_pymstri_ptr_t>
      ("mscomplex","The Morse-Smale complex object",no_init)
      .def("__init__", make_constructor( &new_msc),
           "ctor")
      .def("num_cp",&mscomplex_t::get_num_critpts,
           "Number of Critical Points")
      .def("fn",&mscomplex_t::fn,
           "Function value at critical point i")
      .def("index",&mscomplex_t::index,
           "Morse index od critical point i")
      .def("pair_idx",&mscomplex_t::pair_idx,
           "Index of the cp that is paired with i (-1 if it is not paired)")
      .def("is_boundry",&mscomplex_t::is_boundry,
           "If the cp is on the boundary or not")
      .def("vertid",&mscomplex_t::vertid,
           "vertex id of maximal vertex of critical cell")
      .def("cellid",&mscomplex_t::cellid,
           "cell id of critical cell")
      .def("load",&mscomplex_t::load,
           "Load mscomplex from file")
      .def("save",&mscomplex_t::save,
           "Save mscomplex to file")
      .def("num_canc",&mscomplex_num_canc,
           "Number of cancellation pairs")
      .def("canc",&mscomplex_canc,
           "The ith cancellation pair")
      .def("frange",&mscomplex_frange,
           "Range of function values")
      .def("asc",&mscomplex_conn<ASC>,
           "List of ascending cps connected to a given critical point i")
      .def("des",&mscomplex_conn<DES>,
           "List of descending cps connected to a given critical point i")
      .def("asc_geom",&mscomplex_geom<ASC>,
           "Ascending manifold geometry of a given critical point i")
      .def("des_geom",&mscomplex_geom<DES>,
           "Descending manifold geometry of a given critical point i")
      .def("cps",&mscomplex_cps,
           "Returns a list of surviving critical cps")
      .def("gen_pers_hierarchy",&mscomplex_gen_pers_pairs,
           "Generates the persistence hierarchy using topo simplification")
      .def("get_tri_cc",mscomplex_get_tri_cc,
           "Get the underlying triangulation object")
      .def("compute_off",&mscomplex_compute_off,
           "Compute the Mscomplex from a triagulation given in the .off format\n"\
           "\n"\
           "Note:The scalar function is assumed to be provided as \n"\
           "     the fourth coordinateof each vertex in the off file.\n"\
           "\n"\
           "Note: This only computes the combinatorial structure\n"\
           "     Call collect_mfold(s) to extract geometry\n"
           )
      .def("compute_off_bin",&mscomplex_compute_off_bin,
           "Compute the Mscomplex from a triagulation given in the .off format\n"\
           "and scalar function from the given bin file in the given bin format\n"\
           "\n"\
           "Parameters: \n"\
           "    off_file: the off file containing the triangulation.\n"\
           "    bin_file: the bin file containing the scalar function.\n"\
           "    bin_fmt: binary format .\n"\
           "             Acceptable values = (\"float32\",\"float64\")\n"
           "\n"\
           "Note: This only computes the combinatorial structure\n"\
           "     Call collect_mfold(s) to extract geometry\n"
           )
      .def("collect_geom",&mscomplex_collect_mfolds,
           "Collect the geometry of all survivng critical points\n"\
           "\n"\
           "Note: This must be called only after any of the compute functions are called. \n"\
           )

      .def("simplify_pers",&mscomplex_t::simplify,
           "Simplify the Morse-Smale complex using topological persistence\n"\
           "Parameters:\n"\
           "    tresh: persistence threshold\n"\
           "    is_nrm: is the threshold normalized to [0,1] or not.\n"\
           "            if not then thresh is in same scale as input function\n"\
           "    req_nmax,req_nmin: num maxima/minima that should be retained\n"\
           "                       set to 0 to ignore"
           )


//      .def("arc_geom",&mscomplex_arc_geom,
//           "seqence of cells of arc connecting 2 cps (empty list if no arc exists)")
      ;
}

/*****************************************************************************/
/******** Python wrapped Tri Cell Complex class                       ********/
/*****************************************************************************/



tri_cc_geom_ptr_t tcc_from_off(std::string f)
{
  tri_idx_list_t tris;
  tri_cc_geom_t::vertex_list_t  verts;

  read_off(f,verts,tris);
  tri_cc_geom_ptr_t ret(new tri_cc_geom_t);
  ret->init(tris,verts);
  return ret;
}

inline bp::tuple dcell_range(tri_cc_geom_ptr_t tcc, int d)
{
  ASSERT(is_in_range(d,0,3));

  int b = 0;
  int e = 0;

  switch(d)
  {
    case 2:
      b += tcc->get_tri_cc()->edge_ct();
      e += tcc->get_tri_cc()->tri_ct();

    case 1:
      b += tcc->get_tri_cc()->vert_ct();
      e += tcc->get_tri_cc()->edge_ct();

    case 0:
      e += tcc->get_tri_cc()->vert_ct();
      break;

    default:
      throw std::runtime_error("unrecogized dim");
  };

  return bp::make_tuple(b,e);
}

inline bp::tuple tri_cc_vrts(tri_cc_geom_ptr_t tcc, cellid_t i)
{
  cellid_t pts[4];

  tcc->get_cell_points(i,pts);

  switch(tcc->get_cell_dim(i))
    {
    case 0:
      return bp::make_tuple(i);
    case 1:
      return bp::make_tuple(pts[0],pts[1]);
    case 2:
      return bp::make_tuple(pts[0],pts[1],pts[2]);
    default:
      throw std::runtime_error("unrecogized dim");

    }

  return bp::make_tuple();
}

struct vert_to_tup
{
  static PyObject* convert(tri_cc_geom_t::vertex_t v)
  {
    return boost::python::incref(bp::make_tuple(v[0],v[1],v[2]).ptr());
  }
};


void wrap_tet_cc_t()
{
  boost::python::to_python_converter<tri_cc_geom_t::vertex_t,vert_to_tup>();

  class_<tri_cc_geom_t,tri_cc_geom_ptr_t >("tri_cc_geom",no_init)
      .def("__init__", make_constructor( &tcc_from_off) )
      .def("num_dcells",&tri_cc_geom_t::get_num_cells_dim)
      .def("num_cells",&tri_cc_geom_t::get_num_cells)
      .def("dim",&tri_cc_geom_t::get_cell_dim)
      .def("dcell_range",&dcell_range)
//      .def("get_tris_tetvers",get_tris_tetvers)
      .def("point",&tri_cc_geom_t::get_cell_position)
      .def("cell_vert",&tri_cc_vrts)
//      .def("centroid",&tet_geom_cc_t::get_centroid<0>)
      ;

}

/*****************************************************************************/
/******** Define the pymstet module                                   ********/
/*****************************************************************************/

BOOST_PYTHON_MODULE(pymstri)
{
  #if BOOST_VERSION < 106500
        numeric::array::set_module_and_type("numpy", "ndarray");
  #endif

  wrap_tet_cc_t();

  wrap_mscomplex_t();

}

}/****************************************************************************/


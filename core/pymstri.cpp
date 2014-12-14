#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>


#include <tri_edge.h>
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

/*****************************************************************************/
/******** Python wrapped Morse-Smale complex class                    ********/
/*****************************************************************************/

mscomplex_ptr_t read_msc(std::string f)
{
  mscomplex_ptr_t msc(new mscomplex_t);

  msc->load(f);

  return msc;
}

int mscomplex_num_canc(mscomplex_ptr_t msc)
{
  return msc->m_canc_list.size();
}

bp::tuple mscomplex_canc(mscomplex_ptr_t msc,int i)
{
  ASSERT(is_in_range(i,0,msc->m_canc_list.size()));
  return bp::make_tuple(msc->m_canc_list[i].first,msc->m_canc_list[i].second);
}

bp::tuple mscomplex_frange(mscomplex_ptr_t msc)
{
  return bp::make_tuple(msc->m_fmin,msc->m_fmax);
}

bp::list mscomplex_asc(mscomplex_ptr_t msc,int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(int c,msc->m_asc_conn[i])
  {
    r.append(c);
  }
  return r;
}

bp::list mscomplex_des(mscomplex_ptr_t msc,int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(int c,msc->m_des_conn[i])
  {
    r.append(c);
  }
  return r;
}

bp::list mscomplex_cps(mscomplex_ptr_t msc)
{
  bp::list r;

  for(int i = 0 ; i < msc->get_num_critpts(); ++i)
    if(msc->is_not_paired(i))
      r.append(i);

  return r;
}

void mscomplex_gen_pers_pairs(mscomplex_ptr_t msc)
{
  int lastVer = msc->get_multires_version();
  msc->simplify(1.0,true);
  msc->set_multires_version(lastVer);
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


  class_<mscomplex_t,mscomplex_ptr_t>("mscomplex","The Morse-Smale complex object",no_init)
      .def("__init__", make_constructor( &read_msc),
           "build the MSC object using the graph.bin file")
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
      .def("save",&mscomplex_t::save,
           "Save mscomplex to file")
      .def("num_canc",&mscomplex_num_canc,
           "Number of cancellation pairs")
      .def("canc",&mscomplex_canc,
           "The ith cancellation pair")
      .def("frange",&mscomplex_frange,
           "Range of function values")
      .def("asc",&mscomplex_asc,
           "List of ascending cps connected to i")
      .def("des",&mscomplex_des,
           "List of descending cps connected to i")
      .def("cps",&mscomplex_cps,
           "Returns a list of surviving critical cps")
      .def("gen_pers_hierarchy",&mscomplex_gen_pers_pairs,
           "Generates the persistence hierarchy using topo simplification")

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

//numeric::array get_tris_tetvers(tet_geom_cc_ptr_t tcc,
//                   const numeric::array &tris)
//{
//  const bp::tuple & shape = extract<bp::tuple>(tris.attr("shape"));

//  int ntris = extract<int>(shape[0]);

//  bp::list lst;

//  for(int i = 0 ; i < ntris; ++i)
//  {
//    object ti_obj = tris[bp::make_tuple(i,0)];

//    int tri = boost::python::extract<int>(ti_obj.attr("__int__")());

//    ASSERT(is_in_range(tri,tcc->m_cel_off[2],tcc->m_cel_off[3]));

//    int d1 = tcc->m_cel_drt[tri];
//    int d2 = tcc->beta2(d1);

//    int v1 = (d1%12)/3 , v2 = (d2%12)/3;

//    d1 /=24;
//    if (d2 != -1)
//      d2/=24;

//    lst.append(bp::make_tuple(d1,v1,d2,v2));
//  }

//  return numeric::array(lst);
//}

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
  numeric::array::set_module_and_type("numpy", "ndarray");

  wrap_tet_cc_t();

  wrap_mscomplex_t();

}

}/****************************************************************************/


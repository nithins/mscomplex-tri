#ifndef MSCOMPLEX_SIMP_INCLUDED
#define MSCOMPLEX_SIMP_INCLUDED

#include <cmath>
#include <queue>
#include <limits>

#include <boost/foreach.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/range/adaptors.hpp>

#include <trimesh_mscomplex.h>
#include <trimesh_dataset.h>

namespace trimesh
{

template <eGDIR dir>
inline double get_area(const tri_cc_geom_t &tcc, cellid_t c);

template <> double inline get_area<ASC>(const tri_cc_geom_t &tcc, cellid_t c)
{return tcc.get_vert_area(c);}

template <> double inline get_area<DES>(const tri_cc_geom_t &tcc, cellid_t c)
{return tcc.get_tri_area(c);}

template <eGDIR dir>
double get_harea(const tri_cc_geom_t &tcc, const dataset_ptr_t &ds ,  cellid_t c);

template <> double get_harea<ASC>(const tri_cc_geom_t &tcc,
                                  const dataset_ptr_t &ds ,
                                  cellid_t cid)
{
  cellid_t st[40];

  uint st_ct = tcc.get_vert_star(cid,st);

  double area = 0;

  la::dvec4_t b = la::make_vec<double>
      (tcc.get_cell_position(cid),ds->fn<dataset_t::CFI_AVE>(cid));

  for(uint i = 1; i < st_ct; i++)
  {
    la::dvec4_t a = la::make_vec<double>
        (tcc.get_cell_position(st[i-1]),ds->fn<dataset_t::CFI_AVE>(st[i-1]));

    la::dvec4_t c = la::make_vec<double>
        (tcc.get_cell_position(st[i]),ds->fn<dataset_t::CFI_AVE>(st[i]));

    area += la::tri_area<double,4>(a,b,c);
  }

  return area;
}

template <> double get_harea<DES>(const tri_cc_geom_t &tcc,
                                  const dataset_ptr_t &ds ,
                                  cellid_t id)
{
  cellid_t pts[20];

  tcc.get_cell_points(id,pts);

  la::dvec4_t a = la::make_vec<double>
      (tcc.get_cell_position(pts[0]),ds->fn(pts[0]));

  la::dvec4_t b = la::make_vec<double>
      (tcc.get_cell_position(pts[1]),ds->fn(pts[1]));

  la::dvec4_t c = la::make_vec<double>
      (tcc.get_cell_position(pts[2]),ds->fn(pts[2]));

  return la::tri_area<double,4>(a,b,c);
}

template <eGDIR dir>
double get_mfold_area
  (mscomplex_ptr_t msc,
   tri_cc_geom_ptr_t tcc,
   const int_pair_t &pr)
{
  const int EIDX = (dir == DES)? (0):(1);

  mfold_t &mfold = msc->mfold<dir>(pr[EIDX]);

  double area = 0;

  BOOST_FOREACH(cellid_t c,mfold)
  {
    area += get_area<dir>(*tcc,c);
  }

  return area;
}

enum eApCancType {AP_AREA_WEIGHTED_PERSISTENCE,
                  AP_AREA_BEFORE_PERSISTENCE};

template <eApCancType apct>
struct ap_edge
{
  int_pair_t   edge;
  double       pers;
  double       area;

  ap_edge(int_pair_t e,double a,double p):edge(e),area(a),pers(p){}

  static double s_get_val(double pers,double area);

  double get_val() const
  {return s_get_val(pers,area);}

  bool operator < (const ap_edge & ap_e) const;

};

template<>
double ap_edge<AP_AREA_WEIGHTED_PERSISTENCE>::s_get_val(double pers, double area)
{return pow(pers*area,1.0/10.0);}

template<>
bool ap_edge<AP_AREA_WEIGHTED_PERSISTENCE>::operator < (const ap_edge & ap_e) const
{return ap_e.get_val() < get_val() ;}



template<>
double ap_edge<AP_AREA_BEFORE_PERSISTENCE>::s_get_val(double pers, double area)
{return pow(area,1.0/5.0);}

template<>
bool ap_edge<AP_AREA_BEFORE_PERSISTENCE>::operator < (const ap_edge & ap_e) const
{
  if(ap_e.area != area)
    return ap_e.area < area;
  else
    return ap_e.pers < pers;
}


template<eApCancType apct>
void simplify_ap(mscomplex_ptr_t msc, tri_cc_geom_ptr_t tcc,double tresh)
{
  typedef ap_edge<apct> edge_t;

  using namespace std;

  namespace br = boost::range;

  priority_queue<edge_t> pq;
  vector<double>  ex_mfold_area(msc->get_num_critpts());
  double total_area=0;
  double frange  = (*br::max_element(msc->m_cp_fn) - *br::min_element(msc->m_cp_fn));

  for( int i = 0 ; i < msc->get_num_critpts(); ++i)
  {
    BOOST_FOREACH(int j,msc->m_conn[0][i])
    {
      int_pair_t e = la::make_vec<int>(i,j);

      double area = 0;
      int ex_idx=-1;

      if (msc->index(i) == 2)
      {
        area = get_mfold_area<DES>(msc,tcc,e);
        ex_idx = i;
      }
      else
      {
        area = get_mfold_area<ASC>(msc,tcc,e);
        ex_idx = j;
      }

      double pers = abs<double>(msc->fn(i) -msc->fn(j));

      pq.push(edge_t(e,area,pers));

      ex_mfold_area[ex_idx] = area;

      total_area += area;
    }
  }

  while (pq.size() != 0 )
  {
    edge_t ap_e = pq.top(); pq.pop();

    int ex = ap_e.edge[0],sd = ap_e.edge[1];
    if(msc->index(ex) == 1) swap(ex,sd);

    if(is_valid_canc_edge(*msc,ap_e.edge) == false)
      continue;

    if(ap_e.area != ex_mfold_area[ex])
    {
      pq.push(edge_t(ap_e.edge,ex_mfold_area[ex],ap_e.pers));
      continue;
    }

    if(ap_e.get_val() > edge_t::s_get_val(frange,total_area)*tresh)
      break;

    msc->cancel_pair(ap_e.edge);
    msc->m_canc_list.push_back(ap_e.edge);
    msc->m_canc_pers.push_back(ap_e.get_val()/edge_t::s_get_val(frange,total_area));

    int surv_ex = *msc->m_conn[((msc->index(ex) == 2)?(ASC):(DES))][sd].begin();
    ex_mfold_area[surv_ex] += ex_mfold_area[ex];
    ex_mfold_area[ex]       = 0;

    BOOST_FOREACH(int i,msc->m_des_conn[ap_e.edge[0]])
    {
      BOOST_FOREACH(int j,msc->m_asc_conn[ap_e.edge[1]])
      {
        int_pair_t e = la::make_vec<int>(i,j);

        int e_ex = i,e_sd = j;
        if(msc->index(e_ex) == 1) swap(e_ex,e_sd);

        if(is_valid_canc_edge(*msc,e))
          pq.push(edge_t(e,ex_mfold_area[e_ex],abs<double>(msc->fn(i) - msc->fn(j))));
      }
    }
  }

  for(int i = 0 ; i < msc->m_canc_list.size() ; ++i)
  {
    cout<<"canc -->"
        <<msc->m_canc_list[i].transpose()<<" :: "
        <<msc->m_canc_pers[i]<<endl;
  }

}

inline void simplify_awp(mscomplex_ptr_t msc,tri_cc_geom_ptr_t tcc,double tresh)
{simplify_ap<AP_AREA_WEIGHTED_PERSISTENCE>(msc,tcc,tresh);}

inline void simplify_abp(mscomplex_ptr_t msc,tri_cc_geom_ptr_t tcc,double tresh)
{simplify_ap<AP_AREA_BEFORE_PERSISTENCE>(msc,tcc,tresh);}


}
#endif //MSCOMPLEX_SIMP_INCLUDED

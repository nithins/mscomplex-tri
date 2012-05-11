#include <cmath>
#include <queue>
#include <limits>

#include <boost/foreach.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/range/adaptors.hpp>

#include <trimesh_mscomplex.h>
#include <trimesh_dataset.h>

using namespace std;
namespace br    = boost::range;
namespace badpt = boost::adaptors;

namespace trimesh
{
inline std::string edge_to_string(mscomplex_t *msc,int_pair_t e)
{
  std::stringstream ss;

  ss<<utls::to_string(msc->cellid(e[0]))<<"----"<<utls::to_string(msc->cellid(e[0]));

  return ss.str();
}

mscomplex_t::mscomplex_t()
  :m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),
    m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]){}

mscomplex_t::~mscomplex_t(){clear();}

void mscomplex_t::set_critpt(int i, cellid_t c, char idx, fn_t f, cellid_t v, bool b)
{
  m_cp_cellid[i]     = c;
  m_cp_vertid[i]     = v;
  m_cp_index[i]      = idx;
  m_cp_fn[i]         = f;
  m_cp_is_boundry[i] = b;
}

void  mscomplex_t::resize(int i)
{
  m_cp_cellid.resize(i,invalid_cellid);
  m_cp_vertid.resize(i,invalid_cellid);
  m_cp_index.resize(i,-1);
  m_cp_pair_idx.resize(i,-1);
  m_cp_is_boundry.resize(i,false);
  m_cp_fn.resize(i);
  m_des_conn.resize(i);
  m_asc_conn.resize(i);

  m_des_mfolds.resize(i);
  m_asc_mfolds.resize(i);
}

void mscomplex_t::clear()
{
  m_cp_cellid.clear();
  m_cp_vertid.clear();
  m_cp_pair_idx.clear();
  m_cp_index.clear();
  m_cp_is_boundry.clear();
  m_cp_fn.clear();
  m_des_conn.clear();
  m_asc_conn.clear();
  m_des_mfolds.clear();
  m_asc_mfolds.clear();
}

void mscomplex_t::connect_cps(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);

  // if a d-cp hits a d+-1 cp and the d+-1 cp is paired
  // then the connection is useful iff the dimension of the pair is d

  ASSERT(!(is_paired(p) && index(pair_idx(p))!= index(q)));
  ASSERT(!(is_paired(q) && index(pair_idx(q))!= index(p)));
  ASSERT(m_des_conn[p].count(q) == m_asc_conn[q].count(p));

  if(m_des_conn[p].count(q) == 2)
    return;

  m_des_conn[p].insert(q);
  m_asc_conn[q].insert(p);
}

void mscomplex_t::dir_connect_cps(int p, int q)
{
  ASSERT(is_paired(p) != is_paired(q));
  ASSERT(abs(index(p)-index(q)) == 1);

  if(is_paired(q))
    std::swap(p,q);

  conn_t &conn = (index(p) > index(q))?(m_des_conn[p]):(m_asc_conn[p]);

  if(conn.count(q) == 0)
    conn.insert(q);
}

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q);
  ASSERT(pair_idx(q) == p);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);

  conn_iter_t i,j;

  m_des_conn[p].erase(q);
  m_asc_conn[q].erase(p);

  // cps in lower of u except l
  for(i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
    for(j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
    {
      ASSERT(is_paired(*i) == false);
      ASSERT(is_paired(*j) == false);

      connect_cps(*i,*j);
    }

  for(j = m_des_conn[p].begin();j != m_des_conn[p].end();++j)
    m_asc_conn[*j].erase(p);

  for(j = m_asc_conn[p].begin();j != m_asc_conn[p].end();++j)
    m_des_conn[*j].erase(p);

  for(j = m_des_conn[q].begin();j != m_des_conn[q].end();++j)
    m_asc_conn[*j].erase(q);

  for(j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
    m_des_conn[*j].erase(q);

  m_asc_conn[p].clear();
  m_des_conn[q].clear();
}

void mscomplex_t::uncancel_pair(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q && pair_idx(q) == p);

  conn_iter_t i,j;

  for(int d = 0 ; d <2 ; ++d)
  {
    int ed = (d == 0)?(p):(q);

    conn_t old_conn(m_conn[d][ed].begin(),m_conn[d][ed].end());

    m_conn[d][ed].clear();

    for(i = old_conn.begin();i != old_conn.end() ; ++i)
    {
      if(is_paired(*i) == false)
      {
        dir_connect_cps(ed,*i);
        continue;
      }

      int r = pair_idx(*i);

      if(index(ed) != index(r))
        continue;

      ASSERT(abs(index(*i) - index(r)) == 1);
      ASSERT(pair_idx(r) == int(*i) && pair_idx(*i) ==  r);

      for(j = m_conn[d][r].begin(); j!= m_conn[d][r].end() ; ++j )
        dir_connect_cps(ed,*j);
    }
  }
}

inline bool is_valid_canc_edge(const mscomplex_t &msc,int_pair_t e )
{
  order_pr_by_cp_index(msc,e[0],e[1]);

  if(msc.is_paired(e[0]) || msc.is_paired(e[1]))
    return false;

  if(msc.is_boundry(e[0]) != msc.is_boundry(e[1]))
    return false;

  ASSERT(msc.m_des_conn[e[0]].count(e[1]) == msc.m_asc_conn[e[1]].count(e[0]));

  if(msc.m_des_conn[e[0]].count(e[1]) != 1)
    return false;

  return true;
}

inline bool is_epsilon_persistent(const mscomplex_t &msc,int_pair_t e )
{
  return (msc.vertid(e[0]) == msc.vertid(e[1]));
}

inline int get_num_new_edges(const mscomplex_t &msc, int_pair_t pr)
{
  order_pr_by_cp_index(msc,pr[0],pr[1]);

  int pd = msc.m_conn[DES][pr[0]].size();
  int pa = msc.m_conn[ASC][pr[0]].size();

  int qd = msc.m_conn[DES][pr[1]].size();
  int qa = msc.m_conn[ASC][pr[1]].size();

  return (pd - 1)*(qa - 1) - (pd + qd + pa +qa -1);
}

inline bool fast_persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
{
  bool eps_p0 = is_epsilon_persistent(msc,p0);
  bool eps_p1 = is_epsilon_persistent(msc,p1);

  if( eps_p0 != eps_p1)
    return eps_p0;

  return (get_num_new_edges(msc,p0) < get_num_new_edges(msc,p1));
}

inline fn_t get_persistence(const mscomplex_t & msc,int_pair_t e)
{
  return std::abs(msc.fn(e[0]) - msc.fn(e[1]));
}

inline bool persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
{/*
  bool eps_p0 = is_epsilon_persistent(msc,p0);
  bool eps_p1 = is_epsilon_persistent(msc,p1);

  if( eps_p0 != eps_p1)
    return eps_p0;

  fn_t d0 = get_persistence(msc,p0);
  fn_t d1 = get_persistence(msc,p1);

  if(d0 != d1)
    return d0 < d1;

  return p0 < p1;*/

  order_pr_by_cp_index(msc,p0[0],p0[1]);
  order_pr_by_cp_index(msc,p1[0],p1[1]);

  cellid_t v00 = msc.vertid(p0[0]);
  cellid_t v01 = msc.vertid(p0[1]);
  cellid_t v10 = msc.vertid(p1[0]);
  cellid_t v11 = msc.vertid(p1[1]);

  cellid_t c00 = msc.cellid(p0[0]);
  cellid_t c01 = msc.cellid(p0[1]);
  cellid_t c10 = msc.cellid(p1[0]);
  cellid_t c11 = msc.cellid(p1[1]);

  if( (v00 == v01 ) != (v10 == v11))
    return (v00 == v01 );

  if( (v00 == v01 ) &&(v10 == v11))
  {
    if(v00 == v10)
    {
      if(c00 != c10)
        return c00 < c10;
      else
        return c01 < c11;
    }
    else
    {
      return (v00 < v10);
    }
  }

  fn_t f00 = msc.fn(p0[0]);
  fn_t f01 = msc.fn(p0[1]);
  fn_t f10 = msc.fn(p1[0]);
  fn_t f11 = msc.fn(p1[1]);

  fn_t d1 = std::abs(f01-f00);
  fn_t d2 = std::abs(f11-f10);

  if(d1 != d2)
    return d1 < d2;

  if(c00 != c10)
    return c00 < c10;

  return c01 < c11;
}

inline bool is_within_treshold(const mscomplex_t & msc,int_pair_t e,fn_t t)
{
  return (is_epsilon_persistent(msc,e) || get_persistence(msc,e) < t);
}

void mscomplex_t::simplify(double f_tresh)
{
  BOOST_AUTO(cmp,bind(persistence_lt,boost::cref(*this),_2,_1));

  priority_queue<int_pair_t,int_pair_list_t,decltype(cmp)> pq(cmp);

  double f_range = *max_element(m_cp_fn.begin(),m_cp_fn.end()) -
      *min_element(m_cp_fn.begin(),m_cp_fn.end());

  f_tresh *= f_range;

  for(int i = 0 ;i < get_num_critpts();++i)
  {
    for(conn_iter_t j = m_des_conn[i].begin();j != m_des_conn[i].end() ;++j)
    {
      int_pair_t pr;

      pr[0] = i;
      pr[1] = *j;

      if(is_valid_canc_edge(*this,pr) && is_within_treshold(*this,pr,f_tresh))
        pq.push(pr);
    }
  }

  uint num_cancellations = 0;

  while (pq.size() !=0)
  {
    int_pair_t pr = pq.top();

    pq.pop();

    if(is_valid_canc_edge(*this,pr) == false)
      continue;

    int p = pr[0],q = pr[1];

    order_pr_by_cp_index(*this,p,q);

    pair_cps(p,q);

    cancel_pair(p,q);

    num_cancellations++;

    m_canc_list.push_back(pr);

    for(conn_iter_t i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
      for(conn_iter_t j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
      {
        int_pair_t npr;
        npr[0] = *i;
        npr[1] = *j;

        if(is_valid_canc_edge(*this,npr) && is_within_treshold(*this,npr,f_tresh))
          pq.push(npr);
      }
  }
}

void mscomplex_t::un_simplify()
{
  for(auto it = m_canc_list.rbegin();it != m_canc_list.rend() ; ++it)
  {
    uncancel_pair((*it)[0],(*it)[1]);
  }
}

void mscomplex_t::get_mfolds(dataset_ptr_t ds)
{
  vector<cellid_list_t> contrib_cells[GDIR_CT];

  contrib_cells[0].resize(get_num_critpts());
  contrib_cells[1].resize(get_num_critpts());

  for(int i = 0 ; i < get_num_critpts(); ++i)
  {
    if(is_paired(i) == false)
    {
      contrib_cells[DES][i].push_back(cellid(i));
      contrib_cells[ASC][i].push_back(cellid(i));
      continue;
    }

    ASSERT(is_paired(i) && is_paired(pair_idx(i))== true);
    ASSERT(abs(index(i)- index(pair_idx(i))) == 1);
    ASSERT(pair_idx(pair_idx(i))  == i);

    int dir = (index(i) > index(pair_idx(i)))?(0):(1);

    BOOST_FOREACH(int j, m_conn[dir][i])
    {
      if(is_paired(j) == false)
        contrib_cells[dir^1][j].push_back(cellid(pair_idx(i)));
    }
  }

  auto rng = cp_range()|
      badpt::filtered(bind(&mscomplex_t::is_not_paired,this,_1));

  for(auto b = boost::begin(rng); b != boost::end(rng); ++b)
  {
    if(index(*b)!=0) ds->get_mfold<DES>(mfold<DES>(*b),contrib_cells[DES][*b]);
    if(index(*b)!=2) ds->get_mfold<ASC>(mfold<ASC>(*b),contrib_cells[ASC][*b]);
  }
}

void mscomplex_t::clear_mfolds()
{
  m_mfolds[0].clear();
  m_mfolds[1].clear();

  m_mfolds[0].resize(get_num_critpts());
  m_mfolds[1].resize(get_num_critpts());
}

template<typename T>
inline void bin_write_vec(std::ostream &os, std::vector<T> &v)
{os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));}

template<typename T>
inline void bin_write(std::ostream &os, const T &v)
{os.write((const char*)(const void*)&v,sizeof(T));}

template<typename T>
inline void bin_read_vec(std::istream &is, std::vector<T> &v,int n)
{v.resize(n);is.read((char*)(void*)v.data(),n*sizeof(T));}

inline void bin_read_conn(std::istream &is,conn_t &conn,int n)
{
  vector<int> vec(n);
  is.read((char*)(void*)vec.data(),n*sizeof(int));
  conn.insert(vec.begin(),vec.end());
}


template<typename T>
inline void bin_read(std::istream &is, const T &v)
{is.read((char*)(void*)&v,sizeof(T));}

void mscomplex_t::save(std::ostream &os)
{
  int N = get_num_critpts();

  bin_write(os,N);

  bin_write_vec(os,m_cp_cellid);
  bin_write_vec(os,m_cp_vertid);
  bin_write_vec(os,m_cp_pair_idx);
  bin_write_vec(os,m_cp_index);
  bin_write_vec(os,m_cp_is_boundry);
  bin_write_vec(os,m_cp_fn);

  int_list_t nconn(2*N),nmfold(2*N);
  int_list_t adj;

  for(int i = 0 ; i < N; ++i)
  {
    nconn[2*i]     = m_des_conn[i].size();
    nconn[2*i+1]   = m_asc_conn[i].size();

    nmfold[2*i]    = m_des_mfolds[i].size();
    nmfold[2*i+1]  = m_asc_mfolds[i].size();

    br::copy(m_des_conn[i],back_inserter(adj));
    br::copy(m_asc_conn[i],back_inserter(adj));
  }

  bin_write_vec(os,nconn);
  bin_write_vec(os,adj);

  bin_write(os,int(m_canc_list.size()));
  bin_write_vec(os,m_canc_list);

  bin_write_vec(os,nmfold);

  for(int i = 0 ; i < N; ++i)
  {
    bin_write_vec(os,m_des_mfolds[i]);
    bin_write_vec(os,m_asc_mfolds[i]);
  }
}

void mscomplex_t::load(std::istream &is)
{
  clear();

  int N;
  int_list_t  nconn;
  int_list_t &nmfolds=nconn;

  bin_read(is,N);
  bin_read_vec(is,m_cp_cellid,N);
  bin_read_vec(is,m_cp_vertid,N);
  bin_read_vec(is,m_cp_pair_idx,N);
  bin_read_vec(is,m_cp_index,N);
  bin_read_vec(is,m_cp_is_boundry,N);
  bin_read_vec(is,m_cp_fn,N);

  bin_read_vec(is,nconn,2*N);
  m_des_conn.resize(N);
  m_asc_conn.resize(N);
  for(int i = 0 ; i < N; ++i)
  {
    bin_read_conn(is,m_des_conn[i],nconn[2*i]);
    bin_read_conn(is,m_asc_conn[i],nconn[2*i+1]);
  }


  int NC;
  bin_read(is,NC);
  bin_read_vec(is,m_canc_list,NC);


  bin_read_vec(is,nmfolds,2*N);
  m_des_mfolds.resize(N);
  m_asc_mfolds.resize(N);
  for(int i = 0 ; i < N; ++i)
  {
    bin_read_vec(is,m_des_mfolds[i],nmfolds[2*i]);
    bin_read_vec(is,m_asc_mfolds[i],nmfolds[2*i+1]);
  }

}

//template<typename T>
//inline void bin_write_vec(std::ostream &os, std::vector<T> &v)
//{os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));}

//template<typename T>
//inline void bin_write(std::ostream &os, const T &v)
//{os.write((const char*)(const void*)&v,sizeof(T));}

inline int get_header_size(int num_cps)
{
  return sizeof(int)              + // num_cps
         sizeof(cellid_t)*num_cps + // cellids
         sizeof(char)*num_cps     + // cell indices
         sizeof(int)*(2*num_cps);   // numcells
}

void write_header
(std::ostream & os,
 int_list_t & nmcells,
 cellid_list_t &cps,
 char_list_t &cp_inds,
 ios::off_type hoff= 0)
{
  os.seekp(hoff,ios::beg);

  bin_write(os,(int)cps.size());
  bin_write_vec(os,cps);
  bin_write_vec(os,cp_inds);
  bin_write_vec(os,nmcells);
}


void  mscomplex_t::save_mfolds(std::ostream &os,dataset_ptr_t ds)
{
  auto rng = cp_range()|badpt::filtered
      (boost::bind(&mscomplex_t::is_not_paired,this,_1));

  int num_cps = utls::count(boost::begin(rng),boost::end(rng));
  int hoff    = os.tellp();

  os.seekp(get_header_size(num_cps),ios::cur);

  int_list_t nmcells;

  for(auto i = boost::begin(rng);i!=boost::end(rng);++i)
  {
    int_list_t    desop,ascop;

    auto desop_bi = back_inserter(desop);
    auto ascop_bi = back_inserter(ascop);

    BOOST_FOREACH(cellid_t c,m_mfolds[0][*i]) ds->m_tcc.cellid_to_output(c,desop_bi);
    BOOST_FOREACH(cellid_t c,m_mfolds[1][*i]) ds->m_tcc.cellid_to_output(c,ascop_bi);

    nmcells.push_back(desop.size());
    nmcells.push_back(ascop.size());

    bin_write_vec(os,desop);
    bin_write_vec(os,ascop);
  }

  cellid_list_t cps;
  char_list_t   cp_idxs;

  br::copy(rng,back_inserter(cps));
  br::transform(rng,back_inserter(cp_idxs),bind(&mscomplex_t::index,this,_1));

  br::copy(cps,ostream_iterator<cellid_t>(cout," ")); cout<<endl;

  write_header(os,nmcells,cps,cp_idxs,hoff);
}


void mscomplex_t::save_ascii(const std::string &f)
{
  fstream os(f.c_str(),ios::out);

  os<<"# Num Cps"<<std::endl;

  os<<get_num_critpts()<<std::endl;

  os<<"# SL.No  cpIdx isPaired pair_idx vertNo fn "<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<i<<" ";
    os<<(int)index(i)<<" ";
    os<<(int)is_paired(i)<<" ";
    os<<(int)vertid(i)<<" ";
    os<<(fn_t)fn(i)<<" ";

    os<<std::endl;
  }

  os<<"#slno numDes numAsc connList"<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<(int)m_conn[0][i].size()<<" ";
    os<<(int)m_conn[1][i].size()<<" ";

    br::copy(m_conn[0][i],ostream_iterator<int>(os," "));
    br::copy(m_conn[1][i],ostream_iterator<int>(os," "));

    os<<std::endl;
  }

  os<<"#slno numDes numAsc mfoldCellids"<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<(int)m_mfolds[0][i].size()<<" ";
    os<<(int)m_mfolds[1][i].size()<<" ";

    br::copy(m_mfolds[0][i],ostream_iterator<int>(os," "));
    br::copy(m_mfolds[1][i],ostream_iterator<int>(os," "));

    os<<std::endl;
  }
}

}

#include <tuple>

namespace trimesh
{

typedef int_pair_t edge_t;

template <typename T>
inline edge_t mk_edge(const T & u,const T & v)
{edge_t e; e[0] = u;e[1] = v; return e;}


template <eGDIR dir>
double get_hyp_volume(mscomplex_ptr_t msc,dataset_ptr_t ds,const int_pair_t &pr)
{
  const int EIDX = (dir == DES)? (0):(1);

  mfold_t &mfold = msc->mfold<dir>(pr[EIDX]);

  cellid_t p[40];

  double sum = 0;

  BOOST_FOREACH(cellid_t c,mfold)
  {
    cellid_t *pb = p,*pe = p+ds->get_points<dir>(c,p);

    fn_t max_f = ds->fn(*pb);
    fn_t min_f = ds->fn(*pb);

    for (++pb ;pb != pe; ++pb)
    {
      max_f = max(max_f,ds->fn(*pb));
      min_f = min(min_f,ds->fn(*pb));
    }

    sum += (max_f - min_f);
  }

  return sum;
}

struct hv_edge
{
  edge_t     edge;
  double     hv_pers;

  hv_edge(edge_t e,double p):edge(e),hv_pers(p){}

  bool operator < (const hv_edge & hve) const
  { return hve.hv_pers < hv_pers ;}
};

void mscomplex_t::simplify_hypervolume(dataset_ptr_t ds, double tresh)
{
  priority_queue<hv_edge> pq;
  vector<double>  extrema_pers(get_num_critpts());
  double total_extrema_pers=0;

  for( int i = 0 ; i < get_num_critpts(); ++i)
  {
    BOOST_FOREACH(int j,m_conn[0][i])
    {
      edge_t e = mk_edge(i,j);

      double hv_pers = 0;
      int ex_idx=-1;

      if (index(i) == 2)
      {
        hv_pers = get_hyp_volume<DES>(shared_from_this(),ds,e);
        ex_idx = i;
      }
      else
      {
        hv_pers = get_hyp_volume<ASC>(shared_from_this(),ds,e);
        ex_idx = j;
      }

      pq.push(hv_edge(e,hv_pers));

      extrema_pers[ex_idx] = hv_pers;

      total_extrema_pers += hv_pers;
    }
  }

  while (pq.size() != 0 )
  {
    hv_edge hve = pq.top(); pq.pop();

    int ex = hve.edge[0],sd = hve.edge[1];
    if(index(ex) == 1) swap(ex,sd);

    if(is_valid_canc_edge(*this,hve.edge) == false)
      continue;

    if(hve.hv_pers != extrema_pers[ex])
    {
      pq.push(hv_edge(hve.edge,extrema_pers[ex]));
      continue;
    }

    if(hve.hv_pers > total_extrema_pers*tresh)
      break;

    pair_cps(hve.edge);
    cancel_pair(hve.edge[0],hve.edge[1]);
    m_canc_list.push_back(hve.edge);


    int surv_ex = *m_conn[((index(ex) == 2)?(ASC):(DES))][sd].begin();
    extrema_pers[surv_ex] += extrema_pers[ex];
    extrema_pers[ex]       = 0;

    BOOST_FOREACH(int i,m_des_conn[hve.edge[0]])
    {
      BOOST_FOREACH(int j,m_asc_conn[hve.edge[1]])
      {
        edge_t e = mk_edge(i,j);

        int e_ex = i,e_sd = j;
        if(index(e_ex) == 1) swap(e_ex,e_sd);

        if(is_valid_canc_edge(*this,e))
          pq.push(hv_edge(e,extrema_pers[e_ex]));
      }
    }
  }



}

}



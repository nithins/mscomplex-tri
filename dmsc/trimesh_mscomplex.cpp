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
    m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1])
{
  m_num_surv_dcps[0] = 0;
  m_num_surv_dcps[1] = 0;
  m_num_surv_dcps[2] = 0;
}

mscomplex_t::~mscomplex_t(){clear();}

void mscomplex_t::set_critpt(int i, cellid_t c, char idx, fn_t f, cellid_t v, bool b)
{
  m_cp_cellid[i]     = c;
  m_cp_vertid[i]     = v;
  m_cp_index[i]      = idx;
  m_cp_fn[i]         = f;
  m_cp_is_boundry[i] = b;

  ++m_num_surv_dcps[idx];
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

  m_num_surv_dcps[0] = 0;
  m_num_surv_dcps[1] = 0;
  m_num_surv_dcps[2] = 0;
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

void mscomplex_t::cancel_pair ( int_pair_t pr)
{
  int p = pr[0];
  int q = pr[1];

  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(is_paired(p) == false);
  ASSERT(is_paired(q) == false);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);

  pair_cps(p,q);

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

  --m_num_surv_dcps[index(p)];
  --m_num_surv_dcps[index(q)];
}

void mscomplex_t::uncancel_pair ( int_pair_t pr)
{
  int p = pr[0];
  int q = pr[1];

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

template<bool is_nrm> inline bool is_within_treshold
(const mscomplex_t & msc,int_pair_t e,fn_t t,fn_t r);

template<> inline bool is_within_treshold<true>
(const mscomplex_t & msc,int_pair_t e,fn_t t,fn_t r)
{
  // convention .. if threshold is in [0,1] then use less then
  // if > 1 then cancel till there are as many maxima
  // if < 0 then cancel till there are as many minima
  if( is_epsilon_persistent(msc,e))
    return true;

  if(t < 0)
    return (int(-t) < msc.m_num_surv_dcps[0]);

  if(t > 1)
    return (int(t) < msc.m_num_surv_dcps[2]);

  return get_persistence(msc,e) < t*r;
}

template<> inline bool is_within_treshold<false>
(const mscomplex_t & msc,int_pair_t e,fn_t t,fn_t r)
{
  if( is_epsilon_persistent(msc,e))
    return true;

  return get_persistence(msc,e) < t;
}

template<bool is_tnrm>
void mscomplex_t::simplify_impl(double f_tresh)
{
  BOOST_AUTO(cmp,bind(persistence_lt,boost::cref(*this),_2,_1));

  priority_queue<int_pair_t,int_pair_list_t,decltype(cmp)> pq(cmp);

  double f_range = *br::max_element(m_cp_fn) - *br::min_element(m_cp_fn);

  for(int i = 0 ;i < get_num_critpts();++i)
  {
    BOOST_FOREACH(int j, m_des_conn[i])
    {
      int_pair_t pr = la::make_vec(i,j);

      if(is_valid_canc_edge(*this,pr))
        pq.push(pr);
    }
  }

  while (pq.size() !=0)
  {
    int_pair_t pr = pq.top();

    pq.pop();

    order_pr_by_cp_index(*this,pr[0],pr[1]);

    if(is_valid_canc_edge(*this,pr) == false)
      continue;

    if(is_within_treshold<is_tnrm>(*this,pr,f_tresh,f_range) == false)
      break;

    cancel_pair(pr);

    m_canc_list.push_back(pr);
    m_canc_pers.push_back(get_persistence(*this,pr)/f_range);

    BOOST_FOREACH(int i, m_des_conn[pr[0]])
    BOOST_FOREACH(int j, m_asc_conn[pr[1]])
    {
      int_pair_t npr = la::make_vec(i,j);

      if(is_valid_canc_edge(*this,npr))
        pq.push(npr);
    }
  }
}

void mscomplex_t::simplify(double thresh, bool is_normalized)
{
  if(is_normalized)
    simplify_impl<true>(thresh);
  else
    simplify_impl<false>(thresh);
}

void mscomplex_t::un_simplify()
{
  for(auto it = m_canc_list.rbegin();it != m_canc_list.rend() ; ++it)
  {
    uncancel_pair(*it);
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

  m_cel_off[0] = 0;
  m_cel_off[1] = ds->m_tcc->get_num_cells_max_dim(0);
  m_cel_off[2] = ds->m_tcc->get_num_cells_max_dim(1);
  m_cel_off[3] = ds->m_tcc->get_num_cells_max_dim(2);
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
inline void bin_read_vec(std::istream &is, std::vector<T> &v,int n)
{v.resize(n);is.read((char*)(void*)v.data(),n*sizeof(T));}


template<>
inline void bin_write_vec(std::ostream &os, int_pair_list_t &v)
{
  int_list_t vlist;

  BOOST_FOREACH(int_pair_t p,v)
  {
    vlist.push_back(p[0]);
    vlist.push_back(p[1]);
  }

  os.write((const char*)(const void*)vlist.data(),vlist.size()*sizeof(int));
}

template<>
inline void bin_read_vec(std::istream &is, int_pair_list_t &v,int n)
{
  int_list_t vlist;

  bin_read_vec(is,vlist,n*2);

  v.resize(n);

  for(int i = 0 ; i < n;++i)
  {
    v[i][0] = vlist[2*i+0];
    v[i][1] = vlist[2*i+1];
  }
}

template<typename T>
inline void bin_write(std::ostream &os, const T &v)
{os.write((const char*)(const void*)&v,sizeof(T));}

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
  bin_write(os,m_num_surv_dcps[0]);
  bin_write(os,m_num_surv_dcps[1]);
  bin_write(os,m_num_surv_dcps[2]);

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
  bin_write_vec(os,m_canc_pers);

  bin_write(os,m_cel_off[0]);
  bin_write(os,m_cel_off[1]);
  bin_write(os,m_cel_off[2]);
  bin_write(os,m_cel_off[3]);

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
  bin_read(is,m_num_surv_dcps[0]);
  bin_read(is,m_num_surv_dcps[1]);
  bin_read(is,m_num_surv_dcps[2]);

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
  bin_read_vec(is,m_canc_pers,NC);

  bin_read(is,m_cel_off[0]);
  bin_read(is,m_cel_off[1]);
  bin_read(is,m_cel_off[2]);
  bin_read(is,m_cel_off[3]);

  bin_read_vec(is,nmfolds,2*N);
  m_des_mfolds.resize(N);
  m_asc_mfolds.resize(N);

  for(int i = 0 ; i < N; ++i)
  {
    bin_read_vec(is,m_des_mfolds[i],nmfolds[2*i]);
    bin_read_vec(is,m_asc_mfolds[i],nmfolds[2*i+1]);

//    if(index(i) == 2)
//      br::for_each(m_des_mfolds[i],[&](cellid_t &a){a -=m_cel_off[2];});
  }

}

void mscomplex_t::save_ascii(const std::string &f)
{
  fstream os(f.c_str(),ios::out);

  os<<"# Num Cps"<<std::endl;

  os<<get_num_critpts()<<std::endl;

  os<<"#Cp info: SL.No  cpIdx pair_idx vertNo fn "<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<i<<" ";
    os<<(int)index(i)<<" ";
    os<<(int)is_paired(i)<<" ";
    os<<(int)vertid(i)<<" ";
    os<<(fn_t)fn(i)<<" ";

    os<<std::endl;
  }

  os<<"#Cp connections: slno numDes numAsc connList"<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<(int)i<<" ";
    
    os<<(int)m_conn[0][i].size()<<" ";
    os<<(int)m_conn[1][i].size()<<" ";

    br::copy(m_conn[0][i],ostream_iterator<int>(os," "));
    br::copy(m_conn[1][i],ostream_iterator<int>(os," "));

    os<<std::endl;
  }

  os<<"#Cp geometry: slno numDes numAsc mfoldCellids"<<std::endl;

  for(uint i = 0 ; i < get_num_critpts();++i)
  {
    os<<(int)i<<" ";
    
    os<<(int)m_mfolds[0][i].size()<<" ";
    os<<(int)m_mfolds[1][i].size()<<" ";

    br::copy(m_mfolds[0][i],ostream_iterator<int>(os," "));
    br::copy(m_mfolds[1][i],ostream_iterator<int>(os," "));

    os<<std::endl;
  }

  os<<"#Cancellation sequence:slno p q pIndex qIndex pers(nrm to [0,1])"<<std::endl;

  for(uint i = 0 ; i < m_canc_list.size();++i)
  {
    os<<(int)i<<"\t";

    os<<(int)m_canc_list[i][0]<<"\t";
    os<<(int)m_canc_list[i][1]<<"\t";

    os<<(int)index(m_canc_list[i][0])<<"\t";
    os<<(int)index(m_canc_list[i][1])<<"\t";

    os<<(fn_t)m_canc_pers[i]<<"\t";

    os<<std::endl;
  }
}

}

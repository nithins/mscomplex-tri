#include <cmath>
#include <queue>
#include <limits>

#include <boost/typeof/typeof.hpp>

#include <trimesh_mscomplex.h>


using namespace std;

namespace trimesh
{
inline std::string edge_to_string(mscomplex_t *msc,int_pair_t e)
{
  std::stringstream ss;

  ss<<utls::to_string(msc->cellid(e[0]))<<"----"<<utls::to_string(msc->cellid(e[0]));

  return ss.str();
}

mscomplex_t::mscomplex_t()
  :m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]){}

mscomplex_t::~mscomplex_t(){clear();}

void mscomplex_t::set_critpt(int i, cellid_t c, char idx, cell_fn_t f, cellid_t v, bool b)
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
  m_cp_is_cancelled.resize(i,false);
  m_cp_is_boundry.resize(i,false);
  m_cp_fn.resize(i);
  m_des_conn.resize(i);
  m_asc_conn.resize(i);
}

void mscomplex_t::clear()
{
  m_cp_cellid.clear();
  m_cp_vertid.clear();
  m_cp_pair_idx.clear();
  m_cp_index.clear();
  m_cp_is_cancelled.clear();
  m_cp_is_boundry.clear();
  m_cp_fn.clear();
  m_des_conn.clear();
  m_asc_conn.clear();
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

void mscomplex_t::pair_cps(int p, int q)
{
  m_cp_pair_idx[p] = q;
  m_cp_pair_idx[q] = p;
}

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q);
  ASSERT(pair_idx(q) == p);
  ASSERT(is_canceled(p) == false);
  ASSERT(is_canceled(q) == false);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);

  conn_iter_t i,j;

  m_des_conn[p].erase(q);
  m_asc_conn[q].erase(p);

  // cps in lower of u except l
  for(i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
    for(j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
    {
      ASSERT(is_canceled(*i) == false);
      ASSERT(is_canceled(*j) == false);

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

  m_cp_is_cancelled[p] = true;
  m_cp_is_cancelled[q] = true;

  m_asc_conn[p].clear();
  m_des_conn[q].clear();
}

void mscomplex_t::uncancel_pair(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(is_canceled(p) == true && is_canceled(q) == true);
  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q && pair_idx(q) == p);

  m_cp_is_cancelled[p] = false;
  m_cp_is_cancelled[q] = false;

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

      ASSERT(is_canceled(*i) ==false && is_canceled(r) ==false);
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

  if(msc.is_canceled(e[0])||msc.is_canceled(e[1]))
    return false;

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

  int pd = msc.m_conn[GDIR_DES][pr[0]].size();
  int pa = msc.m_conn[GDIR_ASC][pr[0]].size();

  int qd = msc.m_conn[GDIR_DES][pr[1]].size();
  int qa = msc.m_conn[GDIR_ASC][pr[1]].size();

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

inline cell_fn_t get_persistence(const mscomplex_t & msc,int_pair_t e)
{
  return std::abs(msc.fn(e[0]) - msc.fn(e[1]));
}

inline bool persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
{/*
  bool eps_p0 = is_epsilon_persistent(msc,p0);
  bool eps_p1 = is_epsilon_persistent(msc,p1);

  if( eps_p0 != eps_p1)
    return eps_p0;

  cell_fn_t d0 = get_persistence(msc,p0);
  cell_fn_t d1 = get_persistence(msc,p1);

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

  cell_fn_t f00 = msc.fn(p0[0]);
  cell_fn_t f01 = msc.fn(p0[1]);
  cell_fn_t f10 = msc.fn(p1[0]);
  cell_fn_t f11 = msc.fn(p1[1]);

  cell_fn_t d1 = std::abs(f01-f00);
  cell_fn_t d2 = std::abs(f11-f10);

  if(d1 != d2)
    return d1 < d2;

  if(c00 != c10)
    return c00 < c10;

  return c01 < c11;
}

inline bool is_within_treshold(const mscomplex_t & msc,int_pair_t e,cell_fn_t t)
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
      int_pair_t pr(i,*j);

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
        int_pair_t npr(*i,*j);

        if(is_valid_canc_edge(*this,npr) && is_within_treshold(*this,npr,f_tresh))
          pq.push(npr);
      }
  }
}

void mscomplex_t::un_simplify()
{
  typedef int_pair_list_t::const_reverse_iterator revit_t;

  for(revit_t it = m_canc_list.rbegin();it != m_canc_list.rend() ; ++it)
    uncancel_pair((*it)[0],(*it)[1]);

  m_canc_list.clear();
}

void mscomplex_t::invert_for_collection()
{
  for(int i = 0 ; i < get_num_critpts(); ++i)
  {
    if(is_paired(i)== true)
      continue;

    m_des_conn[i].clear();
    m_asc_conn[i].clear();
    continue;
  }

  for(int i = 0 ; i < get_num_critpts(); ++i)
  {
    if(is_paired(i) == false)
      continue;

    ASSERT(is_paired(i) && is_paired(pair_idx(i))== true);
    ASSERT(abs(index(i)- index(pair_idx(i))) == 1);
    ASSERT(pair_idx(pair_idx(i))  == i);

    int dir = (index(i) > index(pair_idx(i)))?(0):(1);

    for(conn_iter_t j  = m_conn[dir][i].begin(); j != m_conn[dir][i].end(); ++j)
    {
      ASSERT(is_paired(*j) == false);

      m_conn[dir^1][*j].insert(i);
    }

    m_conn[dir][i].clear();
  }
}

template<typename T>
inline void bin_write_vec(std::ostream &os, std::vector<T> &v,bool purge = true)
{os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));if(purge) v.clear();}

template<typename T>
inline void bin_write(std::ostream &os, const T &v)
{os.write((const char*)(const void*)&v,sizeof(T));}

template<typename T>
inline void bin_read_vec(std::istream &is, std::vector<T> &v,int n)
{v.resize(n);is.read((char*)(void*)v.data(),n*sizeof(T));}

template<typename T>
inline void bin_read(std::istream &is, const T &v)
{is.read((char*)(void*)&v,sizeof(T));}

void mscomplex_t::stow(std::ostream &os, bool purge_data)
{
  int N = get_num_critpts();

  bin_write(os,N);

  bin_write_vec(os,m_cp_cellid,purge_data);
  bin_write_vec(os,m_cp_vertid,purge_data);
  bin_write_vec(os,m_cp_pair_idx,purge_data);
  bin_write_vec(os,m_cp_index,purge_data);
  bin_write_vec(os,m_cp_is_cancelled,purge_data);
  bin_write_vec(os,m_cp_is_boundry,purge_data);
  bin_write_vec(os,m_cp_fn,purge_data);

  int_list_t nconn(2*N);
  int_list_t adj;

  for(int i = 0 ; i < N; ++i)
  {
    nconn[2*i]   = m_des_conn[i].size();
    nconn[2*i+1] = m_asc_conn[i].size();

    std::copy(m_des_conn[i].begin(),m_des_conn[i].end(),back_inserter(adj));
    std::copy(m_asc_conn[i].begin(),m_asc_conn[i].end(),back_inserter(adj));
  }

  bin_write(os,(int)adj.size());
  bin_write_vec(os,nconn);
  bin_write_vec(os,adj);

  if(purge_data)
  {
    m_conn[0].clear();
    m_conn[1].clear();
  }

  bin_write(os,int(m_canc_list.size()));
  bin_write_vec(os,m_canc_list,purge_data);
}

void mscomplex_t::load(std::istream &is)
{
  clear();

  int N,NC;
  int_list_t nconn,adj;

  bin_read(is,N);

  bin_read_vec(is,m_cp_cellid,N);
  bin_read_vec(is,m_cp_vertid,N);
  bin_read_vec(is,m_cp_pair_idx,N);
  bin_read_vec(is,m_cp_index,N);
  bin_read_vec(is,m_cp_is_cancelled,N);
  bin_read_vec(is,m_cp_is_boundry,N);
  bin_read_vec(is,m_cp_fn,N);

  bin_read(is,NC);
  bin_read_vec(is,nconn,2*N);
  bin_read_vec(is,adj,NC);

  int_list_t::iterator a,b,c = adj.begin();

  m_des_conn.resize(N);
  m_asc_conn.resize(N);

  for(int i = 0 ; i < N; ++i)
  {
    a = c;
    b = a + (nconn[2*i]);
    c = b + (nconn[2*i+1]);

    m_des_conn[i].insert(a,b);
    m_asc_conn[i].insert(b,c);
  }

  bin_read(is,NC);
  bin_read_vec(is,m_canc_list,NC);
}


//  void mscomplex_t::save(std::ostream & os)
//  {
//    os<<"# Num Cps"<<std::endl;

//    os<<m_cps.size()<<std::endl;

//    os<<"# SL.No  cpIdx isPaired pair_idx vertNo fn "<<std::endl;

//    for(uint i = 0 ; i < m_cps.size();++i)
//    {
//      critpt_t * cp = m_cps[i];

//      os<<i<<" ";
//      os<<(int)cp->index<<" ";
//      os<<(bool)cp->is_paired<<" ";
//      os<<(int)cp->pair_idx<<" ";
//      os<<cp->vert_idx<<" ";
//      os<<(cell_fn_t)cp->fn<<" ";

//      os<<std::endl;
//    }

//    os<<"#slno numDes numAsc connList"<<std::endl;

//    for(uint i = 0 ; i < m_cps.size();++i)
//    {
//      critpt_t * cp = m_cps[i];

//      os<<(int)i<<" ";

//      if(cp->is_paired == true)
//      {
//        os<<"0 0"<<std::endl;
//        continue;
//      }
//      os<<(int)cp->conn[0].size()<<" ";
//      os<<(int)cp->conn[1].size()<<" ";

//      for(uint dir = 0 ; dir <2 ;++dir)
//      {
//        conn_set_t &conn = cp->conn[dir];

//        for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
//        {
//          os<<*it<<" ";
//        }
//      }
//      os<<std::endl;
//    }
//  }

}



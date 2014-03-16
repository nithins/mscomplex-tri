#include <cmath>
#include <queue>
#include <limits>

#include <boost/foreach.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/set.hpp>

#include <trimesh_mscomplex.h>
#include <trimesh_dataset.h>

using namespace std;
namespace br    = boost::range;
namespace badpt = boost::adaptors;

namespace trimesh{

/*****************************************************************************/

/// \brief Data structure to extract MSC geometry at desired multires versions
class merge_dag_t
{
public:
  /// \brief For the given critical point (cp) , get the cps that
  /// that contribute to its geometry at the given resolution
  void get_contrib_cps(std::vector<int> &l, eGDIR dir,int cp,int res) const;

  /// \brief Build the merge dag for the given dim and dir
  template <eGDIR dir,int dim>
  void build(mscomplex_ptr_t msc);

  /// \brief Serialization
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);

//private:

  struct node_t
  {
    int base;
    int other;
    int canc_no;

    node_t():base(-1),other(-1),canc_no(-1){}
    node_t(int b,int o,int c):base(b),other(o),canc_no(c){}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version);
  };

  std::vector<node_t>  m_nodes;
  std::vector<int>     m_cp_geom[2];

  inline const node_t& get_node(int i) const;
  inline int  get_ncps() const;
};

/*---------------------------------------------------------------------------*/

inline const merge_dag_t::node_t& merge_dag_t::get_node(int i) const
{
  static node_t defnode;

  if(i < 0)
    return defnode;
  return m_nodes[i];
}

/*---------------------------------------------------------------------------*/

inline int merge_dag_t::get_ncps() const
{
  return m_cp_geom[0].size();
}

/*---------------------------------------------------------------------------*/

template <int dim,int odim>
inline bool is_dim_pair(mscomplex_ptr_t msc,int_pair_t pr)
{return (msc->index(pr.first) == dim) && (msc->index(pr.second) == odim);}

template <eGDIR dir,int dim>
void merge_dag_t::build(mscomplex_ptr_t msc)
{
  const eGDIR odir = (dir == ASC)?(DES):(ASC);
  const int   odim = (dir == ASC)?(dim +1):(dim -1);

  int version = msc->m_multires_version;

  int ncps = msc->get_num_critpts();

  m_cp_geom[dir].resize(ncps);

  br::copy(boost::counting_range(-ncps,0)|badpt::reversed,m_cp_geom[dir].begin());

  msc->set_multires_version(msc->m_canc_list.size()+1);

  BOOST_FOREACH(int_pair_t pr, msc->m_canc_list
                |badpt::sliced(version,msc->m_canc_list.size())
                |badpt::transformed(bind(order_by_dir_index<dir>,msc,_1))
                |badpt::filtered(bind(is_dim_pair<dim,odim>,msc,_1)))
  {
    int p = pr.first, q = pr.second;

    int pnode = m_cp_geom[dir][p],pcancno = msc->m_cp_cancno[p];

    ensure(get_node(pnode).canc_no < pcancno,
           "earlier pnode has formed from a later cancellation");

    if( pnode >=  ncps || msc->m_mfolds[dir][pnode].size() > 0)
    {
      BOOST_FOREACH(int r,msc->m_conn[odir][q])
      {
        int rnode = m_cp_geom[dir][r];
        m_cp_geom[dir][r] = m_nodes.size();
        m_nodes.push_back(node_t(rnode,pnode,pcancno));
        ensure(get_node(rnode).canc_no < pcancno,
               "earlier rnode has formed from a later cancellation");
      }
    }
  }
  msc->set_multires_version(version);
}

/*---------------------------------------------------------------------------*/

void merge_dag_t::get_contrib_cps
(std::vector<int> &l, eGDIR dir,int cp,int res) const
{
  int g = m_cp_geom[dir][cp];

  while(get_node(g).canc_no >= res )
    g = get_node(g).base;

  std::set<int>    visited;

  std::stack<int> stk;

  stk.push(g);

  while(stk.size() != 0 )
  {
    int g = stk.top();
    stk.pop();

    if(visited.count(g) != 0 )
      continue;

    visited.insert(g);

    if( g < 0)
      l.push_back(-g-1);

    node_t gnode = get_node(g);

    if(gnode.base != -1)
    {
      stk.push(gnode.base);
      stk.push(gnode.other);
    }
  }
}

/*---------------------------------------------------------------------------*/

template<class Archive>
void merge_dag_t::node_t::serialize(Archive & ar, const unsigned int version)
{
  ar& BOOST_SERIALIZATION_NVP(base);
  ar& BOOST_SERIALIZATION_NVP(other);
  ar& BOOST_SERIALIZATION_NVP(canc_no);
}

/*---------------------------------------------------------------------------*/

template<class Archive>
void merge_dag_t::serialize(Archive & ar, const unsigned int version)
{
  ar& BOOST_SERIALIZATION_NVP(m_nodes);
  ar& BOOST_SERIALIZATION_NVP(m_cp_geom);
}

/*****************************************************************************/




/*****************************************************************************/

mscomplex_t::mscomplex_t()
  :m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),
    m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
    m_multires_version(0),
    m_fmax(std::numeric_limits<fn_t>::min()),
    m_fmin(std::numeric_limits<fn_t>::max()),
    m_merge_dag(new merge_dag_t){}

/*---------------------------------------------------------------------------*/

mscomplex_t::mscomplex_t(std::string fname)
  :m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]),
    m_des_mfolds(m_mfolds[0]),m_asc_mfolds(m_mfolds[1]),
    m_multires_version(0),
    m_fmax(std::numeric_limits<fn_t>::min()),
    m_fmin(std::numeric_limits<fn_t>::max()),
    m_merge_dag(new merge_dag_t)
{
  load(fname);
}

/*---------------------------------------------------------------------------*/

mscomplex_t::~mscomplex_t(){clear();}

/*---------------------------------------------------------------------------*/

void mscomplex_t::set_critpt
(int i, cellid_t c, char idx, fn_t f, cellid_t v, bool b)
{
  m_cp_cellid[i]     = c;
  m_cp_vertid[i]     = v;
  m_cp_index[i]      = idx;
  m_cp_fn[i]         = f;
  m_cp_is_boundry[i] = b;

  m_fmax = std::max(f,m_fmax);
  m_fmin = std::min(f,m_fmin);
}

/*---------------------------------------------------------------------------*/

void  mscomplex_t::resize(int i)
{
  m_cp_cellid.resize(i,invalid_cellid);
  m_cp_vertid.resize(i,invalid_cellid);
  m_cp_index.resize(i,-1);
  m_cp_pair_idx.resize(i,-1);
  m_cp_cancno.resize(i,i);
  m_cp_is_boundry.resize(i,false);
  m_cp_fn.resize(i);
  m_des_conn.resize(i);
  m_asc_conn.resize(i);
  m_des_mfolds.resize(i);
  m_asc_mfolds.resize(i);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::clear()
{
  m_cp_cellid.clear();
  m_cp_vertid.clear();
  m_cp_pair_idx.clear();
  m_cp_index.clear();
  m_cp_cancno.clear();
  m_cp_is_boundry.clear();
  m_cp_fn.clear();
  m_des_conn.clear();
  m_asc_conn.clear();
  m_des_mfolds.clear();
  m_asc_mfolds.clear();

  m_fmax = std::numeric_limits<fn_t>::min();
  m_fmin = std::numeric_limits<fn_t>::max();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::connect_cps(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(is_not_paired(p) && is_not_paired(q));
  ASSERT(m_des_conn[p].count(q) == m_asc_conn[q].count(p));

  m_des_conn[p].insert(q);
  m_asc_conn[q].insert(p);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::disconnect_cps(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(is_not_paired(p) && is_not_paired(q));
  ASSERT(m_des_conn[p].count(q) == m_asc_conn[q].count(p));

  m_des_conn[p].erase(m_des_conn[p].find(q));
  m_asc_conn[q].erase(m_asc_conn[q].find(p));
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ensure(m_multires_version == m_canc_list.size(),
         "Cannot cancel pair !! Ms complex resolution is not coarsest.");
  ensure(index(p) == index(q)+1,
         "indices do not differ by 1");
  ensure(m_cp_pair_idx[p] == -1 && m_cp_pair_idx[q] == -1,
         "p/q has already been paired");
  ensure(m_des_conn[p].count(q)  == m_asc_conn[q].count(p),
         "p is not connected to q");
  ensure(m_des_conn[p].count(q) == 1,
         "p and q are multiply connected");

  m_cp_cancno[p] = m_canc_list.size();
  m_cp_cancno[q] = m_canc_list.size();
  m_canc_list.push_back(make_pair(p,q));

  cancel_pair();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::cancel_pair()
{
  try
  {
    ensure(is_in_range(m_multires_version,0,m_canc_list.size()),
           "invalid cancellation position");

    int p = m_canc_list[m_multires_version].first;
    int q = m_canc_list[m_multires_version].second;

    m_multires_version++;

    ASSERT(index(p) == index(q)+1);
    ASSERT(m_cp_pair_idx[p] == -1);
    ASSERT(m_cp_pair_idx[q] == -1);
    ASSERT(m_des_conn[p].count(q) == 1);
    ASSERT(m_asc_conn[q].count(p) == 1);
//    ASSERT(m_des_conn[p][q] == 1);
//    ASSERT(m_asc_conn[q][p] == 1);

    m_cp_pair_idx[p] = q;
    m_cp_pair_idx[q] = p;

    m_des_conn[p].erase(q);
    m_asc_conn[q].erase(p);

    // cps in lower of u except l
    BOOST_FOREACH(int u,m_des_conn[p])
    BOOST_FOREACH(int v,m_asc_conn[q])
    {
      ASSERT(is_paired(u) == false);
      ASSERT(is_paired(v) == false);

      connect_cps(u,v);
    }

    BOOST_FOREACH(int pr,m_des_conn[p]) m_asc_conn[pr].erase(p);
    BOOST_FOREACH(int pr,m_asc_conn[p]) m_des_conn[pr].erase(p);
    BOOST_FOREACH(int pr,m_des_conn[q]) m_asc_conn[pr].erase(q);
    BOOST_FOREACH(int pr,m_asc_conn[q]) m_des_conn[pr].erase(q);
  }
  catch (std::exception &ex)
  {
//    ex.push(_FFL);
    throw;
  }
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::anticancel_pair()
{
  m_multires_version--;

  ensure(is_in_range(m_multires_version,0,m_canc_list.size()),
         "invalid cancellation position");

  int p = m_canc_list[m_multires_version].first;
  int q = m_canc_list[m_multires_version].second;

  ASSERT(index(p) == index(q)+1);
  ASSERT(m_cp_pair_idx[p] == q);
  ASSERT(m_cp_pair_idx[q] == p);

  BOOST_FOREACH(int pr,m_des_conn[p]) m_asc_conn[pr].insert(p);
  BOOST_FOREACH(int pr,m_asc_conn[p]) m_des_conn[pr].insert(p);
  BOOST_FOREACH(int pr,m_des_conn[q]) m_asc_conn[pr].insert(q);
  BOOST_FOREACH(int pr,m_asc_conn[q]) m_des_conn[pr].insert(q);

  // cps in lower of u except l
  BOOST_FOREACH(int u,m_des_conn[p])
      BOOST_FOREACH(int v,m_asc_conn[q])
  {
    disconnect_cps(u,v);
  }

  //    ASSERT(m_des_conn[p][q] == 1);
  //    ASSERT(m_asc_conn[q][p] == 1);

  m_cp_pair_idx[p] = -1;
  m_cp_pair_idx[q] = -1;

  connect_cps(p,q);

  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);
}

/*---------------------------------------------------------------------------*/

inline bool is_epsilon_persistent(const mscomplex_t &msc,int_pair_t e )
{return (msc.vertid(e.first) == msc.vertid(e.second));}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

inline int get_num_new_edges(const mscomplex_t &msc, int_pair_t pr)
{
  order_pr_by_cp_index(msc,pr.first,pr.second);

  int pd=msc.m_des_conn[pr.first].size();
  int pa=msc.m_asc_conn[pr.first].size();
  int qd=msc.m_des_conn[pr.second].size();
  int qa=msc.m_asc_conn[pr.second].size();

  return (pd - 1)*(qa - 1) - (pd + qd + pa + qa -1);
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

inline fn_t get_persistence(const mscomplex_t & msc,int_pair_t e)
{return std::abs(msc.fn(e.first) - msc.fn(e.second));}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

inline bool persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
{
  order_pr_by_cp_index(msc,p0.first,p0.second);
  order_pr_by_cp_index(msc,p1.first,p1.second);

  cellid_t v00 = msc.vertid(p0.first);
  cellid_t v01 = msc.vertid(p0.second);
  cellid_t v10 = msc.vertid(p1.first);
  cellid_t v11 = msc.vertid(p1.second);

  cellid_t c00 = msc.cellid(p0.first);
  cellid_t c01 = msc.cellid(p0.second);
  cellid_t c10 = msc.cellid(p1.first);
  cellid_t c11 = msc.cellid(p1.second);

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

  fn_t f00 = msc.fn(p0.first);
  fn_t f01 = msc.fn(p0.second);
  fn_t f10 = msc.fn(p1.first);
  fn_t f11 = msc.fn(p1.second);

  fn_t d1 = std::abs(f01-f00);
  fn_t d2 = std::abs(f11-f10);

  if(d1 != d2)
    return d1 < d2;

  if(c00 != c10)
    return c00 < c10;

  return c01 < c11;
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

inline bool is_within_treshold(const mscomplex_t & msc,int_pair_t e,fn_t t)
{return (is_epsilon_persistent(msc,e) || get_persistence(msc,e) < t);}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

void mscomplex_t::simplify(double f_tresh, bool is_nrm, int req_nmin, int req_nmax)
{
  BOOST_AUTO(cmp,bind(persistence_lt,boost::cref(*this),_2,_1));

  priority_queue<int_pair_t,int_pair_list_t,BOOST_TYPEOF(cmp)> pq(cmp);

  double f_rng = m_fmax - m_fmin;

  if (is_nrm) f_tresh *= f_rng;

  ensure(f_tresh >=0 && f_tresh <= f_rng,"tresh out of range");

  for(int i = 0 ;i < get_num_critpts();++i)
  {
    BOOST_FOREACH(int j,m_des_conn[i])
    {
      int_pair_t pr(i,j);

      if(is_valid_canc_edge(*this,pr) && is_within_treshold(*this,pr,f_tresh))
        pq.push(pr);
    }
  }

  int nmin = boost::distance
      (cp_range()|badpt::filtered(bind(&mscomplex_t::is_not_paired,this,_1))
       |badpt::filtered(bind(&mscomplex_t::is_minima,this,_1)));

  int nmax = boost::distance
      (cp_range()|badpt::filtered(bind(&mscomplex_t::is_not_paired,this,_1))
       |badpt::filtered(bind(&mscomplex_t::is_maxima,this,_1)));

  while (pq.size() !=0 && req_nmin <= nmin && req_nmax <= nmax )
  {
    int_pair_t pr = pq.top();

    pq.pop();

    if(is_valid_canc_edge(*this,pr) == false)
      continue;

    int p = pr.first,q = pr.second;

    order_pr_by_cp_index(*this,p,q);

    cancel_pair(p,q);

    BOOST_FOREACH(int i,m_des_conn[p])
    BOOST_FOREACH(int j,m_asc_conn[q])
    {
      int_pair_t npr(i,j);

      if(is_valid_canc_edge(*this,npr) && is_within_treshold(*this,npr,f_tresh))
        pq.push(npr);
    }

    if(index(p) == 2) --nmax; else --nmin;
  }
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::set_multires_version(int version)
{
  for(int i = m_multires_version ; i > version && i>0; i--)
  {
    anticancel_pair();
  }

  for( int i = m_multires_version; i <version&&i<m_canc_list.size(); ++i)
  {
    cancel_pair();
  }
}

/*---------------------------------------------------------------------------*/

bool is_t_lt_abs_diff(const mscomplex_t * msc, double t,int i)
{
  ASSERT(is_in_range(i,0,msc->m_canc_list.size()+1));

  if( i == msc->m_canc_list.size())
    return false;
  else
    return t < std::abs
        (msc->fn(msc->m_canc_list[i].first) -
         msc->fn(msc->m_canc_list[i].second));
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

int mscomplex_t::get_multires_version_for_thresh(double t, bool is_nrm) const
{
  if(is_nrm)
    t *= (m_fmax - m_fmin);

  return *std::upper_bound
      (boost::counting_iterator<int>(0),
       boost::counting_iterator<int>(m_canc_list.size()+1),
       t,bind(is_t_lt_abs_diff,this,_1,_2));
}

/*---------------------------------------------------------------------------*/

typedef std::map<int,int_list_t> contrib_t;

template<eGDIR dir,int dim>
void get_contrib_cps(mscomplex_ptr_t msc,contrib_t& contrib)
{
  const eGDIR odir = (dir == ASC)?(DES):(ASC);
  const int   odim = (dir == ASC)?(dim +1):(dim -1);

  contrib_t contrib_map;

  // first, obtain a sequence of cancellations so that
  // a) each pair is ordered by dir ..
  //    i.e if dir is DES then look at a (2,1) as (2,1) and not (1,2)
  // b) pairs that have not yet been cancelled are removed
  // c) pairs whose first element have index == dim

  BOOST_AUTO(canc_rng,msc->m_canc_list
             |badpt::sliced(0,msc->m_multires_version)
             |badpt::transformed(bind(order_by_dir_index<dir>,msc,_1))
             |badpt::filtered(bind(is_dim_pair<dim,odim>,msc,_1))
             );

  // This part computes for each cancelled critical point,
  // the surviving critical points to which it contributes its
  // finest resolution geometry .
  BOOST_FOREACH(int_pair_t pr,canc_rng|badpt::reversed)
  {
    int p = pr.first,q = pr.second;

    int_list_t & pcontrib = contrib_map[p];

    // for each qa in the asc conn of q:
    BOOST_FOREACH(int qa, msc->m_conn[odir][q])
    {
      // a) if qa is not paired ..
      if(msc->is_not_paired(qa))
      {
        // .. p contributes to qa.
        pcontrib.push_back(qa);
      }
      // b) if qa is paired and qa's pair and q have same index ..
      else if(msc->index(q) == msc->index(msc->pair_idx(qa)))
      {
        // .. then foreach qaqa that qa contributes to  ..
        BOOST_FOREACH(int qaqa,contrib_map[qa])
        {
          // .. p contributes to qaqa.
          pcontrib.push_back(qaqa);

          // pdpd has to be a surviving cp
          ASSERT(msc->is_not_paired(qaqa));
        }
      }
    }
  }

  // This part just reverses info computed earlier.
  // i.e each surviving critical point will have a
  // list of cancelled critical points that contribute
  // their finest res geometry to it.

  BOOST_FOREACH(int_pair_t pr, canc_rng)
  {
    int p = pr.first;

    BOOST_FOREACH(int p_, contrib_map[p])
    {
      contrib[p_].push_back(p);
    }
  }

  // Sanity checks.
  BOOST_FOREACH(contrib_t::value_type pr, contrib)
  {
    ensure(msc->index(pr.first) == dim,
           "Computed the contribution of index != dim");
    ensure(msc->is_not_paired(pr.first),
           "Computed the contribution to a cancelled cp");

    BOOST_FOREACH(int ccp, pr.second)
    {
      ensure(msc->index(ccp) == dim,
             "other dim cp is attempting to contribute");
      ensure(msc->is_paired(ccp),
             "an unpaired cp is attempting to contribute");
      ensure(msc->index(msc->pair_idx(ccp)) == odim,
             "wrong cancellation type cp is attempting to contribute");
    }
  }

  // finally add each cp to its own list of contributors
  BOOST_FOREACH(int cp, msc->cp_range()
                |badpt::filtered(bind(&mscomplex_t::is_not_paired,msc,_1))
                |badpt::filtered(bind(&mscomplex_t::is_index_i_cp<dim>,msc,_1)))
  {
    contrib[cp].push_back(cp);
  }
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

template <eGDIR dir, int dim>
inline void __collect_mfolds(mscomplex_ptr_t msc, dataset_ptr_t ds)
{
  contrib_t contrib;
  get_contrib_cps<dir,dim>(msc,contrib);

  BOOST_FOREACH(contrib_t::value_type pr,contrib)
  {
    ensure(msc->m_mfolds[dir][pr.first].size() == 0,
        "Geom for cp has already been collected");

    ds->get_mfold<dir>(msc->m_mfolds[dir][pr.first],
        pr.second|badpt::transformed(bind(&mscomplex_t::cellid,msc,_1)));

  }

//  msc->m_merge_dag->build<dir,dim>(msc);
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

void mscomplex_t::collect_mfolds(eGDIR dir, int dim, dataset_ptr_t ds)
{
  if(dir == ASC && dim == 0 ) __collect_mfolds<ASC,0>(shared_from_this(),ds);
  if(dir == ASC && dim == 1 ) __collect_mfolds<ASC,1>(shared_from_this(),ds);
  if(dir == DES && dim == 1 ) __collect_mfolds<DES,1>(shared_from_this(),ds);
  if(dir == DES && dim == 2 ) __collect_mfolds<DES,2>(shared_from_this(),ds);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::get_mfold(eGDIR dir, int cp,cellid_list_t &mfold,int ver)
{
  if (ver == -1 )
    ver = m_multires_version;

  int_list_t l;
  m_merge_dag->get_contrib_cps(l,dir,cp,ver);

  BOOST_FOREACH(int c,l)
  {
    br::copy(m_mfolds[dir][c],back_inserter(mfold));
  }
}

/*---------------------------------------------------------------------------*/

template<class Archive>
void mscomplex_t::serialize(Archive & ar, const unsigned int version)
{
  ar& BOOST_SERIALIZATION_NVP(m_cp_cellid);
  ar& BOOST_SERIALIZATION_NVP(m_cp_vertid);
  ar& BOOST_SERIALIZATION_NVP(m_cp_pair_idx);
  ar& BOOST_SERIALIZATION_NVP(m_cp_index);
  ar& BOOST_SERIALIZATION_NVP(m_cp_cancno);
  ar& BOOST_SERIALIZATION_NVP(m_cp_is_boundry);
  ar& BOOST_SERIALIZATION_NVP(m_cp_fn);
  ar& BOOST_SERIALIZATION_NVP(m_canc_list);
  ar& BOOST_SERIALIZATION_NVP(m_des_conn);
  ar& BOOST_SERIALIZATION_NVP(m_asc_conn);
  ar& BOOST_SERIALIZATION_NVP(m_des_mfolds);
  ar& BOOST_SERIALIZATION_NVP(m_asc_mfolds);
  ar& BOOST_SERIALIZATION_NVP(m_multires_version);
  ar& BOOST_SERIALIZATION_NVP(m_fmax);
  ar& BOOST_SERIALIZATION_NVP(m_fmin);
  ar& BOOST_SERIALIZATION_NVP(*m_merge_dag);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::save_bin(std::ostream &os) const
{
  boost::archive::binary_oarchive oa(os);
  oa << BOOST_SERIALIZATION_NVP(*this);
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::load_bin(std::istream &is)
{
  boost::archive::binary_iarchive ia(is);
  ia >> BOOST_SERIALIZATION_NVP(*this);
}

/*---------------------------------------------------------------------------*/

template<>
int_pair_t order_by_dir_index<DES>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr.first) < msc->index(pr.second))std::swap(pr.first,pr.second);return pr;}

template<>
int_pair_t order_by_dir_index<ASC>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr.first) > msc->index(pr.second))std::swap(pr.first,pr.second); return pr;}


/*****************************************************************************/
}// namespace trimesh



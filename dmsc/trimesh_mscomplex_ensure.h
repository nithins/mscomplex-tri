#ifndef TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED
#define TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED

#include <trimesh_mscomplex.h>
#include <boost/bind.hpp>

namespace trimesh
{

inline void order_pr_by_cp_index(const mscomplex_t &msc,int &p,int &q)
{if(msc.index(p) < msc.index(q))std::swap(p,q);}

template<eGDIR dir>
int_pair_t order_by_dir_index(mscomplex_ptr_t msc,int_pair_t pr);

inline int  mscomplex_t::get_num_critpts() const
{return m_cp_cellid.size();}

inline const char& mscomplex_t::index(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_index.size()));
  ASSERT(is_in_range(m_cp_index[i],0,3));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_index[i];
}

inline const int& mscomplex_t::pair_idx(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_pair_idx.size()));
  ASSERT(is_in_range(m_cp_pair_idx[i],0,(int)m_cp_pair_idx.size()));
  ASSERT(i == m_cp_pair_idx[m_cp_pair_idx[i]]);}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_pair_idx[i];
}

inline bool mscomplex_t::is_paired(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_pair_idx.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return (m_cp_pair_idx[i] != -1);
}

inline bool mscomplex_t::is_not_paired(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_pair_idx.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return (m_cp_pair_idx[i] == -1);
}

inline const char& mscomplex_t::is_boundry(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_is_boundry.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_is_boundry[i];
}

inline const cellid_t& mscomplex_t::cellid(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_cellid.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_cellid[i];
}

inline const cellid_t& mscomplex_t::vertid(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_vertid.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_vertid[i];
}

inline const fn_t& mscomplex_t::fn(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_fn.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_fn[i];
}

inline bool mscomplex_t::is_extrema(int i) const
{return (index(i) == 0 || index(i) == 2);}

inline bool mscomplex_t::is_saddle(int i) const
{return (index(i) == 1);}

inline bool mscomplex_t::is_maxima(int i) const
{return (index(i) == 2);}

inline bool mscomplex_t::is_minima(int i) const
{return (index(i) == 0);}

template<> inline bool mscomplex_t::is_index_i_cp<0>(int i) const
{return (index(i) == 0);}

template<> inline bool mscomplex_t::is_index_i_cp<1>(int i) const
{return (index(i) == 1);}

template<> inline bool mscomplex_t::is_index_i_cp<2>(int i) const
{return (index(i) == 2);}




inline std::string mscomplex_t::cp_info (int cp_no) const
{
  std::stringstream ss;

  ss<<std::endl;
  ss<<"cp_no        ::"<<cp_no<<std::endl;
  ss<<"cellid       ::"<<cellid(cp_no)<<std::endl;
//    ss<<"vert cell    ::"<<vertid(cp_no)<<std::endl;
  ss<<"index        ::"<<(int)index(cp_no)<<std::endl;
//      ss<<"fn           ::"<<fn(cp_no)<<std::endl;
//    ss<<"is_cancelled ::"<<is_canceled(cp_no)<<std::endl;
//    ss<<"is_paired    ::"<<is_paired(cp_no)<<std::endl;
  ss<<"pair_idx     ::"<<pair_idx(cp_no)<<std::endl;
  return ss.str();
}

inline std::string mscomplex_t::cp_conn (int i) const
{
  std::stringstream ss;

  ss<<std::endl<<"des = ";

  for(const_conn_iter_t it = m_des_conn[i].begin(); it != m_des_conn[i].end(); ++it)
    ss<<cellid(*it);

  ss<<std::endl<<"asc = ";

  for(const_conn_iter_t it = m_asc_conn[i].begin(); it != m_asc_conn[i].end(); ++it)
    ss<<cellid(*it);

  ss<<std::endl;

  return ss.str();
}

inline bool is_valid_canc_edge(const mscomplex_t &msc,int_pair_t e )
{
  order_pr_by_cp_index(msc,e.first,e.second);

  if(msc.is_paired(e.first) || msc.is_paired(e.second))
    return false;

  if(msc.is_boundry(e.first) != msc.is_boundry(e.second))
    return false;

  ASSERT(msc.m_des_conn[e.first].count(e.second) == msc.m_asc_conn[e.second].count(e.first));

  if(msc.m_des_conn[e.first].count(e.second) != 1)
    return false;

  return true;
}

inline void mscomplex_t::save(const std::string &f) const
{
  std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);
  ensure(fs.is_open(),"file not found!!");
  save_bin(fs);
}
inline void mscomplex_t::load(const std::string &f)
{
  std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);
  ensure(fs.is_open(),"file not found!!");
  load_bin(fs);
}


}
#endif

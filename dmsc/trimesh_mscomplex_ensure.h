#ifndef TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED
#define TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED

#include <trimesh_mscomplex.h>
#include <boost/bind.hpp>

namespace trimesh
{

inline void order_pr_by_cp_index(const mscomplex_t &msc,int &p,int &q)
{if(msc.index(p) < msc.index(q))std::swap(p,q);}

inline mscomplex_t::iterator mscomplex_t::begin() const
{return iterator(0);}

inline mscomplex_t::iterator mscomplex_t::end() const
{return iterator(get_num_critpts());}

inline mscomplex_t::id_iterator mscomplex_t::id_begin() const
{return id_iterator(shared_from_this(),begin());}

inline mscomplex_t::id_iterator mscomplex_t::id_end() const
{return id_iterator(shared_from_this(),end());}

inline mscomplex_t::fiterator mscomplex_t::fbegin(mscomplex_t::filter_t f) const
{return boost::make_filter_iterator(f,begin(),end());}

inline mscomplex_t::fiterator mscomplex_t::fend(mscomplex_t::filter_t f) const
{return boost::make_filter_iterator(f,end(),end());}

inline mscomplex_t::fiterator mscomplex_t::fbegin(mscomplex_t::memb_filter_t f) const
{return fbegin(boost::bind(f,this,_1));}

inline mscomplex_t::fiterator mscomplex_t::fend(mscomplex_t::memb_filter_t f) const
{return fend(boost::bind(f,this,_1));}

inline mscomplex_t::id_fiterator mscomplex_t::id_fbegin(mscomplex_t::filter_t f) const
{return id_fiterator(shared_from_this(),fbegin(f));}

inline mscomplex_t::id_fiterator mscomplex_t::id_fend(mscomplex_t::filter_t f) const
{return id_fiterator(shared_from_this(),fend(f));}

inline int  mscomplex_t::get_num_critpts() const
{return m_cp_cellid.size();}

inline char& mscomplex_t::index(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_index.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_index[i];
}

inline const char& mscomplex_t::index(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_index.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_index[i];
}

inline int& mscomplex_t::pair_idx(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_pair_idx.size()));
  /*ASSERT(is_in_range(m_cp_pair_idx[i],0,(int)m_cp_pair_idx.size()));
  ASSERT(i == m_cp_pair_idx[m_cp_pair_idx[i]]);*/}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_pair_idx[i];
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

inline char& mscomplex_t::is_canceled(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_is_cancelled.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_is_cancelled[i];
}

inline const char& mscomplex_t::is_canceled(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_is_cancelled.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_is_cancelled[i];
}

inline char& mscomplex_t::is_boundry(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_is_boundry.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_is_boundry[i];
}

inline const char& mscomplex_t::is_boundry(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_is_boundry.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_is_boundry[i];
}

inline cellid_t& mscomplex_t::_lv_cellid(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_cellid.size()));}
  catch(assertion_error e)
  {e.push(_FFL).push(SVAR(i)).push(SVAR(m_cp_cellid.size()));throw;}

  return m_cp_cellid[i];
}

inline const cellid_t& mscomplex_t::_rv_cellid(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_cellid.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_cellid[i];
}

inline const cellid_t& mscomplex_t::cellid(int i) const{return _rv_cellid(i);}
inline cellid_t& mscomplex_t::cellid(int i) {return _lv_cellid(i);}

inline cellid_t& mscomplex_t::vertid(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_vertid.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_vertid[i];
}

inline const cellid_t& mscomplex_t::vertid(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_vertid.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_vertid[i];
}

inline cell_fn_t& mscomplex_t::fn(int i)
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_fn.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_fn[i];
}

inline const cell_fn_t& mscomplex_t::fn(int i) const
{
  try{ASSERT(is_in_range(i,0,(int)m_cp_fn.size()));}
  catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

  return m_cp_fn[i];
}

inline bool mscomplex_t::is_extrema(int i) const
{return (index(i) == 0 || index(i) == 2);}

inline bool mscomplex_t::is_saddle(int i) const
{return (index(i) == 1);}

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
}
#endif

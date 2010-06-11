#ifndef TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED
#define TRIMESH_MSCOMPLEX_ENSURE_H_INCLUDED

#include <trimesh_mscomplex.h>

// bunch of predicates that throw when I suspect something could be logically
// wrong with the state of the MS complex.. to be disabled in release builds

//#ifndef NDEBUG
#define USE_ENSURE_PREDICATES
//#endif

namespace trimesh
{

  inline void ensure_index_one_separation(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    if((msc->m_cps[e[0]]->index +1 != msc->m_cps[e[1]]->index)&&
       (msc->m_cps[e[1]]->index +1 != msc->m_cps[e[0]]->index))
      throw std::logic_error("index one separation violated");
#endif
  }

  inline void ensure_ordered_index_one_separation(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    if(msc->m_cps[e[1]]->index +1 != msc->m_cps[e[0]]->index)
      throw std::logic_error("ordered index one separation violated");
#endif
  }

  inline void ensure_not_cancelled(mscomplex_t *msc,uint c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(msc->m_cps[c]->isCancelled == true)
        throw std::logic_error("cancellation state violated");
#endif

  }

  inline std::string dir_to_string(uint dir)
  {
    return (dir == 0)?("des"):("asc");
  }

  inline std::string idx_to_string(mscomplex_t *msc,uint i)
  {
    std::stringstream ss;

    ss<<msc->m_cps[i]->cellid;

    return ss.str();
  }

  template <typename T>
  inline void log_line(const T& t)
  {
    std::cout<<t<<"\n";
  }

  template <typename T1,typename T2>
  inline void log_line(const T1& t1,const T2& t2)
  {
    std::cout<<t1<<" "<<t2<<"\n";
  }

  inline std::string edge_to_string(mscomplex_t *msc,uint_pair_t e)
  {
    std::stringstream ss;

    ss<<idx_to_string(msc,e[0])<<"----"<<idx_to_string(msc,e[1]);

    return ss.str();
  }

  inline void ensure_max_two_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES

    order_pr_by_cp_index(msc,e);

    if(msc->m_cps[e[1]]->index +1 == msc->m_cps[e[0]]->index)
    {
      if(msc->m_cps[e[0]]->conn[0].count(e[1]) > 2 ||
         msc->m_cps[e[1]]->conn[1].count(e[0]) > 2)
      {

        std::stringstream ss;

        ss<<"failed to ensure max two connectivity\n"<<
            edge_to_string(msc,e)<<"\n";

        throw std::logic_error(ss.str());
      }
    }


#endif

  }


  inline bool is_multiple_edge(mscomplex_t *msc,uint_pair_t e)
  {
    order_pr_by_cp_index(msc,e);

    if(msc->m_cps[e[1]]->index +1 == msc->m_cps[e[0]]->index)
    {
      if(msc->m_cps[e[0]]->conn[0].count(e[1]) > 1 ||
         msc->m_cps[e[1]]->conn[1].count(e[0]) > 1)
        return true;
    }

    return false;
  }

  inline bool is_saddle(mscomplex_t *msc ,uint c)
  {
    if(msc->m_cps[c]->index == 1 || msc->m_cps[c]->index == 2)
      return true;

    return false;

  }

  inline void ensure_ordered_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    ensure_ordered_index_one_separation(msc,e);

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      if(msc->m_cps[e[dir]]->conn[dir].count(e[dir^1]) == 0)
      {

        std::stringstream ss;
        ss<<"connectivity voilated \n";
        ss<<dir_to_string(dir);
        ss<<msc->m_cps[e[dir]]->cellid<<" does not contain "<<
            msc->m_cps[e[dir^1]]->cellid<<"\n";

        throw std::logic_error(ss.str(  ));
      }
#endif

  }

  inline void ensure_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    order_pr_by_cp_index(msc,e);

    ensure_ordered_connectivity(msc,e);
#endif
  }


  inline void ensure_single_connectivity(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    ensure_ordered_index_one_separation(msc,e);

    for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      if(msc->m_cps[e[dir]]->conn[dir].count(e[dir^1]) != 1)
        throw std::logic_error("connectivity violated");
#endif
  }

  inline void ensure_cellid_critical(mscomplex_t * msc,cellid_t c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(msc->m_id_cp_map.count(c) == 0)
      throw std::logic_error("cellid not entered as critical in msc");
#endif
  }

  inline void ensure_cp_is_cancelled(mscomplex_t *msc,uint i)
  {
#ifdef USE_ENSURE_PREDICATES
    if(!msc->m_cps[i]->isCancelled)
      throw std::logic_error("failed to ensure cp is not canceled");
#endif
  }

  inline void ensure_cp_is_not_cancelled(mscomplex_t *msc,uint i)
  {
#ifdef USE_ENSURE_PREDICATES
    if(msc->m_cps[i]->isCancelled)
      throw std::logic_error("failed to ensure cp is canceled");
#endif
  }

  inline void ensure_pairing(mscomplex_t *msc,uint_pair_t e)
  {
#ifdef USE_ENSURE_PREDICATES
    if(!msc->m_cps[e[0]]->is_paired||
       !msc->m_cps[e[1]]->is_paired||
       msc->m_cps[e[1]]->pair_idx != e[0]||
       msc->m_cps[e[0]]->pair_idx != e[1])
      throw std::logic_error("failed to ensure that edge forms a sane pairing ");
#endif
  }

  inline void ensure_cp_is_paired(mscomplex_t *msc,uint c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(!msc->m_cps[c]->is_paired)
      throw std::logic_error("failed to ensure cell is paired ");

    ensure_pairing(msc,uint_pair_t(c,msc->m_cps[c]->pair_idx));
#endif
  }

  inline void ensure_cp_is_not_paired(mscomplex_t *msc,uint c)
  {
#ifdef USE_ENSURE_PREDICATES
    if(msc->m_cps[c]->is_paired)
      throw std::logic_error("failed to ensure cell is not paired ");
#endif
  }
}
#endif

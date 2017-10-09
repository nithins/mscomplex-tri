#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include <trimesh_dataset.h>
#include <trimesh_mscomplex.h>

using namespace std;

namespace br    = boost::range;
namespace badpt = boost::adaptors;

namespace trimesh
{
  dataset_t::dataset_t
  (const fn_list_t &vert_fns, const tri_idx_list_t &trilist):
    m_vert_fns(vert_fns),m_tcc(new tri_cc_t)
  {
    m_tcc->init(trilist,vert_fns.size());

    int N = m_tcc->get_num_cells();

    m_cell_own.resize(N,invalid_cellid);
    m_cell_mxfct.resize(N,invalid_cellid);
    m_cell_pairs.resize(N,invalid_cellid);
  }

  dataset_t::dataset_t
  (const fn_list_t &vert_fns,tri_cc_ptr_t tcc):
    m_vert_fns(vert_fns),m_tcc(tcc)
  {
    int N = m_tcc->get_num_cells();

    m_cell_own.resize(N,invalid_cellid);
    m_cell_mxfct.resize(N,invalid_cellid);
    m_cell_pairs.resize(N,invalid_cellid);
  }


  dataset_t::~dataset_t ()
  {
    m_cell_own.clear();
    m_cell_mxfct.clear();
    m_cell_pairs.clear();
//    m_tcc->clear();
  }

  template <int dim,typename Titer>
  inline void assign_max_facets(dataset_t &ds,Titer b,Titer e)
  {
    cellid_t f[10];

    BOOST_AUTO(cmp,bind(&dataset_t::compare_cells<dim-1>,&ds,_1,_2));

    for(;b!=e;++b)
      ds.max_fct(*b) = *max_element(f,f+ds.get_cets<DES>(*b,f),cmp);
  }

  inline cellid_t * filter_elst(cellid_t *b,cellid_t *e, cellid_t *r, cellid_t c,const dataset_t &ds)
  {for(;b!=e; ++b) if( ds.max_fct(*b) == c) *r++ = *b; return r;}

  template <int dim,typename Titer>
  inline void assign_pairs(dataset_t &ds,Titer b,Titer e)
  {
	const int size = 20;
	cellid_t cf[size], *cfe;

    BOOST_AUTO(cmp,bind(&dataset_t::compare_cells<dim+1>,&ds,_1,_2));

    for(;b!=e;++b)
    {
      uint pos = ds.get_cets<ASC>(*b, cf);
	  ENSURES(pos < size)
			  <<"A maximum of 20 cofacets are assumed. Adjust const above!!";
      cfe = cf + pos;

      cfe = filter_elst(cf,cfe,cf,*b,ds);

      cellid_t *mcf = min_element(cf,cfe,cmp);

      if( mcf != cfe && ds.is_boundry(*mcf) == ds.is_boundry(*b))
        ds.pair(*b,*mcf);
    }
  }

  inline cellid_t * filter_elst2(cellid_t *b,cellid_t *e, cellid_t *r, cellid_t c,const dataset_t &ds)
  {
    for(;b!=e; ++b)
      if( !ds.is_paired(*b) && ds.max_fct(*b) != c && ds.max_fct(ds.max_fct(*b)) == ds.max_fct(c))
        *r++ = *b;

    return r;
  }

  template <int dim,typename Titer>
  inline void assign_pairs2(dataset_t &ds,Titer b,Titer e)
  {
	const int size = 20;
	cellid_t cf[size],*cfe;

    BOOST_AUTO(cmp,bind(&dataset_t::compare_cells<dim+1>,&ds,_1,_2));

    int np = 0;

    for(;b!=e;++b)
      if(!ds.is_paired(*b))
      {
		uint pos = ds.get_cets<ASC>(*b, cf);
		ENSURES(pos < size)
				  <<"A maximum of 20 cofacets are assumed. Adjust const above!!";
		cfe = cf + pos;
		cfe = filter_elst2(cf,cfe,cf,*b,ds);

        cellid_t *mcf = min_element(cf,cfe,cmp);

        if( mcf != cfe && ds.is_boundry(*mcf) == ds.is_boundry(*b))
        {
          ds.pair(*b,*mcf); np++;
        }
      }
  }

  template<typename Toi,typename Tii>
  inline Toi collect_cps(const dataset_t &ds,Tii b,Tii e,Toi r)
  {
    for(; b!=e; ++b)
      if(!ds.is_paired(*b))
        *r++ = *b;

    return r;
  }

  template <eGDIR dir>
  inline void bfs_owner_extrema(dataset_t &ds,cellid_t s)
  {
    const int dim = (dir == DES)?(gc_max_cell_dim):(0);

    std::stack<cellid_t> stk;
    stk.push(s);
    cellid_t f[20],*fe,*fb;

    ASSERT(ds.cell_dim(s) == dim);

    while(!stk.empty())
    {
      cellid_t c = stk.top(); stk.pop();

      ds.owner(c) = s;

      fb = f; fe = f + ds.get_cets<dir>(c,f);

      for (; fb != fe; ++fb )
        if( ds.is_paired(*fb))
        {
          cellid_t p = ds.pair(*fb);
          if(p != c && ds.cell_dim(p) == dim)
            stk.push(p);
        }
    }
  }

  inline void make_connections(mscomplex_t &msc,const cellid_list_t &ccells,const dataset_t &ds)
  {
    msc.resize(ccells.size());

    map<cellid_t,int> id_cp_map;

    for( int i = 0 ; i < ccells.size() ; ++i)
    {
      cellid_t c = ccells[i];
      msc.set_critpt(i,c,ds.cell_dim(c),ds.fn<dataset_t::CFI_MAX>(c),ds.max_vert<-1>(c),ds.is_boundry(c));
      id_cp_map[c] = i;
    }

    cellid_t f[20]; int f_ct;

    for(cellid_list_t::const_iterator b = ccells.begin(),e =ccells.end();b!=e; ++b)
      if(ds.cell_dim(*b) == 1)
      {
        ASSERT(id_cp_map.count(*b) ==1);
        ds.get_cets<DES>(*b,f);

        ASSERT(id_cp_map.count(ds.owner(f[0])) == 1);
        msc.connect_cps(id_cp_map[*b],id_cp_map[ds.owner(f[0])]);

        ASSERT(id_cp_map.count(ds.owner(f[1])) == 1);
        msc.connect_cps(id_cp_map[*b],id_cp_map[ds.owner(f[1])]);

        f_ct = ds.get_cets<ASC>(*b,f);

        ASSERT(id_cp_map.count(ds.owner(f[0])) == 1);
        msc.connect_cps(id_cp_map[*b],id_cp_map[ds.owner(f[0])]);

        if(f_ct == 1) continue;

        ASSERT(id_cp_map.count(ds.owner(f[1])) == 1);
        msc.connect_cps(id_cp_map[*b],id_cp_map[ds.owner(f[1])]);
      }
  }

  void dataset_t::work(mscomplex_ptr_t msc)
  {
    assign_max_facets<1>(*this,m_tcc->begin(1),m_tcc->end(1));
    assign_max_facets<2>(*this,m_tcc->begin(2),m_tcc->end(2));

    assign_pairs<0>(*this,m_tcc->begin(0),m_tcc->end(0));
    assign_pairs<1>(*this,m_tcc->begin(1),m_tcc->end(1));

    assign_pairs2<1>(*this,m_tcc->begin(1),m_tcc->end(1));

    cellid_list_t ccells;

    collect_cps(*this,m_tcc->begin(),m_tcc->end(),back_inserter(ccells));

    for(cellid_list_t::iterator b = ccells.begin(),e =ccells.end();b!=e; ++b)
    {
      if(cell_dim(*b) == 2) bfs_owner_extrema<DES>(*this,*b);
      if(cell_dim(*b) == 0) bfs_owner_extrema<ASC>(*this,*b);
    }

    make_connections(*msc,ccells,*this);
  }

//  template<typename T>
//  inline void bin_write_vec(std::ostream &os, std::vector<T> &v)
//  {os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));}

//  template<typename T>
//  inline void bin_write(std::ostream &os, const T &v)
//  {os.write((const char*)(const void*)&v,sizeof(T));}

//  int get_header_size(int num_cps)
//  {
//    return sizeof(int)              + // num_cps
//           sizeof(cellid_t)*num_cps + // cellids
//           sizeof(char)*num_cps     + // cell indices
//           sizeof(int)*(2*num_cps);   // numcells
//  }

//  void write_header
//  (std::ostream & os,
//   int_list_t & nmcells,
//   cellid_list_t &cps,
//   char_list_t &cp_inds,
//   ios::off_type hoff= 0)
//  {
//    os.seekp(hoff,ios::beg);

//    bin_write(os,(int)cps.size());
//    bin_write_vec(os,cps);
//    bin_write_vec(os,cp_inds);
//    bin_write_vec(os,nmcells);
//  }


//  void  dataset_t::save_manifolds(std::ostream &os,mscomplex_ptr_t msc)
//  {
//    auto rng = msc->cp_range()|badpt::filtered
//        (boost::bind(&mscomplex_t::is_not_paired,msc,_1));

//    int num_cps = utls::count(boost::begin(rng),boost::end(rng));
//    int hoff    = os.tellp();

//    os.seekp(get_header_size(num_cps),ios::cur);

//    int_list_t nmcells;

//    for(auto i = boost::begin(rng);i!=boost::end(rng);++i)
//    {
//      cellid_list_t des,asc;
//      int_list_t    desop,ascop;

//      if(msc->index(*i) != 0) bfs_collect<DES>(*this,*msc,*i,des);
//      if(msc->index(*i) != 2) bfs_collect<ASC>(*this,*msc,*i,asc);

//      BOOST_AUTO(desop_bi,back_inserter(desop));
//      BOOST_AUTO(ascop_bi,back_inserter(ascop));

//      BOOST_FOREACH(cellid_t c,des) m_tcc->cellid_to_output(c,desop_bi);
//      BOOST_FOREACH(cellid_t c,asc) m_tcc->cellid_to_output(c,ascop_bi);

//      nmcells.push_back(desop.size());
//      nmcells.push_back(ascop.size());

//      bin_write_vec(os,desop);
//      bin_write_vec(os,ascop);
//    }

//    cellid_list_t cps;
//    char_list_t   cp_idxs;

//    br::copy(rng,back_inserter(cps));
//    br::transform(rng,back_inserter(cp_idxs),bind(&mscomplex_t::index,msc,_1));

//    write_header(os,nmcells,cps,cp_idxs,hoff);
//  }

}

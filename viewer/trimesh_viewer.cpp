#include <sstream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string_regex.hpp>
#include <boost/static_assert.hpp>

#include <GL/glew.h>

#include <glutils.h>
#include <GLSLProgram.h>

#include <trimesh_viewer.h>
#include <trimesh_mscomplex.h>

#include <config.h>

#ifndef VIEWER_RENDER_AWESOME
double g_max_cp_size  = 8.0;
#else
double g_max_cp_size  = 0.025;
#endif

double g_max_cp_raise = 0.1;

template<typename T> std::string to_string(const T & t)
{
  std::stringstream ss;
  ss<<t;
  return ss.str();
}

using namespace glutils;
using namespace std;

template <typename iter_t>
inline void log_range ( iter_t begin,iter_t end,const char * title = NULL,bool with_ind = false )
{

  if ( title != NULL )
  {
    std::string title_str ( title );
    std::string ln1_str ( title_str.size(),'=' );
    std::string ln2_str ( title_str.size(),'-' );

    std::cout<< ( ln1_str );
    std::cout<< ( title_str );
    std::cout<< ( ln2_str );
  }

  std::stringstream ss;

  unsigned int i = 0;

  for ( iter_t pos = begin;pos !=end; ++pos )
  {
    std::stringstream ss_temp;

    if ( with_ind )
    {
      ss_temp<<i++<<":";
    }

    ss_temp << *pos<<" ";

    ss<<ss_temp.str();
  }
  std::cout<< ( ss.str() );

  if ( title != NULL )
  {
    std::string title_str ( title );
    std::string ln1_str ( title_str.size(),'=' );
    std::cout<< ( ln1_str );
  }

  std::cout<<std::endl;
}

namespace trimesh
{
glutils::color_t g_cp_colors[gc_max_cell_dim+1] =
{
  glutils::color_t(0.0,0.0,1.0),
  glutils::color_t(0.0,1.0,0.0),
  glutils::color_t(1.0,0.0,0.0),
};

glutils::color_t g_grad_colors[gc_max_cell_dim] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
};

glutils::color_t g_disc_colors[GDIR_CT][gc_max_cell_dim+1] =
{
  {
    glutils::color_t(0.15,0.45,0.35 ),
    glutils::color_t(0.85,0.65,0.75 ),
    glutils::color_t(0.0,0.0,0.0 ),
  },

{
    glutils::color_t(0.0,0.0,0.0 ),
    glutils::color_t(0.65,0.95,0.45 ),
    glutils::color_t(0.15,0.25,0.75 ),
  },
};

glutils::color_t g_cp_conn_colors[gc_max_cell_dim] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
};

glutils::color_t g_roiaabb_color = glutils::color_t(0.85,0.75,0.65);

glutils::color_t g_normals_color = glutils::color_t(0.85,0.75,0.65);

viewer_t::viewer_t
    (std::string tf,std::string gf,std::string mf):
    m_data_dia(0),
    m_bShowSurface(false)

{
  m_graphs.push_back(new graph_data_t(tf,gf,mf));
}

viewer_t::~viewer_t()
{
  for ( uint i = 0 ; i < m_graphs.size();i++ )
    delete m_graphs[i];

  m_graphs.clear();
}

int viewer_t::render()
{
  glPushAttrib(GL_ENABLE_BIT);

  glEnable(GL_RESCALE_NORMAL);

  glScalef(1.0/m_data_dia,1.0/m_data_dia,1.0/m_data_dia);

  for ( uint i = 0 ; i < m_graphs.size();i++ )
  {
    m_graphs[i]->render();
  }

  if(m_bShowSurface)
    m_surf_ren->render();

  glPopAttrib();
}

void viewer_t::init()
{
  glutils::init();

  graph_data_t *gd= m_graphs[0];

  gd->init();

  double *e = gd->m_mfold_data.m_extent;

  vertex_t s(e[1]-e[0],e[3]-e[2],e[5]-e[4]);

  m_data_dia = *std::max_element(s.begin(),s.end());
}

configurable_t::data_index_t viewer_t::dim()
{
  return data_index_t(11,m_graphs.size());
}

bool viewer_t::exchange_field(const data_index_t &idx, boost::any &v)
{
  graph_data_t * gd = m_graphs[idx[1]];

  switch(idx[0])
  {
  case 0: return s_exchange_data_ro(std::string("0"),v);
  case 1: return s_exchange_data_rw(gd->m_bShowAllCps,v);
  case 2: return s_exchange_data_rw(gd->m_bShowCps[0],v);
  case 3: return s_exchange_data_rw(gd->m_bShowCps[1],v);
  case 4: return s_exchange_data_rw(gd->m_bShowCps[2],v);
  case 5: return s_exchange_data_rw(gd->m_bShowCpLabels,v);
  case 6: return s_exchange_data_rw(gd->m_bShowMsGraph,v);
  case 7: return s_exchange_data_rw(gd->m_bShowGrad,v);
  case 8: return s_exchange_data_rw(gd->m_bShowCancCps,v);
  case 9: return s_exchange_data_rw(gd->m_bShowCancMsGraph,v);
  case 10: return s_exchange_data_rw(gd->m_bShowCellNormals,v);
  }

  throw std::logic_error("unknown index");
}
configurable_t::eFieldType viewer_t::exchange_header(const int &i, boost::any &v)
{
  switch(i)
  {
  case 0: v = std::string("oct tree piece"); return EFT_DATA_RO;
  case 1: v = std::string("all cps"); return EFT_DATA_RW;
  case 2: v = std::string("minima"); return EFT_DATA_RW;
  case 3: v = std::string("1 saddle"); return EFT_DATA_RW;
  case 4: v = std::string("maxima"); return EFT_DATA_RW;
  case 5: v = std::string("cp labels"); return EFT_DATA_RW;
  case 6: v = std::string("msgraph");return EFT_DATA_RW;
  case 7: v = std::string("gradient");return EFT_DATA_RW;
  case 8: v = std::string("cancelled cps");return EFT_DATA_RW;
  case 9: v = std::string("cancelled cp msgraph");return EFT_DATA_RW;
  case 10: v = std::string("cell normals");return EFT_DATA_RW;
  }
  throw std::logic_error("unknown index");
}

graph_data_t::graph_data_t(std::string tf, std::string gf, std::string mf):
    m_bShowAllCps(false),
    m_bShowCpLabels ( false ),
    m_bShowMsGraph ( false ),
    m_bShowGrad ( false ),
    m_bShowCancCps(false),
    m_bShowCancMsGraph(false),
    m_bShowCellNormals(false),
    m_msc(new mscomplex_t),
    m_mfold_data(mf,tf,m_msc)
{
  m_msc->load(gf);

  m_bShowCps[0] = false;
  m_bShowCps[1] = false;
  m_bShowCps[2] = false;
}

void graph_data_t::init()
{
  m_mfold_data.init();

  init_cps();
  init_ccps();
}

void graph_data_t::init_cps()
{
  using namespace glutils;

  point_idx_list_t  vind[gc_max_cell_dim+1];
  line_idx_list_t   eind[gc_max_cell_dim];

  int N = m_msc->get_num_critpts();

  for(uint i = 0; i < N; ++i)
    if(!m_msc->is_paired(i))
      vind[m_msc->index(i)].push_back(m_msc->cellid(i));

  bufobj_ptr_t cbo = m_mfold_data.m_cell_pos_bo;

  for(uint i = 0 ; i < N; ++i)
  if(!m_msc->is_paired(i))
    for(conn_iter_t b = m_msc->m_des_conn[i].begin(),e = m_msc->m_des_conn[i].end(); b != e; ++b)
      eind[m_msc->index(i)-1].push_back(line_idx_t(m_msc->cellid(*b),m_msc->cellid(i)));

  ren_cp[0].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[0])));
  ren_cp[1].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[1])));
  ren_cp[2].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[2])));

  ren_cp_conns[0].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[0])));
  ren_cp_conns[1].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[1])));
}

void graph_data_t::init_ccps()
{
  using namespace glutils;

  point_idx_list_t  vind[gc_max_cell_dim+1];
  line_idx_list_t   eind[gc_max_cell_dim];

  int N = m_msc->get_num_critpts();

  for(uint i = 0; i < N; ++i)
    if(m_msc->is_paired(i))
        vind[m_msc->index(i)].push_back(m_msc->cellid(i));

  for(uint i = 0 ; i < N; ++i)
    if(m_msc->is_paired(i))
    {
      int idx  = m_msc->index(i), pidx = m_msc->index(m_msc->pair_idx(i));
      int eidx = min(idx,pidx);
      eGDIR d  = (idx > pidx)?(GDIR_DES):(GDIR_ASC);

      for(conn_iter_t b = m_msc->m_conn[d][i].begin(),e = m_msc->m_conn[d][i].end();b!=e ;++b)
        eind[eidx].push_back(line_idx_t(m_msc->cellid(*b),m_msc->cellid(i)));
    }

  bufobj_ptr_t cbo = m_mfold_data.m_cell_pos_bo;

  ren_canc_cp[0].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[0])));
  ren_canc_cp[1].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[1])));
  ren_canc_cp[2].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[2])));

  ren_canc_cp_conns[0].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[0])));
  ren_canc_cp_conns[1].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[1])));
}

void graph_data_t::render()
{
  glPushMatrix();
  glPushAttrib ( GL_ENABLE_BIT );

  glDisable ( GL_LIGHTING );

#ifndef VIEWER_RENDER_AWESOME
  glPointSize ( g_max_cp_size );

  glEnable(GL_POINT_SMOOTH);
#else
  g_sphere_shader->use();

  g_sphere_shader->sendUniform("g_wc_radius",float(g_max_cp_size*m_viewer->m_data_dia));
#endif

  for(uint i = 0 ; i < gc_max_cell_dim+1;++i)
  {
    if(ren_cp[i]&& (m_bShowCps[i]||m_bShowAllCps))
    {
      glColor3dv(g_cp_colors[i].data());

      ren_cp[i]->render();

      if(ren_cp_labels[i] && m_bShowCpLabels)
        ren_cp_labels[i]->render();
    }
  }

#ifdef VIEWER_RENDER_AWESOME
  g_sphere_shader->sendUniform("g_wc_radius",(float)g_max_cp_size*2/3);
#endif

  if ( m_bShowCancCps)
  {
    for(uint i = 0 ; i < gc_max_cell_dim+1;++i)
    {
      if(ren_canc_cp[i])
      {
        glColor3dv(g_cp_colors[i].data());

        ren_canc_cp[i]->render();

        if(ren_canc_cp_labels[i] &&
           m_bShowCpLabels)
          ren_canc_cp_labels[i]->render();
      }
    }
  }

#ifdef VIEWER_RENDER_AWESOME
  g_sphere_shader->disable();
#endif

  if (m_bShowMsGraph)
  {
    for(uint i = 0 ; i < gc_max_cell_dim;++i)
    {
      if(ren_cp_conns[i])
      {
        glColor3dv(g_cp_conn_colors[i].data());

        ren_cp_conns[i]->render();
      }
    }
  }

  if (m_bShowCancMsGraph)
  {
    for(uint i = 0 ; i < gc_max_cell_dim;++i)
    {
      if(ren_canc_cp_conns[i])
      {
        glColor3dv(g_cp_conn_colors[i].data());

        ren_canc_cp_conns[i]->render();
      }
    }
  }

  glPopAttrib();
  glPopMatrix();

  m_mfold_data.render();
}

mfold_data_t::mfold_data_t(std::string mf, std::string tf, mscomplex_ptr_t msc)
  :m_mfold_file(mf),m_msc(msc),m_tri_file(tf),m_tcc(new tri_cc_geom_t),m_bNeedUpdate(false)
{}

template<typename T>
inline void bin_read_vec(std::istream &is, std::vector<T> &v,int n)
{v.resize(n);is.read((char*)(void*)v.data(),n*sizeof(T));}

template<typename T>
inline void bin_read(std::istream &is, const T &v)
{is.read((char*)(void*)&v,sizeof(T));}

int get_header_size(int num_cps)
{
  return sizeof(int)              + // num_cps
         sizeof(cellid_t)*num_cps + // cellids
         sizeof(int)*(2*num_cps+1); // offsets
}

void read_header(std::istream & is,int_list_t & offsets,cellid_list_t &cps,ios::off_type hoff= 0)
{
  is.seekg(hoff,ios::beg);

  int N;

  bin_read(is,N);
  bin_read_vec(is,cps,N);
  bin_read_vec(is,offsets,2*N+1);
}

void mfold_data_t::__read_mfold(int i,cellid_list_t &mfold)
{
  ASSERT(is_in_range(i,0,m_offsets.size()-1));

  fstream fs(m_mfold_file.c_str(),ios::in|ios::binary);
  ensure(fs.is_open(),"Cannot open file");

  int N = m_offsets[i+1]  - m_offsets[i];

  fs.seekg(get_header_size(m_cellids.size())+m_offsets[i]);
  bin_read_vec(fs,mfold,N);
}

void mfold_data_t::init()
{
  fstream fs(m_mfold_file.c_str(),ios::in|ios::binary);
  ensure(fs.is_open(),"Cannot open file");
  read_header(fs,m_offsets,m_cellids);

  m_scpno_cpno_map.resize(m_cellids.size());

  for( int i = 0,j=0 ; i < m_msc->get_num_critpts(); ++i)
    if(!m_msc->is_paired(i))
    {
      ASSERT(m_cellids[j] == m_msc->cellid(i));
      m_scpno_cpno_map[j++] = i;
    }


  tri_idx_list_t tlist; vertex_list_t vlist;
  read_tri_file(m_tri_file.c_str(),vlist,tlist);
  m_tcc->init(tlist,vlist);
  m_cell_pos_bo = make_buf_obj(m_tcc->get_cell_positions());
  m_cell_nrm_bo = make_buf_obj(m_tcc->get_cell_normals());
  compute_extent(vlist,m_extent);

  m_ren[0].resize(m_cellids.size());
  m_ren[1].resize(m_cellids.size());

  m_ren_show[0].resize(m_cellids.size(),false);
  m_ren_show[1].resize(m_cellids.size(),false);

  m_ren_color[0].resize(m_cellids.size(),g_disc_colors[0][0]);
  m_ren_color[1].resize(m_cellids.size(),g_disc_colors[1][1]);
}

void mfold_data_t::render()
{
  if(m_bNeedUpdate)
  {
    update();
    m_bNeedUpdate = false;
  }

  glPushMatrix();
  glPushAttrib ( GL_ENABLE_BIT );

  for(int d = 0 ; d < 2; ++d)
  for(int i = 0 ; i < m_cellids.size(); ++i)
  if(m_ren_show[d][i] && m_ren[d][i])
  {
    glColor3dv(m_ren_color[d][i].data());
    m_ren[d][i]->render();
  }

  glPopAttrib();
  glPopMatrix();

}

void assign_random_color(glutils::color_t &col)
{
  const uint MAX_RAND = 256;

  col[0] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);
  col[1] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);
  col[2] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);
}

configurable_t::data_index_t mfold_data_t::dim()
{return data_index_t(9,m_cellids.size());}

bool mfold_data_t::exchange_field
    (const data_index_t &idx,boost::any &v)
{
  bool is_read     = v.empty();

  int i = idx[0];

  int scpno = idx[1];
  int cpno = m_scpno_cpno_map[idx[1]];

  switch(i)
  {
  case 0:
    return s_exchange_data_ro((int)m_msc->cellid(cpno),v);
  case 1:
    return s_exchange_data_ro((int)m_msc->index(cpno),v);
  case 2:
  case 3:
    {
      bool show = m_ren_show[i%2][scpno];

      bool need_update = s_exchange_data_rw(show,v);

      m_ren_show[i%2][scpno] = show;

      if(need_update && !is_read)
        m_bNeedUpdate = true;

      return need_update;
    }
  case 4:
  case 5:
    return s_exchange_data_rw(m_ren_color[i%2][scpno],v);
  case 6:
  case 7:
    return s_exchange_action(bind(assign_random_color,boost::ref(m_ren_color[i%2][scpno])),v);
  case 8:
    return s_exchange_data_ro((int)m_msc->vertid(cpno),v);
  };

   throw std::logic_error("octtree_piece_rendata::invalid index");
}

configurable_t::eFieldType mfold_data_t::exchange_header
    (const int &i,boost::any &v)
{
  switch(i)
  {
  case 0: v = std::string("cellid"); return EFT_DATA_RO;
  case 1: v = std::string("index"); return EFT_DATA_RO;
  case 2: v = std::string("des mfold"); return EFT_DATA_RW;
  case 3: v = std::string("asc mfold"); return EFT_DATA_RW;
  case 4: v = std::string("des mfold color"); return EFT_DATA_RW;
  case 5: v = std::string("asc mfold color"); return EFT_DATA_RW;
  case 6: v = std::string("rand des mflod color"); return EFT_ACTION;
  case 7: v = std::string("rand asc mflod color"); return EFT_ACTION;
  case 8: v = std::string("vert no"); return EFT_DATA_RO;

  }
  throw std::logic_error("octtree_piece_rendata::invalid index");
}

template<eGDIR dir>
inline int get_edge_pts(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc);

template<>
inline int get_edge_pts<GDIR_DES>(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc)
{return tcc.get_cell_points(e,pts);}

template<>
inline int get_edge_pts<GDIR_ASC>(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc)
{
  int ncf = tcc.get_cell_co_facets(e,pts);
  pts[2] = pts[1]; pts[1] = e;
  return ncf+1;
}

template<eGDIR dir>
void update_saddle_mfold(mfold_data_t &dp,int scpno)
{
  cellid_list_t mfold;

  dp.read_mfold<dir>(scpno,mfold);

  map<cellid_t,int> pt_idx;

  line_idx_list_t e_idxs;

  for(cellid_list_t::iterator it = mfold.begin(); it!= mfold.end(); ++it)
  {
    cellid_t pt[20];

    dp.m_tcc->get_cell_points(*it,pt);

    int npts = get_edge_pts<dir>(*it,pt,*(dp.m_tcc));

    if( pt_idx.count(pt[0]) == 0) pt_idx[pt[0]] = pt_idx.size()-1;
    if( pt_idx.count(pt[1]) == 0) pt_idx[pt[1]] = pt_idx.size()-1;
    e_idxs.push_back(glutils::line_idx_t(pt_idx[pt[0]],pt_idx[pt[1]]));

    if(dir == GDIR_DES) continue;

    if(npts < 3) continue;

    if(pt_idx.count(pt[2]) == 0) pt_idx[pt[2]] = pt_idx.size()-1;
    e_idxs.push_back(glutils::line_idx_t(pt_idx[pt[1]],pt_idx[pt[2]]));
  }

  vertex_list_t pts(pt_idx.size());

  for(map<cellid_t,int>::iterator it = pt_idx.begin(); it!= pt_idx.end();++it)
    pts[it->second] = dp.m_tcc->get_cell_position(it->first);

  smooth_lines(pts,e_idxs,4);

  dp.m_ren[dir][scpno].reset(create_buffered_lines_ren(make_buf_obj(pts),make_buf_obj(e_idxs)));
}

template<eGDIR dir>
void update_extrema_mfold(mfold_data_t &dp,int scpno)
{
  cellid_list_t mfold;

  dp.read_mfold<dir>(scpno,mfold);

  tri_idx_list_t t_idxs;

  for(cellid_list_t::iterator it = mfold.begin(); it!= mfold.end(); ++it)
  {
    cellid_t st[40];

    if(dir == GDIR_DES)
    {
      dp.m_tcc->get_cell_points(*it,st);
      t_idxs.push_back(tri_idx_t(st[0],st[1],st[2]));
    }
    else
    {
      uint st_ct = dp.m_tcc->get_vert_star(*it,st);

      for(uint i = 1; i < st_ct; i++)
        t_idxs.push_back(glutils::tri_idx_t(st[i-1],*it,st[i]));

      if(st_ct%2 == 0)
        t_idxs.push_back(glutils::tri_idx_t(st[st_ct-1],*it,st[0]));
    }
  }

  dp.m_ren[dir][scpno].reset
      (create_buffered_triangles_ren(dp.m_cell_pos_bo,make_buf_obj(t_idxs),dp.m_cell_nrm_bo));
}

bool mfold_data_t::update()
{
  for( int i = 0 ; i < m_cellids.size(); ++i)
  {
    switch(m_msc->index(m_scpno_cpno_map[i]))
    {
    case 0:
      if(m_ren_show[0][i] && !m_ren[0][i]) update_extrema_mfold<GDIR_ASC>(*this,i);break;
    case 1:
      if(m_ren_show[0][i] && !m_ren[0][i]) update_saddle_mfold<GDIR_DES>(*this,i);
      if(m_ren_show[1][i] && !m_ren[1][i]) update_saddle_mfold<GDIR_ASC>(*this,i);break;
    case 2:
      if(m_ren_show[1][i] && !m_ren[1][i]) update_extrema_mfold<GDIR_DES>(*this,i);break;
    }

    if(!m_ren_show[0][i] && m_ren[0][i]) m_ren[0][i].reset();
    if(!m_ren_show[1][i] && m_ren[1][i]) m_ren[1][i].reset();
  }
}

}

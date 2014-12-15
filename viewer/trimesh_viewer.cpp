#include <sstream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string_regex.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/static_assert.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/foreach.hpp>


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
GLSLProgram * g_cylinder_shader = NULL;
GLSLProgram * g_sphere_shader   = NULL;

extern const char *CYLINDER_VERT_GLSL;
extern const char *CYLINDER_GEOM_GLSL;
extern const char *CYLINDER_FRAG_GLSL;

extern const char *SPHERE_VERT_GLSL;
extern const char *SPHERE_GEOM_GLSL;
extern const char *SPHERE_FRAG_GLSL;
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
namespace br    = boost::range;
namespace badpt = boost::adaptors;
namespace ba    = boost::algorithm;

template <>
    bool configurable_t::s_exchange_data_rw(bool &p_val,boost::any &c_val)
{
  bool ret = false;

  if(c_val.empty())
    c_val = boost::any(p_val);
  else
  {
    bool newval;

    if(c_val.type() == typeid(bool))
      newval = boost::any_cast<bool>(c_val);
    else
      newval = (bool) boost::any_cast<int>(c_val);

    ret   = (p_val != newval);
    p_val = newval;
  }
  return ret;
}


namespace trimesh
{

template <typename T>
inline tri_idx_t mk_tri_idx(const T& a,const T&b,const T& c)
{tri_idx_t t;t[0] = a; t[1] = b; t[2] = c; return t;}

glutils::color_t g_cp_colors[gc_max_cell_dim+1] =
{
  la::make_vec<double>(0.0,0.0,1.0),
  la::make_vec<double>(0.0,1.0,0.0),
  la::make_vec<double>(1.0,0.0,0.0),
};

glutils::color_t g_grad_colors[gc_max_cell_dim] =
{
  la::make_vec<double>(0.0,0.5,0.5 ),
  la::make_vec<double>(0.5,0.0,0.5 ),
};

glutils::color_t g_disc_colors[GDIR_CT][gc_max_cell_dim+1] =
{
  {
    la::make_vec<double>(0.15,0.45,0.35 ),
    la::make_vec<double>(0.85,0.65,0.75 ),
    la::make_vec<double>(0.0,0.0,0.0 ),
  },

{
    la::make_vec<double>(0.0,0.0,0.0 ),
    la::make_vec<double>(0.65,0.95,0.45 ),
    la::make_vec<double>(0.15,0.25,0.75 ),
  },
};

glutils::color_t g_cp_conn_colors[gc_max_cell_dim] =
{
  la::make_vec<double>(0.0,0.5,0.5 ),
  la::make_vec<double>(0.5,0.0,0.5 ),
};

glutils::color_t g_roiaabb_color = la::make_vec<double>(0.85,0.75,0.65);

glutils::color_t g_normals_color = la::make_vec<double>(0.85,0.75,0.65);


glutils::material_properties_t g_surface_material_default =
{
  la::make_vec<float>(0.2,0.2,0.2,1),// ambient
  la::make_vec<float>(0.5,0.5,0.5,1),// diffuse
  la::make_vec<float>(0.8,0.8,0.8,1),      // specular
  la::make_vec<float>(0,0,0,1),      // emission
  66             // shininess
};


// two directional lights by default
// as a convention directional lights are directed along 0,0,1 model coords
// and positional lights are located at 0,0,0 model coords
// another frame is associated with the light to move it around
// and orient it in world coords
glutils::light_properties_t g_lights_default[] =
{
  {
    la::make_vec<float>(0.2f, 0.2f, 0.2f, 1.0f), // ambient
    la::make_vec<float>(0.5f, 0.5f, 0.5f, 1.0f), // diffuse
    la::make_vec<float>(0.7f, 0.7f, 0.7f, 0.7f), // specular
    la::make_vec<float>(0.0f, 1.0f, 1.0f, 0.0f), // position
    1,0,0,                    // attenution c + l + q
    true                      // enabled
  },

  {
    la::make_vec<float>(0.2f, 0.2f, 0.2f, 1.0f), // ambient
    la::make_vec<float>(0.5f, 0.5f, 0.5f, 1.0f), // diffuse
    la::make_vec<float>(0.7f, 0.7f, 0.7f, 0.7f), // specular
    la::make_vec<float>(0.0f, 0.0f, 1.0f, 0.0f), // position
    1,0,0,                    // attenution c + l + q
    true                      // enabled
  }
};

inline color_t get_random_color()
{
  color_t col;

  const uint MAX_RAND = 256;

  col[0] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);
  col[1] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);
  col[2] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);

  return col;
}

viewer_t::viewer_t()
{
}

viewer_t::~viewer_t()
{
}

int viewer_t::render()
{
  glPushAttrib(GL_ENABLE_BIT);

  glEnable(GL_RESCALE_NORMAL);

  glPushMatrix();

  GLfloat lmodel_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);

  g_lights_default[0].render(0);
  g_lights_default[1].render(1);
  g_surface_material_default.render_all(GL_FRONT_FACE);

  glEnable(GL_COLOR_MATERIAL);
  glColorMaterial(GL_FRONT,GL_AMBIENT_AND_DIFFUSE);

  glPopMatrix();

  glPushMatrix();

  for( int i = 0 ; i < m_mscs.size(); ++i)
    m_mscs[i]->render();

  glPopMatrix();

  glPopAttrib();
}

void viewer_t::init()
{
  glutils::init();
  
#ifdef VIEWER_RENDER_AWESOME

  g_cylinder_shader = GLSLProgram::createFromSourceStrings
                      (CYLINDER_VERT_GLSL,CYLINDER_GEOM_GLSL,
                       CYLINDER_FRAG_GLSL);

  g_sphere_shader   = GLSLProgram::createFromSourceStrings
                      (SPHERE_VERT_GLSL,SPHERE_GEOM_GLSL,
                       SPHERE_FRAG_GLSL);
#endif
}

configurable_t::data_index_t viewer_t::dim()
{
  return data_index_t(11,m_mscs.size());
}

bool viewer_t::exchange_field(const data_index_t &idx, boost::any &v)
{
  mscomplex_ren_ptr_t gd = m_mscs[idx[1]];

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
configurable_t::eFieldType viewer_t::exchange_header(const int &i, std::string &v)
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

mscomplex_ren_t::mscomplex_ren_t(std::string tf, std::string mf):
    m_bShowAllCps(false),
    m_bShowCpLabels ( false ),
    m_bShowMsGraph ( false ),
    m_bShowGrad ( false ),
    m_bShowCancCps(false),
    m_bShowCancMsGraph(false),
    m_bShowCellNormals(false),
    m_msc(new mscomplex_t),
    m_tcc(new tri_cc_geom_t),
    m_need_update_geom(false)

{
  m_canc_tresh = 0;

  m_msc->load(mf);

  m_bShowCps[0] = false;
  m_bShowCps[1] = false;
  m_bShowCps[2] = false;

  tri_idx_list_t tlist;
  vertex_list_t  vlist;
  string ext = tf.substr(tf.size()-4,4);

  ba::to_lower(ext);

  if(ext == ".off")
    read_off_file(tf.c_str(),vlist,tlist);
  else if(ext == ".tri")
    read_tri_file(tf.c_str(),vlist,tlist);
  else
    ensure(false,"Cannot read triangulation");

  m_tcc->init(tlist,vlist);
  compute_extent(vlist,m_extent);
  m_center = compute_center(vlist);
}

void mscomplex_ren_t::init()
{
  using namespace glutils;

  m_cell_pos_bo = make_buf_obj(m_tcc->get_cell_positions());
  m_cell_nrm_bo = make_buf_obj(m_tcc->get_cell_normals());

  point_idx_list_t  vind[gc_max_cell_dim+1];
  line_idx_list_t   eind[gc_max_cell_dim];

  for( int i = 0 ; i < m_msc->get_num_critpts(); ++i)
  {
    if(m_msc->is_not_paired(i))
    {
      m_surv_cp_rev[i] = m_surv_cps.size();
      m_surv_cps.push_back(i);
    }
  }

  BOOST_FOREACH(int i,m_surv_cps)
      vind[m_msc->index(i)].push_back(m_msc->cellid(i));

  bufobj_ptr_t cbo = m_cell_pos_bo;

  BOOST_FOREACH(int i,m_surv_cps)
  {
    BOOST_FOREACH(int j,m_msc->m_des_conn[i])
    {
      line_idx_t l = la::make_vec<idx_t>(m_msc->cellid(i),m_msc->cellid(j));
      eind[m_msc->index(j)].push_back(l);
    }
  }

  ren_cp[0].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[0])));
  ren_cp[1].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[1])));
  ren_cp[2].reset(create_buffered_points_ren(cbo,make_buf_obj(vind[2])));

  ren_cp_conns[0].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[0])));
  ren_cp_conns[1].reset(create_buffered_lines_ren(cbo,make_buf_obj(eind[1])));

  m_surv_mfold_rens[0].resize(m_surv_cps.size());
  m_surv_mfold_rens[1].resize(m_surv_cps.size());

  m_surv_mfold_show[0].resize(m_surv_cps.size(),false);
  m_surv_mfold_show[1].resize(m_surv_cps.size(),false);

  m_surv_mfold_color[0].resize(m_surv_cps.size(),g_disc_colors[0][0]);
  m_surv_mfold_color[1].resize(m_surv_cps.size(),g_disc_colors[1][1]);

  m_scale_factor = max(m_extent[1] - m_extent[0],m_extent[3]-m_extent[2]);
  m_scale_factor = max(m_extent[5] - m_extent[4],m_scale_factor);

  m_scale_factor = 2.0/m_scale_factor;

}

inline int get_root(int i,const mscomplex_ren_t::canc_tree_t &canc_tree)
{
  while(canc_tree[i].parent != -1)
    i = canc_tree[i].parent;

  return i;
}

inline int get_ancestor(int i,const mscomplex_ren_t::canc_tree_t &canc_tree,double tresh)
{
  while(canc_tree[i].parent != -1 && canc_tree[i].perst < tresh)
    i = canc_tree[i].parent;

  return i;
}

inline int get_other_ex(int_pair_t pr,mscomplex_ptr_t msc,const mscomplex_ren_t::canc_tree_t &canc_tree)
{
  if(!msc->is_extrema(pr.first)) swap(pr.first,pr.second);

  int dir = (msc->index(pr.first) == 0)?(0):(1);

  auto beg = msc->m_conn[dir][pr.second].begin();
  auto end = msc->m_conn[dir][pr.second].end();

  ASSERT(beg != end);

  int ex1  = *beg++;
  int ex2  = (beg == end)?(-1):(*beg++);

  ASSERT(beg == end);
  ASSERT(is_in_range(ex1,0,msc->get_num_critpts()));

  int ex = get_root(ex1,canc_tree);

  if(ex == pr.first)
  {
    ASSERT(is_in_range(ex2,0,msc->get_num_critpts()));
    ex = get_root(ex2,canc_tree);
  }

  return ex;
}

void mscomplex_ren_t::build_canctree(const int_pair_list_t & canc_list,
                                     const fn_list_t       & canc_pers)
{
  ASSERT(canc_list.size() == canc_pers.size());

  m_canc_tree.resize(m_msc->get_num_critpts());

  double fmax = *br::max_element(m_msc->m_cp_fn);
  double fmin = *br::min_element(m_msc->m_cp_fn);

  for(int i = 0; i < canc_list.size() ;++i)
  {
    int_pair_t pr = canc_list[i];

    ASSERT(is_in_range(pr.first,0,m_msc->get_num_critpts()));
    ASSERT(is_in_range(pr.second,0,m_msc->get_num_critpts()));

    ASSERT(m_msc->index(pr.first) + 1 == m_msc->index(pr.second) ||
           m_msc->index(pr.second) + 1 == m_msc->index(pr.first));

    if(!m_msc->is_extrema(pr.first)) swap(pr.first,pr.second);

    ASSERT(m_msc->is_extrema(pr.first) && m_msc->is_saddle(pr.second));

    fn_t perst = canc_pers[i];

    m_canc_tree[pr.first].perst  = perst;
    m_canc_tree[pr.second].perst  = perst;

    m_canc_tree[pr.first].parent = get_other_ex(pr,m_msc,m_canc_tree);
  }

  for(int i = 0 ; i < m_msc->get_num_critpts(); ++i)
  {
    m_canc_tree[i].color = get_random_color();
  }
}

void mscomplex_ren_t::update_canctree_tresh(double tresh)
{
  if (m_canc_tree.size() == 0)
    return;

  m_canc_tresh = tresh;

  for( int i = 0; i < m_surv_cps.size(); ++i)
  {
    int j = m_surv_cps[i];

    if(m_msc->is_saddle(j))
      continue;

    m_surv_mfold_color[0][i] = m_canc_tree[get_ancestor(j,m_canc_tree,tresh)].color;
    m_surv_mfold_color[1][i] = m_canc_tree[get_ancestor(j,m_canc_tree,tresh)].color;
  }
}

void mscomplex_ren_t::render()
{
  glPushMatrix();
  glPushAttrib ( GL_ENABLE_BIT );

  glEnable(GL_LIGHTING);

  glScalef(m_scale_factor,m_scale_factor,m_scale_factor);
  glTranslated(-m_center[0],-m_center[1],-m_center[2]);

  if(m_need_update_geom)
  {
    update_geom();
    m_need_update_geom = false;
  }

  for(int d = 0 ; d < 2; ++d)
  for(int i = 0 ; i < m_surv_cps.size(); ++i)
  if(m_surv_mfold_show[d][i] && m_surv_mfold_rens[d][i])
  {
    glColor3dv(&m_surv_mfold_color[d][i][0]);

#ifdef VIEWER_RENDER_AWESOME
    if(m_msc->index(m_surv_cps[i]) == 1)
    {
      g_cylinder_shader->use();
      g_cylinder_shader->sendUniform("ug_cylinder_radius",float(0.05/(4*m_scale_factor)));
    }
#endif

    m_surv_mfold_rens[d][i]->render();

#ifdef VIEWER_RENDER_AWESOME
    g_cylinder_shader->disable();
#endif

  }

#ifndef VIEWER_RENDER_AWESOME
    glPointSize ( g_max_cp_size );

    glEnable(GL_POINT_SMOOTH);
#else
    g_sphere_shader->use();

    g_sphere_shader->sendUniform("g_wc_radius",float(g_max_cp_size/m_scale_factor));
#endif

  m_cell_pos_bo->bind_to_vertex_pointer();
  glBegin(GL_POINTS);

  BOOST_FOREACH(int i,m_cp_ren_set)
  {
    glColor3dv(&g_cp_colors[m_msc->index(m_surv_cps[i])][0]);
    glArrayElement(m_msc->cellid(m_surv_cps[i]));
  }

  BOOST_FOREACH(int i,m_surv_cps)
  {
    int indx = m_msc->index(i);
    bool surv = (m_canc_tree.size() == 0) ||
        (get_ancestor(i,m_canc_tree,m_canc_tresh) == i);

    if((m_bShowCps[indx]||m_bShowAllCps) && surv)
    {
      glColor3dv(&g_cp_colors[indx][0]);
      glArrayElement(m_msc->cellid(i));
    }
  }

  glEnd();

  m_cell_pos_bo->unbind_from_vertex_pointer();


//  for(uint i = 0 ; i < gc_max_cell_dim+1;++i)
//  {
//    if(ren_cp[i]&& (m_bShowCps[i]||m_bShowAllCps))
//    {
//      glColor3dv(&g_cp_colors[i][0]);

//      ren_cp[i]->render();

//      if(ren_cp_labels[i] && m_bShowCpLabels)
//        ren_cp_labels[i]->render();
//    }
//  }

#ifdef VIEWER_RENDER_AWESOME
  g_sphere_shader->sendUniform("g_wc_radius",(float)g_max_cp_size*2/3);
#endif

  if ( m_bShowCancCps)
  {
    for(uint i = 0 ; i < gc_max_cell_dim+1;++i)
    {
      if(ren_canc_cp[i])
      {
        glColor3dv(&g_cp_colors[i][0]);

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
        glColor3dv(&g_cp_conn_colors[i][0]);

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
        glColor3dv(&g_cp_conn_colors[i][0]);

        ren_canc_cp_conns[i]->render();
      }
    }
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

configurable_t::data_index_t mscomplex_ren_t::dim()
{return data_index_t(9,m_surv_cps.size());}

}

namespace glutils
{
bool operator!=(const color_t &c1 ,const color_t &c2)
{return !(c1[0] == c2[0] &&c1[1] == c2[1] &&c1[2] == c2[2]);}
}

namespace trimesh
{

bool mscomplex_ren_t::exchange_field
    (const data_index_t &idx,boost::any &v)
{
  bool is_read     = v.empty();

  int i = idx[0];

  int scpno = idx[1];
  int cpno = m_surv_cps[scpno];

  switch(i)
  {
  case 0:
    return s_exchange_data_ro((int)m_msc->cellid(cpno),v);
  case 1:
    return s_exchange_data_ro((int)m_msc->index(cpno),v);
  case 2:
  case 3:
    {
      bool show = m_surv_mfold_show[i%2][scpno];

      bool need_update = s_exchange_data_rw(show,v);

      m_surv_mfold_show[i%2][scpno] = show;

      if(need_update && !is_read)
        m_need_update_geom = true;

      return need_update;
    }
  case 4:
  case 5:
    return s_exchange_data_rw(m_surv_mfold_color[i%2][scpno],v);
  case 6:
  case 7:
    return s_exchange_action(bind(assign_random_color,boost::ref(m_surv_mfold_color[i%2][scpno])),v);
  case 8:
    return s_exchange_data_ro((int)m_msc->vertid(cpno),v);
  };

   throw std::logic_error("octtree_piece_rendata::invalid index");
}

configurable_t::eFieldType mscomplex_ren_t::exchange_header
    (const int &i,std::string &v)
{
  switch(i)
  {
  case 0: v = std::string("cellid"); return EFT_DATA_RO;
  case 1: v = std::string("index"); return EFT_DATA_RO;
  case 2: v = std::string("des mfold"); return EFT_DATA_RW;
  case 3: v = std::string("asc mfold"); return EFT_DATA_RW;
  case 4: v = std::string("des mfold color"); return EFT_DATA_RW;
  case 5: v = std::string("asc mfold color"); return EFT_DATA_RW;
  case 6: v = std::string("rand des mfold color"); return EFT_ACTION;
  case 7: v = std::string("rand asc mfold color"); return EFT_ACTION;
  case 8: v = std::string("vert no"); return EFT_DATA_RO;

  }
  throw std::logic_error("octtree_piece_rendata::invalid index");
}

template<eGDIR dir>
inline int get_edge_pts(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc);

template<>
inline int get_edge_pts<DES>(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc)
{return tcc.get_cell_points(e,pts);}

template<>
inline int get_edge_pts<ASC>(cellid_t e,cellid_t *pts,const tri_cc_geom_t &tcc)
{
  int ncf = tcc.get_cell_co_facets(e,pts);
  pts[2] = pts[1]; pts[1] = e;
  return ncf+1;
}

renderable_ptr_t get_maxima_renderer
(mfold_t &mfold,tri_cc_geom_ptr_t tcc,
 glutils::bufobj_ptr_t cellbo,
 glutils::bufobj_ptr_t nrmbo)
{
  tri_idx_list_t tlist;

  for(cellid_list_t::iterator it = mfold.begin(); it!= mfold.end(); ++it)
  {
    cellid_t st[40];

    tcc->get_cell_points(*it,st);
    tlist.push_back(mk_tri_idx(st[0],st[1],st[2]));
  }

  renderable_ptr_t ren
      (create_buffered_triangles_ren(cellbo,make_buf_obj(tlist),nrmbo));

  return ren;
}

renderable_ptr_t get_minima_renderer
(mfold_t &mfold,tri_cc_geom_ptr_t tcc,
 glutils::bufobj_ptr_t cellbo,
 glutils::bufobj_ptr_t nrmbo)
{
  tri_idx_list_t tlist;

  for(cellid_list_t::iterator it = mfold.begin(); it!= mfold.end(); ++it)
  {
    cellid_t st[40];

    uint st_ct = tcc->get_vert_star(*it,st);

    for(uint i = 1; i < st_ct; i++)
      tlist.push_back(mk_tri_idx(st[i-1],st[i],*it));
  }

  renderable_ptr_t ren
      (create_buffered_triangles_ren(cellbo,make_buf_obj(tlist),nrmbo));

  return ren;
}

template <eGDIR dir>
renderable_ptr_t get_saddle_renderer
(mfold_t &mfold,tri_cc_geom_ptr_t tcc,
 glutils::bufobj_ptr_t cellbo)
{

  map<cellid_t,int> pt_idx;

  line_idx_list_t e_idxs;

  for(cellid_list_t::iterator it = mfold.begin(); it!= mfold.end(); ++it)
  {
    cellid_t pt[20];

    tcc->get_cell_points(*it,pt);

    int npts = get_edge_pts<dir>(*it,pt,*(tcc));

    if( pt_idx.count(pt[0]) == 0) pt_idx[pt[0]] = pt_idx.size()-1;
    if( pt_idx.count(pt[1]) == 0) pt_idx[pt[1]] = pt_idx.size()-1;
    e_idxs.push_back(la::make_vec<uint>(pt_idx[pt[0]],pt_idx[pt[1]]));

    if(dir == DES || npts < 3) continue;


    if(pt_idx.count(pt[2]) == 0) pt_idx[pt[2]] = pt_idx.size()-1;
    e_idxs.push_back(la::make_vec<uint>(pt_idx[pt[1]],pt_idx[pt[2]]));
  }

  vertex_list_t pts(pt_idx.size());

  for(map<cellid_t,int>::iterator it = pt_idx.begin(); it!= pt_idx.end();++it)
    pts[it->second] = tcc->get_cell_position(it->first);

  smooth_lines(pts,e_idxs,4);

  renderable_ptr_t ren(create_buffered_lines_ren
                       (make_buf_obj(pts),make_buf_obj(e_idxs)));

  return ren;
}


void mscomplex_ren_t::update_geom()
{
  for( int i = 0 ; i < m_surv_cps.size(); ++i)
  {
    int j = m_surv_cps[i];

    switch(m_msc->index(j))
    {
    case 0:
      if(m_surv_mfold_show[1][i] && !m_surv_mfold_rens[1][i])
        m_surv_mfold_rens[1][i] = get_minima_renderer
            (m_msc->mfold<ASC>(j),m_tcc,m_cell_pos_bo,m_cell_nrm_bo);
      break;
    case 1:
      if(m_surv_mfold_show[0][i] && !m_surv_mfold_rens[0][i])
      {
        m_surv_mfold_rens[0][i] = get_saddle_renderer<DES>
          (m_msc->mfold<DES>(j),m_tcc,m_cell_pos_bo);
      }

      if(m_surv_mfold_show[1][i] && !m_surv_mfold_rens[1][i])
      {
        m_surv_mfold_rens[1][i] = get_saddle_renderer<ASC>
          (m_msc->mfold<ASC>(j),m_tcc,m_cell_pos_bo);
      }

      break;
    case 2:

      if(m_surv_mfold_show[0][i] && !m_surv_mfold_rens[0][i])
        m_surv_mfold_rens[0][i] = get_maxima_renderer
            (m_msc->mfold<DES>(j),m_tcc,m_cell_pos_bo,m_cell_nrm_bo);
      break;
    }

    if(!m_surv_mfold_show[0][i] && m_surv_mfold_rens[0][i])
      m_surv_mfold_rens[0][i].reset();

    if(!m_surv_mfold_show[1][i] && m_surv_mfold_rens[1][i])
      m_surv_mfold_rens[1][i].reset();
  }
}

}

#include <sstream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string_regex.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/static_assert.hpp>

#include <GL/glew.h>

#include <glutils.h>
#include <GLSLProgram.h>
#include <logutil.h>

#include <trimesh_viewer.h>
#include <trimesh_datamanager.h>
#include <trimesh_mscomplex.h>
#include <trimesh_mscomplex_ensure.h>
#include <trimesh_dataset.h>

template<typename T> std::string to_string(const T & t)
{
  std::stringstream ss;
  ss<<t;
  return ss.str();
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

  glutils::color_t g_disc_colors[GRADDIR_COUNT][gc_max_cell_dim+1] =
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
      (data_manager_t * gdm):
      m_scale_factor(0),
      m_bRebuildRens(true),m_bShowRoiBB(false),m_bCenterToRoi(false),
      m_gdm(gdm)
  {
    for(uint i = 0 ;i < m_gdm->m_pieces.size();++i)
      m_piece_rens.push_back(new octtree_piece_rendata(m_gdm->m_pieces.at(i)));
  }

  viewer_t::~viewer_t()
  {
    for ( uint i = 0 ; i < m_piece_rens.size();i++ )
      delete m_piece_rens[i];

    m_piece_rens.clear();

    delete m_gdm;
  }

  void viewer_t::set_roi_dim_range_nrm(double l,double u,int dim)
  {
    if(!(l<=u && 0.0 <= l && u <=1.0 && 0<=dim && dim < 3))
      return;

    double span = m_extent[dim].span();

    m_roi[dim][0]  = m_extent[dim][0] + (uint)(l*span);
    m_roi[dim][1]  = m_extent[dim][0] + (uint)(u*span);

    m_roi_base_pt  = ((m_roi.upper_corner() +  m_roi.lower_corner())/2);
  }

  int viewer_t::render()
  {
    if(m_bRebuildRens)
    {
      build_rens();

      m_bRebuildRens = false;
    }

    glPushAttrib(GL_ENABLE_BIT);

    glEnable(GL_RESCALE_NORMAL);

    glScalef(m_scale_factor,m_scale_factor,m_scale_factor);

    point_t s = ((m_extent.upper_corner() +  m_extent.lower_corner())/2);

    if(m_bCenterToRoi)
      glTranslatef(-m_roi_base_pt[0],-m_roi_base_pt[1],-m_roi_base_pt[2]);
    else
      glTranslated(-s[0],-s[1],-s[2]);

    if(m_bShowRoiBB)
    {
      glPushAttrib(GL_ENABLE_BIT);

      glDisable(GL_LIGHTING);

      glColor3dv(g_roiaabb_color.data());

      glutils::draw_aabb_line(m_roi.lower_corner(),m_roi.upper_corner());

      glPopAttrib();
    }

    for ( uint i = 0 ; i < m_piece_rens.size();i++ )
    {
      m_piece_rens[i]->render_msgraph_data();
    }

    for ( uint i = 0 ; i < m_piece_rens.size();i++ )
    {
      m_piece_rens[i]->render_dataset_data();
    }

    glPopAttrib();
  }

  void viewer_t::build_rens()
  {
    for ( uint i = 0 ; i < m_piece_rens.size();i++ )
    {
      m_piece_rens[i]->create_cp_rens(m_roi);
      m_piece_rens[i]->create_canc_cp_rens(m_roi);
      m_piece_rens[i]->create_grad_rens(m_roi);
    }
  }

  void viewer_t::init()
  {
    glutils::init();

    if(m_extent.eff_dim() == 0)
      throw std::runtime_error("NULL extent for viewer");

    glutils::tri_idx_list_t tlist;

    glutils::vertex_list_t  vlist;

    glutils::read_tri_file(m_gdm->m_tri_filename.c_str(),vlist,tlist);

    glutils::compute_extent(vlist,m_extent.data()->data());

    m_roi = m_extent;

    point_t s = m_extent.span();

    m_scale_factor =1.0/ *std::max_element(s.begin(),s.end());

    m_piece_rens[0]->tri_cc.reset(new tri_cc_t);
    m_piece_rens[0]->tri_cc->init(tlist,vlist.size());
    m_piece_rens[0]->create_cell_pos_nrm_bo(vlist);
    m_piece_rens[0]->create_disc_rds();

  }

  configurable_t::data_index_t viewer_t::dim()
  {
    return data_index_t(11,m_piece_rens.size());
  }

  bool viewer_t::exchange_field(const data_index_t &idx, boost::any &v)
  {
    octtree_piece_rendata * otprd = m_piece_rens[idx[1]];

    switch(idx[0])
    {
    case 0: return s_exchange_data_ro(otprd->dp->label(),v);
    case 1: return s_exchange_data_rw(otprd->m_bShowAllCps,v);
    case 2: return s_exchange_data_rw(otprd->m_bShowCps[0],v);
    case 3: return s_exchange_data_rw(otprd->m_bShowCps[1],v);
    case 4: return s_exchange_data_rw(otprd->m_bShowCps[2],v);
    case 5: return s_exchange_data_rw(otprd->m_bShowCpLabels,v);
    case 6: return s_exchange_data_rw(otprd->m_bShowMsGraph,v);
    case 7: return s_exchange_data_rw(otprd->m_bShowGrad,v);
    case 8: return s_exchange_data_rw(otprd->m_bShowCancCps,v);
    case 9: return s_exchange_data_rw(otprd->m_bShowCancMsGraph,v);
    case 10: return s_exchange_data_rw(otprd->m_bShowCellNormals,v);
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


  octtree_piece_rendata::octtree_piece_rendata (datapiece_t * _dp):
      m_bShowAllCps(false),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowGrad ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      m_bNeedUpdateDiscRens(false),
      m_bShowCellNormals(false),
      dp(_dp)
  {
    using namespace boost::lambda;

    std::for_each(m_bShowCps,m_bShowCps+gc_max_cell_dim+1,_1 = false);

  }

  void octtree_piece_rendata::create_cell_pos_nrm_bo(const glutils::vertex_list_t &vlist)
  {
    tri_cc_geom_t tcc_geom;

    tcc_geom.init(tri_cc,vlist);

    cell_pos_bo = glutils::make_buf_obj(tcc_geom.get_cell_positions());

    cell_nrm_bo = glutils::make_buf_obj(tcc_geom.get_cell_normals());

    ren_cell_normals.reset(glutils::create_buffered_normals_ren
                           (cell_pos_bo,
                            glutils::make_buf_obj(),
                            glutils::make_buf_obj(),
                            cell_nrm_bo,
                            tcc_geom.get_average_edge_length()/2));
  }

  void  octtree_piece_rendata::create_cp_rens(const rect_t & roi)
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<glutils::point_idx_t>   crit_pt_idxs[gc_max_cell_dim+1];
    std::vector<glutils::line_idx_t>    crit_conn_idxs[gc_max_cell_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      uint index = dp->msgraph->m_cps[i]->index;

      crit_pt_idxs[index].push_back(c);
    }

    for(uint i = 0 ; i < gc_max_cell_dim+1; ++i)
    {
      ren_cp[i].reset(glutils::create_buffered_points_ren
                      (cell_pos_bo,
                       glutils::make_buf_obj(crit_pt_idxs[i]),
                       glutils::make_buf_obj()));
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      uint index = dp->msgraph->m_cps[i]->index;

      for(conn_iter_t it  = dp->msgraph->m_cps[i]->conn[0].begin();
      it != dp->msgraph->m_cps[i]->conn[0].end(); ++it)
      {
        cellid_t c_ = dp->msgraph->m_cps[*it]->cellid;

        crit_conn_idxs[index-1].push_back(glutils::line_idx_t(c,c_));
      }
    }

    for(uint i = 0 ; i < gc_max_cell_dim; ++i)
    {
      ren_cp_conns[i].reset(glutils::create_buffered_lines_ren
                            (cell_pos_bo,
                             glutils::make_buf_obj(crit_conn_idxs[i]),
                             glutils::make_buf_obj()));
    }

  }

  void  octtree_piece_rendata::create_canc_cp_rens(const rect_t & roi)
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<glutils::point_idx_t>   canc_cp_idxs[gc_max_cell_dim+1];
    std::vector<glutils::line_idx_t>    canc_cp_conn_idxs[gc_max_cell_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(!dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      uint index = dp->msgraph->m_cps[i]->index;

      canc_cp_idxs[index].push_back(c);

    }

    for(uint i = 0 ; i < gc_max_cell_dim+1; ++i)
    {
      ren_canc_cp[i].reset(glutils::create_buffered_points_ren
                           (cell_pos_bo,
                            glutils::make_buf_obj(canc_cp_idxs[i]),
                            glutils::make_buf_obj()));
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(!dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      uint index = dp->msgraph->m_cps[i]->index;

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        for(conn_iter_t it  = dp->msgraph->m_cps[i]->conn[dir].begin();
        it != dp->msgraph->m_cps[i]->conn[dir].end(); ++it)
        {
          cellid_t c_ = dp->msgraph->m_cps[*it]->cellid;

          canc_cp_conn_idxs[index-(dir^1)].push_back
              (glutils::line_idx_t(c,c_));
        }
      }
    }

    for(uint i = 0 ; i < gc_max_cell_dim; ++i)
    {
      ren_canc_cp_conns[i].reset(glutils::create_buffered_lines_ren
                                 (cell_pos_bo,
                                  glutils::make_buf_obj(canc_cp_conn_idxs[i]),
                                  glutils::make_buf_obj()));
    }

  }

  void octtree_piece_rendata::create_grad_rens(const rect_t & roi)
  {
#warning "create_grad_rens not implemented"
  }

  void octtree_piece_rendata::create_disc_rds()
  {
    if(dp->msgraph == NULL)
      return;

    boost::shared_ptr<disc_rendata_t> sptr;

    for(uint i = 0 ; i < dp->msgraph->m_cps.size();++i)
    {
      critpt_t * cp = dp->msgraph->m_cps[i];

      if(cp->is_paired) continue;

      sptr.reset(new disc_rendata_t(cp->cellid,cp->index));

      disc_rds.push_back(sptr);
    }
  }

  void octtree_piece_rendata::update_active_disc_rens()
  {
    if(dp->msgraph == NULL)
      return;

    for(uint i = 0 ; i < disc_rds.size();++i)
    {
      if(disc_rds[i]->update(this))
      {
        active_disc_rens[disc_rds[i]->index].insert(disc_rds[i]);
      }
      else
      {
        active_disc_rens[disc_rds[i]->index].erase(disc_rds[i]);
      }
    }
  }

  void octtree_piece_rendata::render_msgraph_data()
  {
    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_LIGHTING );

    glPointSize ( 4.0 );

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
  }

  void octtree_piece_rendata::render_dataset_data()
  {
    if(m_bNeedUpdateDiscRens)
    {
      update_active_disc_rens();
      m_bNeedUpdateDiscRens = false;
    }

    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_LIGHTING );

    if(m_bShowGrad)
    {
      for(uint i = 0 ; i < gc_max_cell_dim; ++i)
      {
        if(ren_grad[i])
        {
          glColor3dv ( g_grad_colors[i].data() );

          ren_grad[i]->render();
        }
      }
    }

    for(disc_rendata_sp_set_t::iterator it = active_disc_rens[1].begin();
        it != active_disc_rens[1].end() ; ++it)
    {
      (*it)->render();
    }

    glEnable ( GL_LIGHTING );

    for(disc_rendata_sp_set_t::iterator it = active_disc_rens[0].begin();
        it != active_disc_rens[0].end() ; ++it)
    {
      (*it)->render();
    }

    for(disc_rendata_sp_set_t::iterator it = active_disc_rens[2].begin();
        it != active_disc_rens[2].end() ; ++it)
    {
      (*it)->render();
    }

    glColor3dv(g_normals_color.data());

    if(m_bShowCellNormals && ren_cell_normals)
      ren_cell_normals->render();

    glPopAttrib();
    glPopMatrix();

  }

  struct random_color_assigner
  {
    disc_rendata_sp_t m_drd;

    int m_no;

    static const uint MAX_RAND = 256;

    random_color_assigner(disc_rendata_sp_t drd,int no):m_drd(drd),m_no(no){};

    void operator()()
    {
      for(uint c = 0 ; c < 3 ; ++c)
        m_drd->color[m_no][c] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);

    }
  };

  configurable_t::data_index_t octtree_piece_rendata::dim()
  {
    return data_index_t(8,disc_rds.size());
  }

  bool octtree_piece_rendata::exchange_field
      (const data_index_t &idx,boost::any &v)
  {
    bool is_read     = v.empty();

    int i = idx[0];

    disc_rendata_sp_t drd = disc_rds[idx[1]];

    switch(i)
    {
    case 0:
      return s_exchange_data_ro((int)drd->cellid,v);
    case 1:
      return s_exchange_data_ro((int)drd->index,v);
    case 2:
    case 3:
      {
        bool need_update =  s_exchange_data_rw(drd->show[i%2],v);

        if(need_update && is_read == false )
          m_bNeedUpdateDiscRens = true;

        return need_update;
      }
    case 4:
    case 5:
      return s_exchange_data_rw(drd->color[i%2],v);
    case 6:
    case 7:
      return s_exchange_action(random_color_assigner(drd,i%2),v);
    };

     throw std::logic_error("octtree_piece_rendata::invalid index");
  }

  configurable_t::eFieldType octtree_piece_rendata::exchange_header
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
    }
    throw std::logic_error("octtree_piece_rendata::invalid index");
  }

  disc_rendata_t::disc_rendata_t(cellid_t c,uint i):cellid(c),index(i)
  {
    color[0] = g_disc_colors[1][index];
    color[1] = g_disc_colors[0][index];

    show[0] =false; ren[0] =NULL;
    show[1] =false; ren[1] =NULL;
  }

  disc_rendata_t::~disc_rendata_t()
  {
    show[0] =false;
    show[1] =false;

    update(NULL);

  }

  void disc_rendata_t::render()
  {
    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] && ren[dir])
      {
        glColor3dv(color[dir].data());
        ren[dir]->render();
      }
    }
  }

  bool disc_rendata_t::update(octtree_piece_rendata *drd)
  {
    using namespace boost::lambda;

    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] && this->ren[dir] == NULL && drd->dp->msgraph)
      {
        ensure_cellid_critical(drd->dp->msgraph,cellid);

        critpt_t *cp = drd->dp->msgraph->m_cps[drd->dp->msgraph->m_id_cp_map[cellid]];

        std::set<cellid_t> cset;

        for(uint j = 0 ; j < cp->contrib[dir].size();++j)
        {
          critpt_t *cp_contrib = drd->dp->msgraph->m_cps[cp->contrib[dir][j]];

          if(cp_contrib->index != cp->index)
            throw std::logic_error("contrib and cp must have same idx");

          for(uint i = 0; i < cp_contrib->disc[dir].size(); ++i)
          {
            cellid_t c = cp_contrib->disc[dir][i];

            if(cset.count(c) == 0)
              cset.insert(c);
          }
        }

        if(cp->index == 1 && dir == 0)
        {
          glutils::line_idx_list_t e_idxs;

          for(std::set<cellid_t>::iterator it = cset.begin(); it!= cset.end(); ++it)
          {
            cellid_t pt[20];

            drd->tri_cc->get_cell_points(*it,pt);

            e_idxs.push_back(glutils::line_idx_t(pt[0],pt[1]));

          }

          ren[dir] = glutils::create_buffered_lines_ren
                     (drd->cell_pos_bo,
                      glutils::make_buf_obj(e_idxs),
                      glutils::make_buf_obj());

        }

        if(cp->index == 1 && dir == 1)
        {
          glutils::line_idx_list_t e_idxs;

          for(std::set<cellid_t>::iterator it = cset.begin(); it!= cset.end(); ++it)
          {
            cellid_t cf[20];

            uint cf_ct = drd->tri_cc->get_cell_co_facets(*it,cf);

            e_idxs.push_back(glutils::line_idx_t(*it,cf[0]));

            if(cf_ct == 2)
              e_idxs.push_back(glutils::line_idx_t(*it,cf[1]));

          }

          ren[dir] = glutils::create_buffered_lines_ren
                     (drd->cell_pos_bo,
                      glutils::make_buf_obj(e_idxs),
                      glutils::make_buf_obj());

        }

        if(cp->index == 2 && dir == 0)
        {
          glutils::tri_idx_list_t t_idxs;

          for(std::set<cellid_t>::iterator it = cset.begin(); it!= cset.end(); ++it)
          {
            cellid_t pt[20];

            drd->tri_cc->get_cell_points(*it,pt);

            t_idxs.push_back(glutils::tri_idx_t(pt[0],pt[2],pt[1]));

          }

          ren[dir] = glutils::create_buffered_triangles_ren
                     (drd->cell_pos_bo,
                      glutils::make_buf_obj(t_idxs),
                      glutils::make_buf_obj(),
                      drd->cell_nrm_bo);
        }

        if(cp->index == 0 && dir == 1)
        {
          glutils::tri_idx_list_t t_idxs;

          for(std::set<cellid_t>::iterator it = cset.begin(); it!= cset.end(); ++it)
          {
            cellid_t st[40];

            uint st_ct = drd->tri_cc->get_vert_star(*it,st);

            for(uint i = 1; i < st_ct; i++)
            {
              t_idxs.push_back(glutils::tri_idx_t(st[i-1],*it,st[i]));
            }

            if(st_ct%2 == 0)
              t_idxs.push_back(glutils::tri_idx_t(st[st_ct-1],*it,st[0]));
          }

          ren[dir] = glutils::create_buffered_triangles_ren
                     (drd->cell_pos_bo,
                      glutils::make_buf_obj(t_idxs),
                      glutils::make_buf_obj(),
                      drd->cell_nrm_bo);

        }

      }

      if(!show[dir] && this->ren[dir] != NULL )
      {
        delete ren[dir];

        ren[dir] = NULL;

      }
    }

    return (show[0] || show[1]);
  }

}

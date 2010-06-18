#include <QMenu>
#include <QTreeView>
#include <QColorDialog>
#include <QDebug>
#include <QKeyEvent>

#include <boost/typeof/typeof.hpp>

#include <trimesh_viewer.h>
#include <trimesh_viewer_mainwindow.h>
#include <trimesh_datamanager.h>

const int g_roi_show_aabb_time_msec = 5*1000;

namespace trimesh
{

  void glviewer_t::draw()
  {
    m_ren->render();
  }

  void glviewer_t::init()
  {
    setSnapshotFormat("PNG");

    setSnapshotQuality(100);

    glEnable ( GL_CULL_FACE );

    glCullFace ( GL_BACK );

    glPolygonMode ( GL_FRONT, GL_FILL );

    glPolygonMode ( GL_BACK, GL_LINE );

    m_ren->init();
  }

  glviewer_t::glviewer_t(QWidget * par):
      m_is_recording(false),
      m_bf_cull(true),
      m_wireframe(false),
      m_ren(NULL)
  {
    setParent(par);
  }

  void glviewer_t::setup(data_manager_t *dm)
  {
    m_ren = new viewer_t(dm);
  }

  glviewer_t::~glviewer_t()
  {
    if(m_ren != NULL)
      delete m_ren;
  }

  void glviewer_t::keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==Qt::Key_C) && (modifiers==Qt::ControlModifier))
    {
      m_is_recording = !m_is_recording;

      if(m_is_recording)
        connect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));
      else
        disconnect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));
    }
    else if ((e->key()==Qt::Key_R) && (modifiers==Qt::ControlModifier))
    {
      m_ren->m_bShowRoiBB = !m_ren->m_bShowRoiBB;
    }
    else if ((e->key()==Qt::Key_B) && (modifiers==Qt::ControlModifier))
    {
      m_bf_cull = !m_bf_cull;

      if(m_bf_cull)
        glEnable ( GL_CULL_FACE );
      else
        glDisable ( GL_CULL_FACE );
    }
    else if ((e->key()==Qt::Key_W) && (modifiers==Qt::ControlModifier))
    {
      m_wireframe = !m_wireframe;

      if(m_wireframe)
        glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );
      else
      {
        glPolygonMode ( GL_FRONT, GL_FILL );
        glPolygonMode ( GL_BACK, GL_LINE );

      }
    }
    else
    {
      QGLViewer::keyPressEvent(e);
    }

    updateGL();
  }

  QString glviewer_t::helpString() const
  {
    QString text("<h2>MS Complex Viewer</h2>");
    return text;
  }

  void viewer_mainwindow::on_datapiece_view_customContextMenuRequested  ( const QPoint &p )
  {
    QModelIndexList l =  datapiece_view->selectionModel()->selectedIndexes();

    configurable_ctx_menu(glviewer->m_ren,l,datapiece_view->mapToGlobal(p));

    glviewer->updateGL();
  }

  void viewer_mainwindow::on_critpt_view_customContextMenuRequested ( const QPoint &p )
  {

    QModelIndexList l = m_cp_model_proxy->mapSelectionToSource
                        (critpt_view->selectionModel()->selection()).indexes();

    configurable_ctx_menu(glviewer->m_ren->m_piece_rens[m_active_otp_idx],
                          l,critpt_view->mapToGlobal(p));

    glviewer->updateGL();
  }

  void viewer_mainwindow::on_datapiece_view_activated ( const QModelIndex & index  )
  {
    if(m_active_otp_idx == index.row())
      return;

    m_active_otp_idx = index.row();

    m_cp_model->reset_configurable
        (glviewer->m_ren->m_piece_rens[m_active_otp_idx]);
  }

  inline double get_nrm_value(double d_val,double d_min,double d_max)
  {
    return (d_val - d_min)/(d_max - d_min);
  }

  void viewer_mainwindow::update_roi_box(double l,double u,uint dim)
  {
    glviewer->m_ren->set_roi_dim_range_nrm(l,u,dim);

    if (m_clear_roi_aabb_timer->isActive() ||
        glviewer->m_ren->m_bShowRoiBB == false)
    {
      glviewer->m_ren->m_bShowRoiBB = true;

      m_clear_roi_aabb_timer->start();
    }

    glviewer->updateGL();
  }

  void viewer_mainwindow::on_xroi_spanslider_spanChanged(int l , int u )
  {
    BOOST_AUTO(sldr,xroi_spanslider);

    update_roi_box(get_nrm_value(l,sldr->minimum(),sldr->maximum()),
                   get_nrm_value(u,sldr->minimum(),sldr->maximum()),0);

  }

  void viewer_mainwindow::on_yroi_spanslider_spanChanged(int l , int u )
  {
    BOOST_AUTO(sldr,yroi_spanslider);

    update_roi_box(get_nrm_value(l,sldr->minimum(),sldr->maximum()),
                   get_nrm_value(u,sldr->minimum(),sldr->maximum()),1);

  }

  void viewer_mainwindow::on_zroi_spanslider_spanChanged(int l , int u )
  {
    BOOST_AUTO(sldr,zroi_spanslider);

    update_roi_box(get_nrm_value(l,sldr->minimum(),sldr->maximum()),
                   get_nrm_value(u,sldr->minimum(),sldr->maximum()),2);

  }

  void viewer_mainwindow::on_update_roi_pushButton_clicked(bool)
  {
    glviewer->m_ren->m_bRebuildRens = true;

    glviewer->updateGL();
  }

  void viewer_mainwindow::on_center_to_roi_checkBox_clicked(bool state)
  {
    glviewer->m_ren->m_bCenterToRoi = state;

    glviewer->updateGL();
  }

  void viewer_mainwindow::clear_roi_aabb()
  {
    glviewer->m_ren->m_bShowRoiBB = false;

    glviewer->updateGL();
  }

  viewer_mainwindow::viewer_mainwindow
      (data_manager_t * gdm):
      m_active_otp_idx(0)
  {
    setupUi (this);

    glviewer->setup(gdm);

    m_otp_model = new configurable_item_model
                  (glviewer->m_ren,this);

    spinviewer->setScene(new QGraphicsScene(this));

    spinviewer->scene()->addItem(new spin::si_graphics_item_t(glviewer->m_ren));

    datapiece_view->setModel ( m_otp_model );

    m_cp_model = new configurable_item_model
                 (glviewer->m_ren->m_piece_rens[m_active_otp_idx],this);

    m_cp_model_proxy = new QSortFilterProxyModel(this);

    m_cp_model_proxy->setSourceModel(m_cp_model);

    critpt_view->setModel ( m_cp_model_proxy );

    connect(critpt_filter_edit,SIGNAL(textChanged(QString)),
            m_cp_model_proxy,SLOT(setFilterFixedString(QString)));

    m_clear_roi_aabb_timer = new QTimer(this);

    m_clear_roi_aabb_timer->setSingleShot(true);

    m_clear_roi_aabb_timer->setInterval(g_roi_show_aabb_time_msec);

    connect(m_clear_roi_aabb_timer,SIGNAL(timeout()),this,SLOT(clear_roi_aabb()));
  }

  void viewer_mainwindow::showEvent ( QShowEvent * )
  {
    m_cp_model->force_reset();
  }

  viewer_mainwindow::~viewer_mainwindow()
  {
  }

  inline QColor to_qcolor (const glutils::color_t & c)
  {
    return QColor::fromRgbF(c[0],c[1],c[2]);
  }

  QVariant configurable_item_model::data
      ( const QModelIndex &index, int role ) const
  {
    if ( !index.isValid() )
      return QVariant();

    if(index.column() >= m_column_idxs.size())
      return QVariant();

    configurable_t::data_index_t idx(m_column_idxs[index.column()],index.row());

    boost::any val;

    m_conf->exchange_field(idx,val);

    if(role == Qt::DisplayRole)
    {
      if(val.type() == typeid(std::string))
        return QString(boost::any_cast<std::string>(val).c_str());
      else if (val.type() == typeid(bool))
        return boost::any_cast<bool>(val);
      else if (val.type() == typeid(int))
        return boost::any_cast<int>(val);
      else if (val.type() == typeid(glutils::color_t))
        return to_qcolor(boost::any_cast<glutils::color_t>(val));
    }
    else if( role == Qt::DecorationRole)
    {
      if (val.type() == typeid(glutils::color_t))
        return to_qcolor(boost::any_cast<glutils::color_t>(val));
    }

    return QVariant();
  }

  int configurable_item_model::columnCount ( const QModelIndex &parent  ) const
  {
    return m_column_idxs.size();
  }

  QVariant configurable_item_model::headerData
      ( int section, Qt::Orientation orientation,int role ) const
  {
    if ( orientation == Qt::Horizontal &&
         role == Qt::DisplayRole &&
         section < m_column_idxs.size())
    {
      boost::any h;

      m_conf->exchange_header(m_column_idxs[section],h);

      if(h.type() == typeid(std::string))
        return (boost::any_cast<std::string>(h)).c_str();
    }
    return QVariant();
  }

  int configurable_item_model::rowCount ( const QModelIndex &parent ) const
  {
    return m_conf->dim()[1];
  }

  void configurable_item_model::reset_configurable(configurable_t *conf)
  {
    if(conf == m_conf)
      return;
    m_conf = conf;

    force_reset();
  }

  void configurable_item_model::setColumnFilter(int columnFilter)
  {
    m_column_filter = columnFilter;

    m_column_idxs.clear();

    configurable_t::data_index_t idx = m_conf->dim();

    boost::any val;

    for(uint i = 0; i < idx[0] ; ++i)
    {

      bool col_usable = false;

      switch(m_conf->exchange_header(i,val))
      {
      case configurable_t::EFT_DATA_RO:
        col_usable = m_column_filter&CF_EFT_DATA_RO;break;
      case configurable_t::EFT_DATA_RW:
        col_usable = m_column_filter&CF_EFT_DATA_RW;break;
      case configurable_t::EFT_ACTION:
        col_usable = m_column_filter&CF_EFT_ACTION;break;
      };

      if(col_usable)
        m_column_idxs.push_back(i);
    }
  }

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p)
  {

    if(l.size() == 0)
      return;

    QMenu m;

    std::set<int> row_set;

    for(uint i = 0 ; i < l.size(); ++i )
      row_set.insert(l[i].row());

    std::vector<int> rows(row_set.size());

    std::copy(row_set.begin(),row_set.end(),rows.begin());

    for(uint i = 0 ; i < c->dim()[0];++i)
    {
      boost::any h;

      if(c->exchange_header(i,h) == configurable_t::EFT_DATA_RO)
        continue;

      std::string hdr("could not read column header");

      if(h.type() == typeid(std::string))
        hdr = (boost::any_cast<std::string>(h)).c_str();

      QAction * action  = m.addAction ( hdr.c_str());

      std::vector<boost::any> vals;

      for(uint j = 0 ; j < rows.size();++j )
      {
        boost::any val;

        configurable_t::data_index_t idx(i,rows[j]);

        c->exchange_field(idx,val);

        vals.push_back(val);
      }

      if(vals[0].type() == typeid(bool))
      {
        action->setCheckable(true);
        action->setChecked(boost::any_cast<bool>(vals[0]));
      }

      configurable_ctx_menu_sig_collector * coll =
          new configurable_ctx_menu_sig_collector(c,vals,i,rows,&m);

      m.connect(action,SIGNAL ( triggered ( bool ) ),
                    coll,SLOT(triggered ( bool )));
    }

    m.exec(p);
  }

  void configurable_ctx_menu_sig_collector::triggered(bool state)
  {
    boost::any out_val;

    if(m_vals[0].type() == typeid(bool))
      out_val = boost::any(state);

    if(m_vals[0].type() == typeid(glutils::color_t))
    {
      glutils::color_t c = boost::any_cast<glutils::color_t>(m_vals[0]);

      QColor ic = QColor::fromRgbF(c[0],c[1],c[2],1.0);

      QColor qc = QColorDialog::getColor(ic);

      if(qc.isValid())
        out_val = glutils::color_t(qc.redF(),qc.greenF(),qc.blueF());
    }

    if(m_vals[0].type() == typeid(configurable_t::action_callback_t))
    {
      for(uint i = 0 ; i < m_vals.size();++i)
        boost::any_cast<configurable_t::action_callback_t>(m_vals[i])();
    }

    if(out_val.empty())
      return;

    for(uint i = 0 ; i < m_rows.size();++i)
    {
      configurable_t::data_index_t idx;
      idx[1] = m_rows[i];
      idx[0] = m_col;

      m_conf->exchange_field(idx,out_val);
    }
  }

}

namespace spin
{
  QRectF si_graphics_item_t::boundingRect() const
  {
    QRectF r;

    if(!m_si) return r;

    si_point_t si_lc = m_si->m_iextent.lower_corner();
    si_point_t si_uc = m_si->m_iextent.upper_corner();

    r.setBottomLeft(to_qpointf(si_lc));
    r.setTopRight(to_qpointf(si_uc));

    return r;
  }


  QRgb to_qrgb(si_scalar_t c)
  {
    c*= 20;

    return qRgb(c,c,c);
  }

  void si_graphics_item_t::update_image()
  {
    if(m_si == m_viewer->m_spin_image) return;

    m_si = m_viewer->m_spin_image;

    if(!m_si) return;

    si_point_t sz = m_si->m_iextent.span();

    si_point_t lc = m_si->m_iextent.lower_corner();

    m_image = new QImage(QSize(sz[0],sz[1]),QImage::Format_RGB888);

    for(uint y = 0 ; y < sz[1]; ++y)
      for(uint x = 0 ; x < sz[0]; ++x)
        m_image->setPixel( x,y,to_qrgb((*m_si->m_image)(lc + si_ipoint_t(x,y))));

  }

  void si_graphics_item_t::paint
      (QPainter *painter,
       const QStyleOptionGraphicsItem *option,
       QWidget *widget)
  {
    update_image();

    if(!m_si) return;

    si_point_t lc = m_si->m_iextent.lower_corner();

    painter->drawImage(QPoint(lc[0],lc[1]),*m_image);
  }
}

#include <QMenu>
#include <QTreeView>
#include <QColorDialog>
#include <QDebug>
#include <QKeyEvent>

#include <boost/typeof/typeof.hpp>

#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>

const int g_roi_show_aabb_time_msec = 5*1000;

namespace grid
{

  void glviewer_t::draw()
  {
    m_ren->render();
  }

  void glviewer_t::init()
  {
    setSnapshotFormat("PNG");
    setSnapshotQuality(100);

    m_ren->init();
  }

  glviewer_t::glviewer_t(data_manager_t * gdm):m_is_recording(false)
  {
    m_ren = new grid_viewer_t(gdm);
  }

  glviewer_t::~glviewer_t()
  {
    delete m_ren;
  }

  void glviewer_t::keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();

    bool handled = false;
    if ((e->key()==Qt::Key_R) && (modifiers==Qt::ControlModifier))
    {
      m_is_recording = !m_is_recording;

      if(m_is_recording)
        connect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));
      else
        disconnect(this, SIGNAL(drawFinished(bool)),this, SLOT(saveSnapshot(bool)));

      std::cout<<"Recording state::"<<m_is_recording<<"\n";

      handled = true;
    }
    else if ((e->key()==Qt::Key_B) && (modifiers==Qt::ControlModifier))
    {
      m_ren->m_bShowRoiBB = !m_ren->m_bShowRoiBB;
      handled = true;
    }

    if (handled)
      updateGL();
    else
      QGLViewer::keyPressEvent(e);

  }

  QString glviewer_t::helpString() const
  {
    QString text("<h2>MS Complex Viewer</h2>");
    return text;
  }

  void viewer_mainwindow::on_datapiece_view_customContextMenuRequested  ( const QPoint &p )
  {
    QModelIndexList l =  datapiece_view->selectionModel()->selectedIndexes();

    configurable_ctx_menu(m_viewer->m_ren,l,datapiece_view->mapToGlobal(p));

    m_viewer->updateGL();
  }

  void viewer_mainwindow::on_critpt_view_customContextMenuRequested ( const QPoint &p )
  {

    QModelIndexList l = m_cp_model_proxy->mapSelectionToSource
                        (critpt_view->selectionModel()->selection()).indexes();

    configurable_ctx_menu(m_viewer->m_ren->m_grid_piece_rens[m_active_otp_idx],
                          l,critpt_view->mapToGlobal(p));

    m_viewer->updateGL();
  }

  void viewer_mainwindow::on_datapiece_view_activated ( const QModelIndex & index  )
  {
    if(m_active_otp_idx == index.row())
      return;

    m_active_otp_idx = index.row();

    m_cp_model->reset_configurable
        (m_viewer->m_ren->m_grid_piece_rens[m_active_otp_idx]);
  }

  inline double get_nrm_value(double d_val,double d_min,double d_max)
  {
    return (d_val - d_min)/(d_max - d_min);
  }

  void viewer_mainwindow::update_roi_box(double l,double u,uint dim)
  {
    m_viewer->m_ren->set_roi_dim_range_nrm(l,u,dim);

    if (m_clear_roi_aabb_timer->isActive() ||
        m_viewer->m_ren->m_bShowRoiBB == false)
    {
      m_viewer->m_ren->m_bShowRoiBB = true;

      m_clear_roi_aabb_timer->start();
    }

    m_viewer->updateGL();
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
    m_viewer->m_ren->m_bRebuildRens = true;

    m_viewer->updateGL();
  }

  void viewer_mainwindow::on_center_to_roi_checkBox_clicked(bool state)
  {
    m_viewer->m_ren->m_bCenterToRoi = state;

    m_viewer->updateGL();
  }

  void viewer_mainwindow::clear_roi_aabb()
  {
    m_viewer->m_ren->m_bShowRoiBB = false;

    m_viewer->updateGL();
  }

  viewer_mainwindow::viewer_mainwindow
      (data_manager_t * gdm):
      m_active_otp_idx(0)
  {
    setupUi (this);

    m_viewer = new glviewer_t(gdm);

    m_viewer->setParent(glviewer);

    m_viewer->resize(glviewer->size());

    m_otp_model = new configurable_item_model (m_viewer->m_ren,this);

    datapiece_view->setModel ( m_otp_model );

    m_cp_model = new configurable_item_model
                 (m_viewer->m_ren->m_grid_piece_rens[m_active_otp_idx],this);

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

    if(index.column() >= m_conf->columns())
      return QVariant();

    configurable_t::data_index_t idx(index.column(),index.row());

    boost::any val;

    m_conf->exchange_data(idx,val);

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
    return m_conf->columns();
  }

  QVariant configurable_item_model::headerData
      ( int section, Qt::Orientation orientation,int role ) const
  {
    if ( orientation == Qt::Horizontal &&
         role == Qt::DisplayRole &&
         section < m_conf->columns())
    {

      boost::any h  = m_conf->get_header(section);

      if(h.type() == typeid(std::string))
        return (boost::any_cast<std::string>(h)).c_str();
    }
    return QVariant();
  }

  int configurable_item_model::rowCount ( const QModelIndex &parent ) const
  {
    return m_conf->rows();
  }

  void configurable_item_model::reset_configurable(configurable_t *conf)
  {
    if(conf == m_conf)
      return;
    m_conf = conf;

    reset();
  }

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p)
  {

    if(l.size() == 0)
      return;

    uint first_row = l[0].row();

    QMenu m;

    for(uint i = 0 ; i < c->columns();++i)
    {
      boost::any val;

      configurable_t::data_index_t idx(i,first_row);

      bool is_rw = c->exchange_data(idx,val);

      if(is_rw == false) continue;

      boost::any h  = c->get_header(i);

      std::string hdr("could not read column header");

      if(h.type() == typeid(std::string))
        hdr = (boost::any_cast<std::string>(h)).c_str();

      QAction * action  = m.addAction ( hdr.c_str());

      if(val.type() == typeid(bool))
      {
        action->setCheckable(true);
        action->setChecked(boost::any_cast<bool>(val));
      }

      configurable_ctx_menu_sig_collector * coll =
          new configurable_ctx_menu_sig_collector(c,val,i,l,&m);

      m.connect(action,SIGNAL ( triggered ( bool ) ),
                    coll,SLOT(triggered ( bool )));
    }

    m.exec(p);
  }

  void configurable_ctx_menu_sig_collector::triggered(bool state)
  {
    boost::any out_val;

    if(m_val.type() == typeid(bool))
      out_val = boost::any(state);

    if(m_val.type() == typeid(glutils::color_t))
    {

      glutils::color_t c = boost::any_cast<glutils::color_t>(m_val);

      QColor ic = QColor::fromRgbF(c[0],c[1],c[2],1.0);

      QColor qc = QColorDialog::getColor(ic);

      if(qc.isValid())
        out_val = glutils::color_t(qc.redF(),qc.greenF(),qc.blueF());
    }

    if(out_val.empty())
      return;

    for(uint i = 0 ; i < m_rows.size();++i)
    {
      configurable_t::data_index_t idx;
      idx[1] = m_rows[i].row();
      idx[0] = m_col;

      m_conf->exchange_data(idx,out_val);
    }
  }
}

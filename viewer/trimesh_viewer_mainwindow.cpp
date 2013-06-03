#include <QMenu>
#include <QTreeView>
#include <QColorDialog>
#include <QDebug>
#include <QKeyEvent>
#include <QFileDialog>
#include <QMessageBox>
#include <QDir>

#include <boost/typeof/typeof.hpp>
#include <boost/range/algorithm.hpp>

#include <trimesh_viewer.h>
#include <trimesh_mscomplex.h>
#include <trimesh_viewer_mainwindow.h>

const int g_roi_show_aabb_time_msec = 5*1000;

extern double g_max_cp_size;

using namespace std;

namespace trimesh
{

  void glviewer_t::draw()
  {
    m_ren->render();
  }

  bool glviewer_t::saveImageSnapshot(const QString& fileName,int over_sample)
  {
    double oversampling = over_sample;
    QSize subSize(int(this->width()/oversampling),
                  int(this->height()/oversampling));

    QSize finalSize = this->size();

    double aspectRatio = width() / static_cast<double>(height());

    double zNear = camera()->zNear();
    double zFar = camera()->zFar();

    double xMin, yMin;
    if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
    {
      xMin = zNear * tan(camera()->fieldOfView() / 2.0) * aspectRatio;
      yMin = xMin / aspectRatio;
    }
    else
    {
      camera()->getOrthoWidthHeight(xMin, yMin);
      yMin = xMin / aspectRatio;
    }

#if QT_VERSION >= 0x040000
    QImage image(finalSize.width(), finalSize.height(), QImage::Format_ARGB32);
#else
    QImage image(finalSize.width(), finalSize.height(), 32);
#endif

    if (image.isNull())
    {
      QMessageBox::warning(this, "Image saving error",
                           "Unable to create resulting image",
                           QMessageBox::Ok, QMessageBox::NoButton);
      return false;
    }

    double deltaX = 2.0 * xMin * subSize.width() / finalSize.width();
    double deltaY = 2.0 * yMin * subSize.height() / finalSize.height();

    int nbX = finalSize.width() / subSize.width();
    int nbY = finalSize.height() / subSize.height();

    // Extra subimage on the border if needed
    if (nbX * subSize.width() < finalSize.width())
      nbX++;
    if (nbY * subSize.height() < finalSize.height())
      nbY++;

    makeCurrent();

    int count=0;
    for (int i=0; i<nbX; i++)
      for (int j=0; j<nbY; j++)
      {
        preDraw();
        // Change projection matrix
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        if (camera()->type() == qglviewer::Camera::PERSPECTIVE)
          glFrustum(-xMin + i*deltaX, -xMin + (i+1)*deltaX,
                    yMin - (j+1)*deltaY, yMin - j*deltaY,
                    zNear, zFar);
        else
          glOrtho(-xMin + i*deltaX, -xMin + (i+1)*deltaX,
                  yMin - (j+1)*deltaY, yMin - j*deltaY,
                  zNear, zFar);
        glMatrixMode(GL_MODELVIEW);

        draw();
        postDraw();


        QImage snapshot = grabFrameBuffer(true);


#if QT_VERSION >= 0x040000
        QImage subImage = snapshot.scaled(subSize, Qt::IgnoreAspectRatio,
                                          Qt::SmoothTransformation);
#else
# if QT_VERSION >= 0x030000
        QImage subImage = snapshot.scale(subSize, QImage::ScaleFree);
# else
        QImage subImage = snapshot.smoothScale(subSize.width(),
                                               subSize.height());
# endif
#endif

        // Copy subImage in image
        for (int ii=0; ii<subSize.width(); ii++)
        {
          int fi = i*subSize.width() + ii;
          if (fi == image.width())
            break;
          for (int jj=0; jj<subSize.height(); jj++)
          {
            int fj = j*subSize.height() + jj;
            if (fj == image.height())
              break;
            image.setPixel(fi, fj, subImage.pixel(ii,jj));
          }
        }
        count++;
      }

#if QT_VERSION >= 0x040000
    bool saveOK = image.save(fileName, snapshotFormat().toLatin1().constData(),
                             snapshotQuality());
#else
    bool saveOK = image.save(fileName, snapshotFormat(), snapshotQuality());
#endif


    return saveOK;
  }

  void glviewer_t::init()
  {
    restoreStateFromFile();

    setSnapshotFormat("PNG");

    setSnapshotQuality(100);

    glEnable ( GL_CULL_FACE );

    glCullFace ( GL_BACK );

    glPolygonMode ( GL_FRONT, GL_FILL );

    glPolygonMode ( GL_BACK, GL_LINE );

    m_ren->init();

    setBackgroundColor(Qt::white);

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();

  }

  glviewer_t::glviewer_t(QWidget * par):
      m_is_recording(false),
      m_bf_cull(true),
      m_wireframe(false),
      m_ren(NULL)
  {
    setParent(par);
  }

  void glviewer_t::setup()
  {
    m_ren = new viewer_t();
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
    if(m_active_otp_idx <0)
      return;

    QModelIndexList l = m_cp_model_proxy->mapSelectionToSource
                        (critpt_view->selectionModel()->selection()).indexes();

    mscomplex_ren_ptr_t msc_ren = glviewer->m_ren->m_mscs[m_active_otp_idx];

    configurable_ctx_menu(msc_ren.get(),l,critpt_view->mapToGlobal(p));

    glviewer->updateGL();
  }

  void viewer_mainwindow::on_datapiece_view_activated ( const QModelIndex & index  )
  {
    if(m_active_otp_idx == index.row())
      return;

    m_active_otp_idx = index.row();

    mscomplex_ren_ptr_t msc_ren = glviewer->m_ren->m_mscs[m_active_otp_idx];

    m_cp_model->reset_configurable(msc_ren.get());
  }

  void viewer_mainwindow::on_actionLoad_Ms_Complex_triggered(bool)
  {
    QString tri_file = QFileDialog::getOpenFileName
        (this,tr("Select triangulation file"),
         QDir::currentPath(),"2D Tri-Mesh (*.tri *.off)");

    if(tri_file == "") return;

    QString ms_file = QFileDialog::getOpenFileName
        (this,tr("Select mscomplex file"),
         QDir::currentPath(),"Mscomplex (*.mscomplex.bin *.mscomplex.full.bin)");

    if(ms_file == "") return;

    load_mscomplex(tri_file,ms_file);
  }

  void viewer_mainwindow::load_mscomplex(QString tf, QString mf)
  {
    mscomplex_ren_ptr_t msc_ren(new mscomplex_ren_t(tf.toStdString(),
                                                    mf.toStdString()));

    msc_ren->init();

    glviewer->m_ren->m_mscs.push_back(msc_ren);

    if(m_cp_model == NULL)
    {
      m_cp_model = new configurable_item_model
                   (msc_ren.get(),this);

      m_cp_model_proxy = new QSortFilterProxyModel(this);

      m_cp_model_proxy->setSourceModel(m_cp_model);

      critpt_view->setModel ( m_cp_model_proxy );

      connect(critpt_filter_edit,SIGNAL(textChanged(QString)),
              m_cp_model_proxy,SLOT(setFilterFixedString(QString)));

      m_active_otp_idx = 0;
    }

    m_cp_model->force_reset();
    m_otp_model->force_reset();
  }

  void viewer_mainwindow::close_mscomplex()
  {
    if(m_active_otp_idx < 0)
      return;


    glviewer->m_ren->m_mscs.erase(glviewer->m_ren->m_mscs.begin() +
                                  m_active_otp_idx);

    m_active_otp_idx = std::max<int>(glviewer->m_ren->m_mscs.size(),
                                     m_active_otp_idx)-1;

    if(m_active_otp_idx < 0)
    {

      disconnect(critpt_filter_edit,SIGNAL(textChanged(QString)),
              m_cp_model_proxy,SLOT(setFilterFixedString(QString)));

      m_cp_model_proxy->setSourceModel(NULL);
      critpt_view->setModel ( NULL );

      delete m_cp_model;
      m_cp_model =NULL;

      delete m_cp_model_proxy;
      m_cp_model_proxy = NULL;
    }
  }

  void viewer_mainwindow::save_snapshot(QString str, int ms)
  {
    glviewer->updateGL();

    glviewer->saveImageSnapshot(str,ms);
  }

  void viewer_mainwindow::on_actionEval_Script_triggered(bool)
  {
    QString fname = QFileDialog::getOpenFileName
        (this,tr("Select script file"),
         QDir::currentPath(),"python (*.py)");

    if(fname == "")
      return;

    eval_script(fname);

    m_pqt_cons->appendCommandPrompt();
  }

  void viewer_mainwindow::eval_script(QString str)
  {
    m_pqt.evalFile(str);
  }

  QList<int> viewer_mainwindow::get_asc_cps(int i)
  {
    QList<int> l;

    i= get_msc_ren()->m_surv_cps[i];

    if(!get_msc_ren())
      return l;

    boost::range::transform
        (get_msc_ren()->m_msc->m_asc_conn[i],back_inserter(l),
         [&](int a){return get_msc_ren()->m_surv_cp_rev[a];});

    return l;
  }

  QList<int> viewer_mainwindow::get_des_cps(int i)
  {
    QList<int> l;

    if(!get_msc_ren())
      return l;

    i = get_msc_ren()->m_surv_cps[i];

    boost::range::transform
        (get_msc_ren()->m_msc->m_des_conn[i],back_inserter(l),
         [&](int a){return get_msc_ren()->m_surv_cp_rev[a];});

    return l;

  }

  void viewer_mainwindow::render_cp(int i)
  {
    if(!get_msc_ren())
      return;

    get_msc_ren()->m_cp_ren_set.insert(i);
  }

  void viewer_mainwindow::set_cp_size(double s)
  {
    g_max_cp_size = s;
  }

  int viewer_mainwindow::get_asc_mfold_size(int i)
  {
    if(!get_msc_ren())
      return 0;

    i = get_msc_ren()->m_surv_cps[i];

    return get_msc_ren()->m_msc->m_asc_mfolds[i].size();
  }



  void viewer_mainwindow::on_actionLoad_Canc_Tree_triggered(bool)
  {
    if(m_active_otp_idx <0)
      return;

    std::string fname = QFileDialog::getOpenFileName
        (this,tr("Select simplified ms complex to load canc tree from"),
         QDir::currentPath(),"Mscomplex (*.mscomplex.bin)").toStdString();

    if(fname == "")
      return;

    mscomplex_t msc;

    msc.load(fname);

    mscomplex_ren_ptr_t msc_ren = glviewer->m_ren->m_mscs[m_active_otp_idx];

    msc_ren->build_canctree(msc.m_canc_list,msc.m_canc_pers);

    msc_ren->update_canctree_tresh(0.0);

    glviewer->updateGL();
  }

  mscomplex_ren_ptr_t viewer_mainwindow::get_msc_ren()
  {
    if(m_active_otp_idx == -1)
      return mscomplex_ren_ptr_t();

    return glviewer->m_ren->m_mscs[m_active_otp_idx];
  }

  void viewer_mainwindow::on_canc_tree_slider_valueChanged ( int value )
  {
    if(m_active_otp_idx <0)
      return;

    mscomplex_ren_ptr_t msc_ren = glviewer->m_ren->m_mscs[m_active_otp_idx];

    msc_ren->update_canctree_tresh(double(value)/99.0);

    glviewer->updateGL();
  }

  viewer_mainwindow::viewer_mainwindow():
      m_active_otp_idx(-1),m_cp_model_proxy(NULL),
      m_cp_model(NULL),m_otp_model(NULL)
  {
    setupUi (this);

    glviewer->setup();

    m_otp_model = new configurable_item_model
                  (glviewer->m_ren,this);

//    spinviewer->setScene(new QGraphicsScene(this));

//    spinviewer->scene()->addItem(new spin::si_graphics_item_t(glviewer->m_ren));

    datapiece_view->setModel ( m_otp_model );

    PythonQt::init();

    m_pqt = PythonQt::self()->getMainModule();

    m_pqt.addObject("ms_mw", this);

    m_pqt_cons = new PythonQtScriptingConsole(data_views_splitter,m_pqt);
  }

  void viewer_mainwindow::showEvent ( QShowEvent * )
  {
    if(m_cp_model)
      m_cp_model->force_reset();
  }

  viewer_mainwindow::~viewer_mainwindow()
  {
  }

  inline QColor to_qcolor (const glutils::color_t & c)
  {
    return QColor::fromRgbF(c[0],c[1],c[2]);
  }

  QSize viewer_mainwindow::msc_conf_dim()
  {
    if(!get_msc_ren())
      return QSize(0,0);

    return QSize(get_msc_ren()->dim()[0],get_msc_ren()->dim()[1]);
  }

  QString viewer_mainwindow::msc_conf_header(int i)
  {
    if(!get_msc_ren())
      return "";

    string str;

    get_msc_ren()->exchange_header(i,str);

    return str.c_str();
  }

  QStringList viewer_mainwindow::msc_conf_headers()
  {
    QStringList ret;

    if(get_msc_ren())
    {
      int ncols = get_msc_ren()->dim()[0];

      for( int i = 0 ; i < ncols; ++i)
      {
        std::string str;

        get_msc_ren()->exchange_header(i,str);

        ret.push_back(str.c_str());
      }
    }

    return ret;
  }


  inline QVariant any_to_qvariant(const boost::any &val)
  {
    if(val.type() == typeid(std::string))
      return QString(boost::any_cast<std::string>(val).c_str());
    else if (val.type() == typeid(bool))
      return boost::any_cast<bool>(val);
    else if (val.type() == typeid(int))
      return boost::any_cast<int>(val);
    else if (val.type() == typeid(glutils::color_t))
      return to_qcolor(boost::any_cast<glutils::color_t>(val));
    else
      return QVariant();
  }

  inline boost::any qvariant_to_any(const QVariant &val)
  {
    boost::any ret;

    switch(val.type())
    {
    case QVariant::Int:    ret = val.toInt();break;
    case QVariant::String: ret = val.toString().toStdString();break;
    case QVariant::Bool:   ret = val.toBool();break;
    case QVariant::Color:
    {
      QColor qc = val.value<QColor>();

      ret = la::make_vec<double>(qc.redF(),qc.greenF(),qc.blueF());
      break;
    }

    }
    return ret;
  }


  QVariant viewer_mainwindow::msc_conf_get_data(int i,int j)
  {
    if(!get_msc_ren())
      return QVariant();

    boost::any v;

    get_msc_ren()->exchange_field(configurable_t::data_index_t(i,j),v);

    return any_to_qvariant(v);
  }

  void viewer_mainwindow::msc_conf_set_data(int i,int j,QVariant var)
  {
    boost::any a = qvariant_to_any(var);

    get_msc_ren()->exchange_field(configurable_t::data_index_t(i,j),a);
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
      std::string h;

      m_conf->exchange_header(m_column_idxs[section],h);

      return h.c_str();
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

    std::string val;

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
      std::string h;

      if(c->exchange_header(i,h) == configurable_t::EFT_DATA_RO)
        continue;

      QAction * action  = m.addAction ( h.c_str());

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
        out_val = la::make_vec<double>(qc.redF(),qc.greenF(),qc.blueF());
    }

    if(m_vals[0].type() == typeid(configurable_t::action_callback_t))
    {      
      out_val = 1; // just some freakin value so it'll get to exchange field
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

//namespace spin
//{
//  QRectF si_graphics_item_t::boundingRect() const
//  {
//    QRectF r;

//    if(!m_si) return r;

//    si_point_t si_lc = m_si->m_iextent.lower_corner();
//    si_point_t si_uc = m_si->m_iextent.upper_corner();

//    r.setBottomLeft(to_qpointf(si_lc));
//    r.setTopRight(to_qpointf(si_uc));

//    return r;
//  }


//  QRgb to_qrgb(si_scalar_t c)
//  {
//    c*= 20;

//    return qRgb(c,c,c);
//  }

//  void si_graphics_item_t::update_image()
//  {
////    if(m_si == m_viewer->m_spin_image) return;

////    m_si = m_viewer->m_spin_image;

////    if(!m_si) return;

////    si_point_t sz = m_si->m_iextent.span();

////    si_point_t lc = m_si->m_iextent.lower_corner();

////    m_image = new QImage(QSize(sz[0],sz[1]),QImage::Format_RGB888);

////    for(uint y = 0 ; y < sz[1]; ++y)
////      for(uint x = 0 ; x < sz[0]; ++x)
////        m_image->setPixel( x,y,to_qrgb((*m_si->m_image)(lc + si_ipoint_t(x,y))));

//  }

//  void si_graphics_item_t::paint
//      (QPainter *painter,
//       const QStyleOptionGraphicsItem *option,
//       QWidget *widget)
//  {
//    update_image();

//    if(!m_si) return;

//    si_point_t lc = m_si->m_iextent.lower_corner();

//    painter->drawImage(QPoint(lc[0],lc[1]),*m_image);
//  }
//}

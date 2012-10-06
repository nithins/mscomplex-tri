#ifndef TRIMESH_VIEWER_MAINWINDOW_INCLUDED
#define TRIMESH_VIEWER_MAINWINDOW_INCLUDED

#include <trimesh.h>

#include <QDialog>
#include <QMainWindow>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>
#include <QTreeView>
#include <QSortFilterProxyModel>
#include <QTimer>

#include <PythonQt.h>
#include <PythonQtScriptingConsole.h>
#include <QGLViewer/qglviewer.h>

#include <boost/any.hpp>

class configurable_t;

namespace trimesh
{


  class viewer_t;

  class data_manager_t;

  class glviewer_t : public QGLViewer
  {

  public:

    viewer_t *m_ren;

    bool m_is_recording;
    bool m_bf_cull;
    bool m_wireframe;

    glviewer_t(QWidget *par);

    void setup();

    ~glviewer_t();

  protected:

    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
    virtual void keyPressEvent(QKeyEvent *e);

  public:
    bool saveImageSnapshot(const QString& ,int );

  };

  class configurable_item_model;

}

#include <ui_trimesh_viewer_mainwindow.h>

namespace trimesh
{

  class mscomplex_ren_t;
  typedef boost::shared_ptr<mscomplex_ren_t> mscomplex_ren_ptr_t;

  class viewer_mainwindow:
      public QMainWindow,
      public Ui::viewer_mainwindow_Dialog
  {

  Q_OBJECT

  public:

    int                      m_active_otp_idx;

    QSortFilterProxyModel   *m_cp_model_proxy;
    configurable_item_model *m_cp_model;
    configurable_item_model *m_otp_model;

    PythonQtObjectPtr          m_pqt;
    PythonQtScriptingConsole  *m_pqt_cons;
  public:

    viewer_mainwindow();

    ~viewer_mainwindow();

    void update_roi_box(double l,double u,uint dim);

    mscomplex_ren_ptr_t get_msc_ren();




  public:

    virtual void showEvent ( QShowEvent * );

  private slots:
    void on_datapiece_view_customContextMenuRequested ( const QPoint &p );

    void on_critpt_view_customContextMenuRequested ( const QPoint &p );

    void on_datapiece_view_activated ( const QModelIndex & index  );

    void on_actionLoad_Canc_Tree_triggered(bool);

    void on_actionLoad_Ms_Complex_triggered(bool);

    void on_canc_tree_slider_valueChanged ( int value );

    void on_actionEval_Script_triggered(bool);

  public slots:

    void load_mscomplex(QString tf , QString mf);

    void close_mscomplex();

    void eval_script(QString str);

    QSize msc_conf_dim();

    QString msc_conf_header(int i);

    QVariant msc_conf_get_data(int i,int j);

    void  msc_conf_set_data(int i,int j,QVariant var);

    void save_snapshot(QString str,int ms);

    QList<int> get_asc_cps(int i);

    int get_asc_mfold_size(int i);

    QList<int> get_des_cps(int i);

    void render_cp(int i);
  };

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p);

  class configurable_ctx_menu_sig_collector:public QObject
  {
    Q_OBJECT

  public:

    std::vector<boost::any>  m_vals;
    int                      m_col;
    configurable_t *         m_conf;
    const std::vector<int> & m_rows;

    configurable_ctx_menu_sig_collector
        (configurable_t * conf,
         const std::vector<boost::any> & vals,
         const int & col,
         const std::vector<int> & rows,
         QObject *par):
        m_conf(conf),
        m_vals(vals),
        m_col(col),
        m_rows(rows)
    {setParent(par);}

  private slots:
    void triggered(bool state);
  };

  class configurable_item_model : public QAbstractTableModel
  {
    Q_OBJECT

  public:

    enum eColumnFilter
    {CF_EFT_DATA_RO = 1,
     CF_EFT_DATA_RW = 2,
     CF_EFT_ACTION  = 4};

    configurable_item_model ( configurable_t *conf,QObject *parent = 0 ):
        QAbstractTableModel ( parent ),m_conf(conf)
    {
      setColumnFilter(CF_EFT_DATA_RO|CF_EFT_DATA_RW);
    }

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    int columnCount ( const QModelIndex &parent = QModelIndex() ) const;

    void reset_configurable(configurable_t *conf);

    void force_reset(){setColumnFilter(m_column_filter);reset();}

    void setColumnFilter(int columnFilter);

  private:

    int m_column_filter;

    std::vector<int> m_column_idxs;

    configurable_t * m_conf;

  };
}

//#include <QGraphicsItem>
//#include <QPainter>
//namespace spin
//{
//  inline QPointF to_qpointf(const si_point_t &p)
//  {
//    return QPointF(p[0],p[1]);
//  }

//  inline QPoint to_qpoint(const si_ipoint_t &p)
//  {
//    return QPoint(p[0],p[1]);
//  }

//  class si_graphics_item_t : public QGraphicsItem
//  {
//  private:
//    trimesh::viewer_t * m_viewer;

//    QImage             *m_image;
//    spin_image_ptr_t    m_si;

//  public:

//    si_graphics_item_t(trimesh::viewer_t* v):m_viewer(v){m_image = NULL;};

//    void update_image();

//    QRectF boundingRect() const;

//    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
//               QWidget *widget);
//  };

//  typedef boost::shared_ptr<si_graphics_item_t> si_graphics_item_ptr_t;
//}



#endif

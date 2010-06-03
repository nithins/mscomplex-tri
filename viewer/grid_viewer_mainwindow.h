#ifndef GRID_VIEWER_MAINWINDOW_INCLUDED
#define GRID_VIEWER_MAINWINDOW_INCLUDED

#include <grid.h>

#include <QDialog>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>
#include <QTreeView>
#include <QSortFilterProxyModel>
#include <QTimer>

#include <QGLViewer/qglviewer.h>

#include <ui_grid_viewer_mainwindow.h>
#include <boost/any.hpp>

class configurable_t;

namespace grid
{


  class grid_viewer_t;

  class data_manager_t;

  class glviewer_t : public QGLViewer
  {

  public:

    grid_viewer_t *m_ren;

    bool m_is_recording;

    glviewer_t(data_manager_t * p);

    ~glviewer_t();

  protected:

    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
    virtual void keyPressEvent(QKeyEvent *e);

  };

  class configurable_item_model;

  class viewer_mainwindow:
      public QDialog,
      public Ui::grid_viewer_mainwindow_Dialog
  {

  Q_OBJECT

  public:

    glviewer_t              *m_viewer;
    uint                     m_active_otp_idx;

    QSortFilterProxyModel   *m_cp_model_proxy;
    configurable_item_model *m_cp_model;
    configurable_item_model *m_otp_model;
    QTimer                  *m_clear_roi_aabb_timer;

  public:

    viewer_mainwindow(data_manager_t *gdm);

    ~viewer_mainwindow();

    void update_roi_box(double l,double u,uint dim);


  public:

    virtual void showEvent ( QShowEvent * );

  private slots:
    void on_datapiece_view_customContextMenuRequested ( const QPoint &p );

    void on_critpt_view_customContextMenuRequested ( const QPoint &p );

    void on_datapiece_view_activated ( const QModelIndex & index  );

    void on_xroi_spanslider_spanChanged(int l , int u );

    void on_yroi_spanslider_spanChanged(int l , int u );

    void on_zroi_spanslider_spanChanged(int l , int u );

    void on_update_roi_pushButton_clicked(bool);

    void on_center_to_roi_checkBox_clicked(bool);

  private slots:
    void clear_roi_aabb();
  };

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p);

  class configurable_ctx_menu_sig_collector:public QObject
  {
    Q_OBJECT

  public:

    boost::any              m_val;
    int                     m_col;
    configurable_t *        m_conf;
    const QModelIndexList & m_rows;

    configurable_ctx_menu_sig_collector
        (configurable_t * conf,
         const boost::any & val,
         const int & col,
         const QModelIndexList & rows,
         QObject *par):
        m_conf(conf),
        m_val(val),
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

    configurable_item_model ( configurable_t *conf,QObject *parent = 0 ):
        QAbstractTableModel ( parent ),m_conf(conf){}

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    int columnCount ( const QModelIndex &parent = QModelIndex() ) const;

    void reset_configurable(configurable_t *conf);

    void force_reset(){reset();}

  private:

    configurable_t * m_conf;

  };
}



#endif

#ifndef TRIMESH_VIEWER_H_INCLUDED
#define TRIMESH_VIEWER_H_INCLUDED

#include <set>
#include <glutils.h>
#include <cpputils.h>
#include <aabb.h>
#include <configurable.h>

#include <trimesh.h>

namespace trimesh
{
  typedef glutils::renderable_ptr_t  renderable_ptr_t;
  typedef glutils::vertex_list_t     vertex_list_t;

  class mscomplex_ren_t:public configurable_t
  {
  public:
    mscomplex_ptr_t m_msc;

    // set externally to control what is rendered
    bool m_bShowCps[gc_max_cell_dim+1];
    bool m_bShowAllCps;
    bool m_bShowCpLabels;
    bool m_bShowMsGraph;
    bool m_bShowGrad;
    bool m_bShowCancCps;
    bool m_bShowCancMsGraph;
    bool m_bShowCellNormals;

    renderable_ptr_t ren_grad[gc_max_cell_dim];
    renderable_ptr_t ren_cp_labels[gc_max_cell_dim+1];
    renderable_ptr_t ren_cp[gc_max_cell_dim+1];
    renderable_ptr_t ren_cp_conns[gc_max_cell_dim];
    renderable_ptr_t ren_canc_cp_labels[gc_max_cell_dim+1];
    renderable_ptr_t ren_canc_cp[gc_max_cell_dim+1];
    renderable_ptr_t ren_canc_cp_conns[gc_max_cell_dim];
    renderable_ptr_t ren_cell_normals;

    glutils::bufobj_ptr_t m_cell_pos_bo;
    glutils::bufobj_ptr_t m_cell_nrm_bo;

    int_list_t                     m_surv_cps;
    std::vector<renderable_ptr_t>  m_surv_mfold_rens[GDIR_CT];
    bool_list_t                    m_surv_mfold_show[GDIR_CT];
    glutils::color_list_t          m_surv_mfold_color[GDIR_CT];
    bool                           m_need_update_geom;

    double            m_extent[6];
    glutils::vertex_t m_center;
    tri_cc_geom_ptr_t m_tcc;

    mscomplex_ren_t(std::string tf,std::string mf);

    void init();
    void render();

    void update_geom();

    // configurable_t interface
  public:
    virtual data_index_t  dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };

  class viewer_t:
      public glutils::renderable_t,
      public configurable_t
  {
  public:
    mscomplex_ren_t       m_msc_ren;
    renderable_ptr_t      m_surf_ren;

  public:
    double m_data_dia;
    bool   m_bShowSurface;

  public:

    viewer_t(std::string tf,std::string gf);

    ~viewer_t();

  private:
    // renderable_t interface
  public:

    void init();

    int  render();

    // configurable_t interface
  public:
    virtual data_index_t  dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };
}
#endif //VIEWER_H_INCLUDED

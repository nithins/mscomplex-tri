#ifndef TRIMESH_VIEWER_H_INCLUDED
#define TRIMESH_VIEWER_H_INCLUDED

#include <set>
#include <glutils.h>
#include <cpputils.h>
#include <aabb.h>
#include <configurable.h>

#include <trimesh.h>
#include <spin_image.h>

namespace trimesh
{
  typedef aabb::aabb_t<double,3>          rect_t;
  typedef aabb::aabb_t<double,3>::point_t point_t;

  class octtree_piece_rendata;

  class viewer_t;

  class disc_rendata_t
  {
  public:

    cellid_t               cellid;
    uint                   index;
    cellid_t               vert_no;

    glutils::renderable_t *ren[GRADDIR_COUNT];
    bool                   show[GRADDIR_COUNT];
    glutils::color_t       color[GRADDIR_COUNT];

    disc_rendata_t(cellid_t c,uint i,cellid_t vno);
    ~disc_rendata_t();

    void render();
    bool update(octtree_piece_rendata *);

    spin::spin_image_ptr_t compute_spin_image(octtree_piece_rendata *,uint dir);

    static void init();
    static void cleanup();
  };

  typedef boost::shared_ptr<glutils::renderable_t> renderable_sp_t;
  typedef boost::shared_ptr<disc_rendata_t>        disc_rendata_sp_t;
  typedef std::set<disc_rendata_sp_t>              disc_rendata_sp_set_t;

  class octtree_piece_rendata:public configurable_t
  {
  public:

    viewer_t      * m_viewer;
    dataset_ptr_t   m_dataset;
    mscomplex_ptr_t m_msgraph;

    // set externally to control what is rendered
    bool m_bShowCps[gc_max_cell_dim+1];
    bool m_bShowAllCps;
    bool m_bShowCpLabels;
    bool m_bShowMsGraph;
    bool m_bShowGrad;
    bool m_bShowCancCps;
    bool m_bShowCancMsGraph;
    bool m_bShowCellNormals;

    // set externally .. cleared by render
    bool m_bNeedUpdateDiscRens;

    renderable_sp_t ren_grad[gc_max_cell_dim];
    renderable_sp_t ren_cp_labels[gc_max_cell_dim+1];
    renderable_sp_t ren_cp[gc_max_cell_dim+1];
    renderable_sp_t ren_cp_conns[gc_max_cell_dim];
    renderable_sp_t ren_canc_cp_labels[gc_max_cell_dim+1];
    renderable_sp_t ren_canc_cp[gc_max_cell_dim+1];
    renderable_sp_t ren_canc_cp_conns[gc_max_cell_dim];
    renderable_sp_t ren_cell_normals;

    // the triangulation

    tri_cc_geom_ptr_t        tri_cc_geom;

    glutils::bufobj_ptr_t    cell_pos_bo;

    glutils::bufobj_ptr_t    cell_nrm_bo;

    std::vector<disc_rendata_sp_t> disc_rds;

    disc_rendata_sp_set_t    active_disc_rens[gc_max_cell_dim + 1];

    rect_t                              m_extent;

    void create_disc_rds();
    void update_active_disc_rens();

    void create_cell_pos_nrm_bo();
    void create_cp_rens(const rect_t &roi);
    void create_canc_cp_rens(const rect_t &roi);
    void create_grad_rens(const rect_t &roi);

    void render_msgraph_data() ;
    void render_dataset_data() ;

    octtree_piece_rendata(dataset_ptr_t, mscomplex_ptr_t,viewer_t *);
    void init(const tri_idx_list_t &,const vertex_list_t &);

    // configurable_t interface
  public:
    virtual data_index_t  dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };

  typedef octtree_piece_rendata * datapiece_rendata_ptr_t;

  class data_manager_t;

  class viewer_t:
      public glutils::renderable_t,
      public configurable_t
  {
  public:
    std::vector<octtree_piece_rendata * >  m_piece_rens;
    rect_t                                 m_roi;
    rect_t                                 m_extent;
    double                                 m_scale_factor;
    point_t                                m_roi_base_pt;
    renderable_sp_t                        m_surf_ren;

  public:
    bool m_bShowRoiBB;
    bool m_bRebuildRens;
    bool m_bCenterToRoi;
    bool m_bShowSurface;

    data_manager_t *                       m_gdm;

    spin::spin_image_ptr_t                 m_spin_image;
  public:

    viewer_t(data_manager_t * );

    ~viewer_t();

    // ensure normalization of l and u, l < u , dim in {0,1,2}
    void set_roi_dim_range_nrm(double l,double u,int dim);

  private:

    void build_rens();

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

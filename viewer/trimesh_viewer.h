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
  typedef glutils::renderable_ptr_t  renderable_ptr_t;
  typedef glutils::vertex_list_t     vertex_list_t;

  class mfold_data_t:public configurable_t
  {
  public:
    std::string              m_mfold_file;
    cellid_list_t            m_cellids;
    int_list_t               m_offsets;

    mscomplex_ptr_t          m_msc;
    int_list_t               m_scpno_cpno_map;

    std::string              m_tri_file;
    tri_cc_geom_ptr_t        m_tcc;
    glutils::bufobj_ptr_t    m_cell_pos_bo;
    glutils::bufobj_ptr_t    m_cell_nrm_bo;
    double                   m_extent[6];

    bool                          m_bNeedUpdate;
    std::vector<renderable_ptr_t> m_ren[GDIR_CT];
    bool_list_t                   m_ren_show[GDIR_CT];
    glutils::color_list_t         m_ren_color[GDIR_CT];

    mfold_data_t(std::string mf,std::string tf,mscomplex_ptr_t msc);

    template<eGDIR dir> inline void read_mfold(int i,cellid_list_t &mfold);
    void __read_mfold(int i,cellid_list_t &mfold);

    void init();
    void render();
    bool update();

    // configurable_t interface
  public:
    virtual data_index_t  dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };

  template<> inline void mfold_data_t::read_mfold<GDIR_DES>(int i,cellid_list_t &mfold)
  {__read_mfold(2*i,mfold);}
  template<> inline void mfold_data_t::read_mfold<GDIR_ASC>(int i,cellid_list_t &mfold)
  {__read_mfold(2*i+1,mfold);}

  class graph_data_t
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

    mfold_data_t m_mfold_data;

    graph_data_t(std::string tf,std::string gf,std::string mf);

    void init();
    void render();

    void init_cps();
    void init_ccps();
  };

  class viewer_t:
      public glutils::renderable_t,
      public configurable_t
  {
  public:
    std::vector<graph_data_t * >  m_graphs;
    renderable_ptr_t              m_surf_ren;

  public:
    double m_data_dia;
    bool m_bShowSurface;

  public:

    viewer_t(std::string tf,std::string gf,std::string mf);

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

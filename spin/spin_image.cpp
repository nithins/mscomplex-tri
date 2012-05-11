#include <spin_image.h>

namespace spin
{
  spin_image_t::spin_image_t():
      m_image(new si_buffer_t){};

  void spin_image_t::init
      (const oriented_point_t & bp, // base pt
       const si_extent_t &e,        // extent
       const si_point_t &r,         // resoulution
       const double &c)             // cutoff angle
  {
    m_iextent = si_iextent_t(e.lower_corner()/r,std::ceil(e.upper_corner()/r));

    m_base_pt = bp;

    m_cutoff_angle = c;

    m_resolution   = r;

    si_ipoint_t  sz = m_iextent.span();

    m_image->resize(boost::extents[sz[0]][sz[1]]);

    m_image->reindex(m_iextent.lower_corner());

    std::fill_n(m_image->data(),sz[0]*sz[1],0);

  }

  void spin_image_t::clear()
  {
    m_image.reset(new si_buffer_t);
  }

  void spin_image_t::accumulate_point(const point_t &x)
  {
    point_t &p  = get_point(m_base_pt);
    point_t &n  = get_normal(m_base_pt);

    double rho  = boost::numeric::ublas::norm_2(cross_product((x-p),n));
    double   z  = boost::numeric::ublas::inner_prod(x-p,n);

    si_point_t  si_pt(rho,z);

    si_pt /= m_resolution;

    if(!m_iextent.contains(si_pt) )  return;

    si_ipoint_t si_ipt = std::floor(si_pt);

    si_point_t bl_crd = si_pt - si_ipt;

    (*m_image)(si_ipt + si_ipoint_t(0,0)) += (1.0 - bl_crd[0])*(1.0 - bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(1,0)) += (0.0 + bl_crd[0])*(1.0 - bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(0,1)) += (1.0 - bl_crd[0])*(0.0 + bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(1,1)) += (0.0 + bl_crd[0])*(0.0 + bl_crd[1]);
  }

}

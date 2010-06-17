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
    m_extent = e;

    m_base_pt = bp;

    m_cutoff_angle = c;

    m_resolution   = r;

    si_ipoint_t img_sz  = std::ceil(e.span()/r);

    si_ipoint_t img_bp  = (e.lower_corner()/r);

    m_image->resize(boost::extents[img_sz[0]][img_sz[1]]);

    m_image->reindex(img_bp);

    std::fill_n(m_image->data(),img_sz[0]*img_sz[1],0);

  }

  void spin_image_t::clear()
  {
    m_image.reset(new si_buffer_t);
  }

  void spin_image_t::accumulate_point(const point_t &x)
  {
    point_t &p  = get_point(m_base_pt);
    point_t &n  = get_normal(m_base_pt);

    double rho  = euclid_norm(cross_product((x-p),n));
    double   z  = dot_product(x-p,n);

    si_point_t  si_pt(rho,z);

    if(!m_extent.contains(si_pt) )  return;

    si_pt /= m_resolution;

    si_ipoint_t si_ipt = std::floor(si_pt);

    si_point_t bl_crd = si_pt - si_ipt;

    (*m_image)(si_ipt + si_ipoint_t(0,0)) = (1.0 - bl_crd[0])*(1.0 - bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(1,0)) = (0.0 + bl_crd[0])*(1.0 - bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(0,1)) = (1.0 - bl_crd[0])*(0.0 + bl_crd[1]);
    (*m_image)(si_ipt + si_ipoint_t(1,1)) = (0.0 + bl_crd[0])*(0.0 + bl_crd[1]);
  }

  QRectF si_graphics_item_t::boundingRect() const
  {
    QRectF r;

    r.setBottomLeft(to_qpointf(m_si->m_extent.lower_corner()));
    r.setTopRight(to_qpointf(m_si->m_extent.upper_corner()));

    return r;
  }

  void si_graphics_item_t::paint
      (QPainter *painter,
       const QStyleOptionGraphicsItem *option,
       QWidget *widget)
  {
    painter->drawRoundedRect(-10, -10, 20, 20, 5, 5);
  }

}

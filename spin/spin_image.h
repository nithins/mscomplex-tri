#ifndef SPIN_IMAGE_H
#define SPIN_IMAGE_H

#include <cpputils.h>
#include <boost/multi_array.hpp>
#include <boost/smart_ptr.hpp>
#include <boost/tuple/tuple.hpp>

#include <glutils.h>
#include <aabb.h>

namespace spin
{
  typedef glutils::vertex_t                       point_t;

  typedef glutils::vertex_list_t                  point_list_t;

  typedef glutils::vertex_t                       normal_t;

  typedef glutils::vertex_list_t                  normal_list_t;

  typedef boost::tuples::tuple<point_t,normal_t>  oriented_point_t;

  typedef double                                  si_scalar_t;

  typedef boost::multi_array<si_scalar_t,2>       si_buffer_t;

  typedef boost::shared_ptr<si_buffer_t>          si_buffer_ptr_t;

  typedef aabb::aabb_t<si_scalar_t,2>             si_extent_t;

  typedef si_extent_t::point_t                    si_point_t;

  typedef int                                     si_int_t;

  typedef aabb::aabb_t<si_int_t,2>                si_iextent_t;

  typedef si_iextent_t::point_t                   si_ipoint_t;

  class spin_image_t
  {
  public:

    inline point_t &        get_point(oriented_point_t & op);

    inline const point_t &  get_point(const oriented_point_t & op);

    inline normal_t &       get_normal(oriented_point_t & op);

    inline const normal_t & get_normal(const oriented_point_t & op);

  private:

    oriented_point_t m_base_pt;

    si_buffer_ptr_t  m_image;

    double           m_cutoff_angle;

    si_point_t       m_resolution;

    si_iextent_t     m_iextent;

  public:

    spin_image_t();

    void init
        (const oriented_point_t & , // base pt
         const si_extent_t &,       // extent
         const si_point_t &,        // resoulution
         const double &);           // cutoff angle

    void clear();

    void accumulate_point(const point_t &);

    friend class si_graphics_item_t;
  };

  typedef boost::shared_ptr<spin_image_t> spin_image_ptr_t;

  inline point_t & spin_image_t::get_point(oriented_point_t & op)
  {
    return boost::tuples::get<0>(op);
  }

  inline const point_t & spin_image_t::get_point(const oriented_point_t & op)
  {
    return boost::tuples::get<0>(op);
  }

  inline normal_t & spin_image_t::get_normal(oriented_point_t & op)
  {
    return boost::tuples::get<1>(op);
  }

  inline const normal_t & spin_image_t::get_normal(const oriented_point_t & op)
  {
    return boost::tuples::get<1>(op);
  }

}
#endif

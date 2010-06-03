#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <cpputils.h>

//#include <boost/static_assert.hpp>
//#define static_assert BOOST_STATIC_ASSERT

namespace grid
{
  const uint gc_grid_dim = 2;

  template <typename coord_type,coord_type invalid_value>
      class   rectangle_complex
  {

  public:

    typedef n_vector_t<coord_type,gc_grid_dim> point_t;

    struct range_t:public n_vector_t<coord_type,2>
    {
      range_t ( const coord_type &l,const coord_type &u)
      {

        (*this)[0] = std::min ( l,u);
        (*this)[1] = std::max ( l,u);
      }

      range_t ():n_vector_t<coord_type,2>(invalid_value,invalid_value){}

      inline bool isInOpen(const coord_type &c) const
      {
        return (( (*this)[0] < c ) && (  c < (*this)[1] ));
      }

      inline bool isInClosed(const coord_type &c) const
      {
        return (( (*this)[0] <= c ) && (  c <= (*this)[1] ));
      }

      inline bool isOnBndry(const coord_type &c) const
      {
        return (( (*this)[0] == c ) || (  c == (*this)[1] ));
      }

      inline bool contains(const range_t & r) const
      {
        return isInOpen(r[0]) && isInOpen(r[1]);
      }

      inline bool intersects(const range_t & r) const
      {
        return !((r[0] > (*this)[1]) || ((*this)[0] > r[1]));
      }

      inline bool intersection(const range_t & r,range_t & i) const
      {
        i = range_t(std::max(r[0],(*this)[0]),std::min(r[1],(*this)[1]));

        return intersects(r);
      }

      inline range_t range_union(const range_t & r) const
      {
        return range_t(std::min(r[0],(*this)[0]),std::max(r[1],(*this)[1]));
      }

      inline coord_type span()
      {
        return ((*this)[1]-(*this)[0]);
      }
    };

    struct rectangle_t:public n_vector_t<range_t,gc_grid_dim>
    {
      typedef n_vector_t<range_t,gc_grid_dim> base_t;

      rectangle_t
          (
              const coord_type & start_x,
              const coord_type & end_x,
              const coord_type & start_y,
              const coord_type & end_y,
              const coord_type & start_z,
              const coord_type & end_z
              )

      {
        (*this)[0] = range_t(start_x,end_x);
        (*this)[1] = range_t(start_y,end_y);
        (*this)[2] = range_t(start_z,end_z);
      }

      rectangle_t
          (
              const range_t &r1,
              const range_t &r2,
              const range_t &r3
              )
      {
        (*this)[0] = r1;
        (*this)[1] = r2;
        (*this)[2] = r3;
      }

      rectangle_t
          (
              const point_t &p1,
              const point_t &p2
              )
      {
        for(uint i = 0 ; i < gc_grid_dim; ++i)
          (*this)[i] = range_t(p1[i],p2[i]);
      }

      rectangle_t(){}

      inline point_t size() const
      {
        point_t ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i][1]-(*this)[i][0];

        return ret;
      }


      bool isInInterior ( const point_t & p ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].isInOpen(p[i]);

        return ret;
      }

      bool contains ( const point_t & p ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].isInClosed(p[i]);

        return ret;
      }

      bool isOnBoundry ( const point_t & p ) const
      {
        return ( contains ( p ) && !isInInterior ( p ) );
      }

      bool contains ( const rectangle_t &r ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].contains(r[i]);

        return ret;
      }

      bool intersects ( const rectangle_t &r ) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].intersects(r[i]);

        return ret;
      }

      bool intersection(const rectangle_t & r,rectangle_t &ixn) const
      {
        bool ret = true;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret &= (*this)[i].intersection(r[i],ixn[i]);

        return ret;
      }

      rectangle_t bounding_box(const rectangle_t & r) const
      {
        rectangle_t ret;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          ret[i] = (*this)[i].range_union(r[i]);

        return ret;
      }

      point_t lower_corner() const
      {
        point_t c;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          c[i]= (*this)[i][0];

        return c;
      }

      point_t upper_corner() const
      {
        point_t c;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          c[i]= (*this)[i][1];

        return c;
      }

      coord_type eff_dim() const
      {
        coord_type d;

        for(size_t i = 0 ; i < base_t::static_size;++i )
          d += ((*this)[i][1] != (*this)[i][0]) ?(1):(0);

        return d;
      }
    };
  };


  typedef int16_t                              cell_coord_t;
  typedef float                                cell_fn_t;
  typedef rectangle_complex<cell_coord_t,-1>   rect_cmplx_t;
  typedef rect_cmplx_t::rectangle_t          rect_t;
  typedef rect_cmplx_t::point_t              cellid_t;
  typedef rect_cmplx_t::point_t              rect_point_t;
  typedef rect_cmplx_t::point_t              rect_size_t;
  typedef rect_cmplx_t::range_t              rect_range_t;
  typedef std::vector<cellid_t>                cellid_list_t;

  enum eGradientDirection
  {
    GRADDIR_DESCENDING,
    GRADDIR_ASCENDING,
    GRADDIR_COUNT,
  };

}

#endif

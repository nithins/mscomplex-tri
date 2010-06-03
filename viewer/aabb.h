#ifndef RECT_H_INCLUDED
#define RECT_H_INCLUDED

#include <cpputils.h>

namespace aabb
{
  template <typename coord_type>
      struct range_t:public n_vector_t<coord_type,2>
  {
    range_t ( const coord_type &l,const coord_type &u)
    {

      (*this)[0] = std::min ( l,u);
      (*this)[1] = std::max ( l,u);
    }

    range_t ():n_vector_t<coord_type,2>(0,0){}

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

  template <typename coord_type,uint max_dim>
      struct aabb_t: public n_vector_t< range_t<coord_type> ,max_dim>
  {
    typedef range_t<coord_type> aabb_range_t;

    typedef n_vector_t<aabb_range_t,max_dim> base_t;

    typedef n_vector_t<coord_type,max_dim> point_t;

    aabb_t(const aabb_range_t &r1,const aabb_range_t &r2,const aabb_range_t &r3)
    {
      (*this)[0] = r1;
      (*this)[1] = r2;
      (*this)[2] = r3;
    }

    aabb_t(const point_t &p1,const point_t &p2)
    {
      for(uint i = 0 ; i < base_t::static_size; ++i)
        (*this)[i] = aabb_range_t(p1[i],p2[i]);
    }

    aabb_t(){}

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

    bool contains ( const aabb_t &r ) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < base_t::static_size;++i )
        ret &= (*this)[i].contains(r[i]);

      return ret;
    }

    bool intersects ( const aabb_t &r ) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < base_t::static_size;++i )
        ret &= (*this)[i].intersects(r[i]);

      return ret;
    }

    bool intersection(const aabb_t & r,aabb_t &ixn) const
    {
      bool ret = true;

      for(size_t i = 0 ; i < base_t::static_size;++i )
        ret &= (*this)[i].intersection(r[i],ixn[i]);

      return ret;
    }

    aabb_t bounding_box(const aabb_t & r) const
    {
      aabb_t ret;

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

    inline point_t span()
    {
      point_t c;

      for(size_t i = 0 ; i < base_t::static_size;++i )
        c[i]= (*this)[i].span();

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
}

#endif

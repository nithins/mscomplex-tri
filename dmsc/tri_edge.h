/***************************************************************************
 *   Copyright (C) 2009 by Nithin Shivashankar,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef __TRI_EDGE_H_INCLUDED__
#define __TRI_EDGE_H_INCLUDED__

#include <boost/shared_ptr.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <utl.h>

class tri_cc_t
{
public:

  static const int cc_dim = 2;
  static const uint INVALID_VALUE;

  struct tri
  {
    uint v;
    uint e;
    uint fnext;
  };

  typedef int                      cellid_t;

  typedef std::vector<cellid_t>    cellid_list_t;
  typedef std::vector<tri>         tri_list_t;

  typedef unsigned int             idx_t;
  typedef la::uivec3_t             tri_idx_t;
  typedef std::vector<tri_idx_t>   tri_idx_list_t;

  tri_list_t    m_tris;   // all versions of all tri.. 3 of each
  cellid_list_t m_verts;   // index list to tris that contain {v1} = {a}
  cellid_list_t m_edges;   // index list to tris that contain {v1,v2} = {a,b}

  tri_cc_t();
  ~tri_cc_t();

  void init(const tri_idx_list_t &,const uint & num_verts);
  void clear();

  void logTri(const uint &qpos , std::ostream &os = std::cout) const ;
  void logTriSet(const uint &trisetstart, std::ostream &os = std::cout) const;

  inline uint vertIndex ( uint t ) const {ASSERT(is_in_range(t,0,m_tris.size()));return m_tris[t].v;}
  inline uint edgeIndex ( uint t ) const {ASSERT(is_in_range(t,0,m_tris.size()));return m_tris[t].e;}
  inline uint triIndex ( uint t ) const  {ASSERT(is_in_range(t,0,m_tris.size()));return t/3;}

  inline uint fnext ( uint t ) const {ASSERT(is_in_range(t,0,m_tris.size()) && has_fnext(t));return m_tris[t].fnext;}
  bool has_fnext ( uint t ) const    {ASSERT(is_in_range(t,0,m_tris.size()));return ! ( m_tris[t].fnext == INVALID_VALUE );}

  inline int vert_ct() const {return m_verts.size();}
  inline int edge_ct() const {return m_edges.size();}
  inline int tri_ct() const {return m_tris.size()/3;}

  uint get_cell_dim(cellid_t c) const ;
  uint get_cell_points(cellid_t  ,cellid_t   * ) const;
  uint get_cell_tris(cellid_t  ,cellid_t   * ) const;
  uint get_cell_facets(cellid_t  ,cellid_t  * ) const;
  uint get_cell_co_facets(cellid_t  ,cellid_t  * ) const;
  uint get_vert_star(cellid_t  ,cellid_t  * ) const;
  uint get_vert_link_verts(cellid_t  ,cellid_t  * ) const;
  cellid_t get_opp_cell(cellid_t c, cellid_t cf) const;

  bool is_adjacent(cellid_t  ,cellid_t ) const;
  bool is_cell_boundry(cellid_t ) const;

  inline uint get_num_cells () const;
  inline uint get_num_cells_dim (uint dim) const;
  inline uint get_num_cells_max_dim (uint dim) const;

  typedef boost::counting_iterator<int> iterator;
  inline iterator begin() const;
  inline iterator end() const;
  inline iterator begin(int i) const;
  inline iterator end(int i) const;

  template<typename Oi>
  inline Oi cellid_to_output(cellid_t ,Oi o);
};

typedef boost::shared_ptr<tri_cc_t> tri_cc_ptr_t;

class tri_cc_geom_t
{

public:

  typedef tri_cc_t::cellid_t       cellid_t;

  typedef tri_cc_t::tri_idx_list_t tri_idx_list_t;

  typedef tri_cc_t::tri_idx_t      tri_idx_t;

  typedef la::dvec3_t              vertex_t;
  typedef std::vector<vertex_t>    vertex_list_t;

  typedef la::dvec3_t              normal_t;
  typedef std::vector<vertex_t>    normal_list_t;

  static const uint cc_dim =       tri_cc_t::cc_dim;

protected:

  tri_cc_ptr_t  m_tri_cc;

  vertex_list_t m_cell_pos;

  normal_list_t m_cell_normal;

  double        m_average_length;

protected:

  double compute_average_length();

public:

  void init(const tri_idx_list_t &,const vertex_list_t &);

  void init(tri_cc_ptr_t ,const vertex_list_t &);

  void clear();

  inline uint get_cell_dim (cellid_t c) const
  { return m_tri_cc->get_cell_dim(c);  }

  inline tri_cc_ptr_t get_tri_cc()
  {return m_tri_cc;}

  inline uint get_cell_points (cellid_t  c,cellid_t   * cl) const
  { return m_tri_cc->get_cell_points(c,cl);  }

  inline uint get_cell_tris (cellid_t  c,cellid_t   * cl) const
  { return m_tri_cc->get_cell_tris(c,cl);  }

  inline uint get_cell_facets (cellid_t  c,cellid_t  * cl) const
  { return m_tri_cc->get_cell_facets(c,cl);  }

  inline uint get_cell_co_facets (cellid_t  c,cellid_t  * cl) const
  { return m_tri_cc->get_cell_co_facets(c,cl);  }

  inline uint get_vert_star(cellid_t  c,cellid_t  * cl) const
  { return m_tri_cc->get_vert_star(c,cl);  }

  inline bool is_cell_boundry(cellid_t c) const
  { return m_tri_cc->is_cell_boundry(c);  }

  inline uint get_num_cells () const
  { return m_tri_cc->get_num_cells();  }

  inline uint get_num_cells_dim (uint d) const
  { return m_tri_cc->get_num_cells_dim(d);  }

  inline uint get_num_cells_max_dim (uint d) const
  { return m_tri_cc->get_num_cells_max_dim(d);  }

  inline vertex_t get_cell_position(cellid_t c) const
  { return m_cell_pos[c]; }

  inline normal_t get_cell_normal(cellid_t c) const
  { return m_cell_normal[c]; }

  inline const vertex_list_t &get_cell_positions() const
  { return m_cell_pos; }

  inline const normal_list_t &get_cell_normals() const
  { return m_cell_normal; }

  inline double get_average_edge_length() const
  { return m_average_length;}

  double get_tri_area(cellid_t c) const;
  double get_vert_area(cellid_t c) const;

};

typedef boost::shared_ptr<tri_cc_geom_t> tri_cc_geom_ptr_t;

inline uint tri_cc_t::get_num_cells_max_dim (uint dim) const
{
  uint n = 0;

  switch(dim)
  {
  case 2: n += tri_ct();
  case 1: n += edge_ct();
  case 0: n += vert_ct();
    return n;
  }
  ASSERT(false&&"invalid dim");return -1;
}

uint tri_cc_t::get_num_cells_dim (uint dim) const
{
  switch(dim)
  {
  case 2: return tri_ct();
  case 1: return edge_ct();
  case 0: return vert_ct();
  }
  ASSERT(false&&"invalid dim");return -1;
}

inline uint tri_cc_t::get_num_cells () const
{return get_num_cells_max_dim(cc_dim);}

inline tri_cc_t::iterator tri_cc_t::begin() const
{return iterator(0);}

inline tri_cc_t::iterator tri_cc_t::end() const
{return iterator(get_num_cells());}

inline tri_cc_t::iterator tri_cc_t::begin(int i) const
{return iterator(get_num_cells_max_dim(i) - get_num_cells_dim(i));}

inline tri_cc_t::iterator tri_cc_t::end(int i) const
{return iterator(get_num_cells_max_dim(i));}

inline uint enext ( uint t ){return ( 3* ( t/3 ) + ( t+1 ) %3 );}
inline uint eprev ( uint t ){ return ( 3* ( t/3 ) + ( t+2 ) %3 );}

template<typename Oi>
inline Oi tri_cc_t::cellid_to_output(cellid_t c,Oi o)
{
  switch(get_cell_dim(c))
  {
  case 0:
    *o++ = c; break;
  case 1:
  {
    cellid_t t  = m_edges[c-vert_ct()];

    *o++ = vertIndex(t);
    *o++ = vertIndex(enext(t));
    *o++ = triIndex(t);
    *o++ = (has_fnext(t))?(triIndex(fnext(t))):(-1);
    break;
  }
  case 2:
    *o++ = c-vert_ct()-edge_ct();
    break;
  }
  return o;
}



#endif

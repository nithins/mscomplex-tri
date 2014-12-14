/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
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

#include <map>
#include <set>
#include <vector>
#include <queue>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <stdexcept>


#include <utl.h>
#include <tri_edge.h>

using namespace std;

const uint tri_cc_t::INVALID_VALUE = 0xffffffff;

tri_cc_t::tri_cc_t(){}

tri_cc_t::~tri_cc_t(){clear();}

typedef la::ivec2_t edge_t;

template<typename T>
inline edge_t mk_edge(const T& u,const T& v)
{edge_t e; e[0] = u;e[1] = v; return e;}

void check_tlist(const tri_cc_t::tri_idx_list_t &tlist,const uint & N);
void check_verts(const tri_cc_t &te);
void check_tris(const tri_cc_t &te);
void check_edges(const tri_cc_t &te);

struct edge_cmp
{
bool operator()(const edge_t &a,const edge_t &b) const
{
  if( a[0] != b[0]) return a[0] <b[0];
  return a[1] <b[1];
}
};

void tri_cc_t::init(const tri_idx_list_t &tlist,const uint & N)
{
  ENSURE(tlist.size() >0 ," No tris!!!");
  ENSURE(N >0 ," No Verts !!!");

  typedef map<edge_t,int,edge_cmp> edge_map_t;

  check_tlist(tlist,N);

  int T = tlist.size();

  edge_cmp cmp;

  edge_map_t edge_map(cmp);

  m_verts.resize(N,INVALID_VALUE);
  m_tris.resize(T*3);

  for ( int ti = 0 ; ti < T; ++ti )
  {
    tri_idx_t t = tlist[ti];

//    if(edge_map.count(mk_edge(t[0],t[1])) == 1 ||
//       edge_map.count(mk_edge(t[1],t[2])) == 1 ||
//       edge_map.count(mk_edge(t[2],t[0])) == 1)
//      std::swap(t[0],t[1]);

    for( int u = 0,v = 1 ; u < 3 ; ++u, v = (v + 1)%3 )
    {
      int tvi = 3*ti + u;

      m_tris[tvi].v = t[u];
      m_verts[t[u]] = tvi;

      edge_t e = mk_edge(t[u],t[v]);
      edge_t e_= mk_edge(t[v],t[u]);

      ENSUREV2(edge_map.count(e) == 0,"2 edges with same induced orientation found",t[u],t[v]);

      if(edge_map.count(e_) == 1)
      {
        int tvj = edge_map[e_];

        ASSERT(m_tris[tvj].fnext == INVALID_VALUE);

        m_tris[tvi].fnext = tvj;
        m_tris[tvj].fnext = tvi;
        m_tris[tvi].e     = m_tris[tvj].e;

        edge_map.erase(e_);
      }
      else
      {
        m_tris[tvi].fnext = INVALID_VALUE;
        m_tris[tvi].e     = m_edges.size();

        m_edges.push_back(tvi);

        edge_map[e] = tvi;
      }
    }
  }

  for(edge_map_t::iterator it = edge_map.begin(); it != edge_map.end(); ++it)
    m_verts[m_tris[it->second].v] = it->second;

  check_verts(*this);
  check_edges(*this);
  check_tris(*this);
}

void tri_cc_t::logTri ( const uint &t , std::ostream &os) const
{
  os<<"#"<<t<<"\t("<<
  m_tris[t].v<<","<<
  m_tris[enext ( t ) ].v<<","<<
  m_tris[enext ( enext ( t ) ) ].v<<")"<<
  "\te = "   <<m_tris[t].e<<
  "\tfnext ="<<m_tris[t].fnext;
  os<<endl;
}

void tri_cc_t::logTriSet ( const uint &trisetstart , std::ostream &os) const
{
  unsigned qstart = 3* ( trisetstart/3 );

  for ( int i = 0 ; i < 3 ; i++ )
  {
    uint qpos = qstart + i;
    logTri ( qpos,os);
  }
  os<<endl;
}

void tri_cc_t::clear()
{
  m_tris.clear();
  m_edges.clear();
  m_verts.clear();
}

uint tri_cc_t::get_cell_dim (cellid_t c) const
{
  if(c < vert_ct())
    return 0;

  c -= vert_ct();

  if(c < edge_ct())
    return 1;

  c -= edge_ct();

  if(c < tri_ct())
    return 2;

  ASSERT(false&&"cellid out of range");
}

uint tri_cc_t::get_cell_points (cellid_t  c,cellid_t   *p ) const
{
  if(c < vert_ct())
  {
    p[0] = c;

    return 1;
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    uint t = m_edges[c];

    p[0] = vertIndex(t); t = enext(t);
    p[1] = vertIndex(t);

    return 2;
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    uint t = c*3;

    p[0] = vertIndex(t);
    p[1] = vertIndex(t+1);
    p[2] = vertIndex(t+2);

    return 3;
  }

  ASSERT(false&&"cellid out of range");
}

uint tri_cc_t::get_cell_tris(cellid_t  c,cellid_t   *p ) const
{
  int tbias = vert_ct() + edge_ct();

  if(c < vert_ct())
  {
    uint tstart = m_verts[c],t = tstart,ct = 0;

    do
    {
      p[ct++] = tbias + triIndex(t);
      t = eprev(t);

      if(!has_fnext(t)) {break;}

      t = fnext(t);
    }
    while(t != tstart);

    return ct;
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    uint t = m_edges[c];

    p[0] = tbias + triIndex(t);
    if(has_fnext(t) )
    {
      t = fnext(t);
      p[1] = tbias + triIndex(t);
      return 2;
    }
    else
    {
      return 1;
    }
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    p[0] = c;
    return 1;
  }

  ASSERT(false&&"cellid out of range");
}

uint tri_cc_t::get_cell_facets (cellid_t  c,cellid_t  * f) const
{
  if(c < vert_ct())
  {
    return 0;
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    uint t = m_edges[c];

    f[0] = vertIndex(t); t = enext(t);
    f[1] = vertIndex(t);

    return 2;
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    uint t = c*3;

    f[0] = vert_ct() + edgeIndex(t);
    f[1] = vert_ct() + edgeIndex(t+1);
    f[2] = vert_ct() + edgeIndex(t+2);

    return 3;
  }

  ASSERT(false&&"cellid out of range");
}

uint tri_cc_t::get_cell_co_facets (cellid_t c ,cellid_t  * cf) const
{
  if(c < vert_ct())
  {
    uint tstart = m_verts[c],t = tstart,cf_ct = 0;

    do
    {
      cf[cf_ct++] = vert_ct() + edgeIndex(t);
      t = eprev(t);

      if(!has_fnext(t))
      {
        cf[cf_ct++] = vert_ct() + edgeIndex(t);
        break;
      }
      t = fnext(t);
    }
    while(t != tstart);

    return cf_ct;
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    uint cf_ct = 0,tri_id_bias = edge_ct()+vert_ct(),t = m_edges[c] ;

    cf[cf_ct++] = tri_id_bias + triIndex(t);

    if(has_fnext(t))
    {
      t = fnext(t);
      cf[cf_ct++] = tri_id_bias + triIndex(t);
    }

    return cf_ct;
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    return 0;
  }

  ASSERT(false&&"cellid out of range");
}

uint tri_cc_t::get_vert_star(cellid_t  c,cellid_t  * cf) const
{
  if(c < vert_ct())
  {
    uint tstart = m_verts[c],t = tstart,cf_ct = 0, toff = edge_ct() + vert_ct();

    cf[cf_ct++] = vert_ct() + edgeIndex(t);

    do
    {
      t = eprev(t);

      cf[cf_ct++] = toff + triIndex(t);
      cf[cf_ct++] = vert_ct() + edgeIndex(t);

      if(!has_fnext(t))
        break;

      t = fnext(t);
    }
    while (t != tstart);

    return cf_ct;
  }

  throw std::runtime_error("invalid vertex id");
}

uint tri_cc_t::get_vert_link_verts(cellid_t  c,cellid_t  * lv) const
{
  if(c < vert_ct())
  {
    uint tstart = m_verts[c],t = tstart,lv_ct = 0;

    do
    {
      lv[lv_ct++] = vertIndex(enext(t));

      t = eprev(t);

      if(!has_fnext(t))
        break;

      t = fnext(t);
    }
    while (t != tstart);

    return lv_ct;
  }

  throw std::runtime_error("invalid vertex id");
}

bool tri_cc_t::is_adjacent(cellid_t  c,cellid_t p) const
{
  if(c > p)
    std::swap(c,p);

  if(c < vert_ct())
  {
    uint tstart = m_verts[c], t = tstart;

    do
    {
      if (p == vert_ct() + edgeIndex(t))
        return true;

      t = eprev(t);

      if(!has_fnext(t))
      {
        if(p == vert_ct() + edgeIndex(t))
          return true;
        break;
      }
      t = fnext ( t );
    }
    while(t != tstart);

    return false;
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    uint toff = edge_ct() + vert_ct(),t = m_edges[c];

    if(p == toff + triIndex(t))
      return true;

    if(has_fnext(t))
    {
      t = fnext(t);
      if (p == toff + triIndex(t))
        return true;
    }

    return false;
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    return false;
  }

  ASSERT(false&&"cellid out of range");
  return false;
}

bool tri_cc_t::is_cell_boundry(cellid_t c) const
{
  if(c < vert_ct())
  {
    return (!has_fnext(m_verts[c]));
  }

  c -= vert_ct();

  if(c < edge_ct())
  {
    return (!has_fnext(m_edges[c]));
  }

  c -= edge_ct();

  if(c < tri_ct())
  {
    return false;
  }

  ASSERT(false&&"cellid out of range");
}

tri_cc_t::cellid_t tri_cc_t::get_opp_cell(cellid_t c, cellid_t cf) const
{
  ASSERT(get_cell_dim(c)+1 == get_cell_dim(cf));

  if(c < vert_ct())
  {
    cf -= vert_ct();

    cellid_t u = m_tris[m_edges[cf]].v;
    cellid_t v = m_tris[enext(m_edges[cf])].v;

    ASSERT(u == c || v == c);

    return (c == u)?(v):(u);
  }

  c -= vert_ct();

  if (c < edge_ct())
  {
    cf -= vert_ct() + edge_ct();

    cellid_t e1 = m_tris[cf*3].e;
    cellid_t e2 = m_tris[cf*3+1].e;
    cellid_t e3 = m_tris[cf*3+2].e;

    cellid_t u = m_tris[cf*3].v;
    cellid_t v = m_tris[cf*3+1].v;
    cellid_t w = m_tris[cf*3+2].v;

    ASSERT( c == e1 || c == e2 || c == e3);

    return (c == e1)?(w):((c == e2)?(u):(v));
  }

  ASSERT(false && "invalid cellid");

}

double tri_cc_geom_t::compute_average_length()
{
  double avg_el = 0;

  for(cellid_t c = get_num_cells_max_dim(0),c_end = get_num_cells_max_dim(1);c != c_end; ++c)
  {
    cellid_t pts[20];

    get_cell_points(c,pts);

    avg_el += (get_cell_position(pts[0])-get_cell_position(pts[1])).norm();
  }

  return avg_el/get_num_cells_dim(1);
}



void tri_cc_geom_t::init(const tri_idx_list_t &tl,const vertex_list_t &vl)
{
  tri_cc_ptr_t tcc(new tri_cc_t);

  tcc->init(tl,vl.size());

  init(tcc,vl);
}

void tri_cc_geom_t::init(tri_cc_ptr_t tcc,const vertex_list_t &vl)
{
  m_tri_cc = tcc;

  m_cell_pos.assign(vl.begin(),vl.end());

  for(cellid_t c = get_num_cells_max_dim(0),c_end = get_num_cells();c != c_end; ++c)
  {
    cellid_t pts[20];

    uint pt_ct = get_cell_points(c,pts);

    vertex_t v; v[0] = 0; v[1] = 0;v[2] = 0;

    for(uint j = 0 ; j < pt_ct; ++j)
      v += m_cell_pos[pts[j]];

    m_cell_pos.push_back(v/pt_ct);
  }

  m_cell_normal.resize(get_num_cells());

  for(cellid_t c = get_num_cells_max_dim(1),c_end = get_num_cells_max_dim(2);c != c_end; ++c)
  {
    cellid_t pts[20];

    get_cell_points(c,pts);

    la::dvec3_t u = m_cell_pos[pts[0]] -m_cell_pos[pts[1]];
    la::dvec3_t v = m_cell_pos[pts[0]] -m_cell_pos[pts[2]];

    ENSURES(u != v)
        << "Found Degenerate triangles in Mesh." << endl
        << STRMVAR(c) << endl
        << STRMVAR(pts[0]) << STRMVAR(pts[1]) << STRMVAR(pts[2]) << endl
        << STRMVAR(m_cell_pos[pts[0]]) << endl
        << STRMVAR(m_cell_pos[pts[1]]) << endl
        << STRMVAR(m_cell_pos[pts[2]]) << endl;


    m_cell_normal[c]  = u.cross(v);

    m_cell_normal[c].normalize();

 }

  cellid_t c = get_num_cells_max_dim(1);

  do
  {
    --c;

    cellid_t cf[20];

    uint cf_ct = get_cell_co_facets(c,cf);

    normal_t n; n[0] = 0; n[1] = 0; n[2] = 0;

    for(uint i = 0; i < cf_ct;++i)
      n += m_cell_normal[cf[i]];

    m_cell_normal[c] = n/cf_ct;
  }
  while(c != 0);

  m_average_length = compute_average_length();
}

void tri_cc_geom_t::clear()
{
  m_tri_cc.reset();

  m_cell_pos.clear();

  m_cell_normal.clear();
}

double tri_cc_geom_t::get_tri_area(cellid_t c) const
{
  ASSERT(get_cell_dim(c) == 2);

  cellid_t pts[20];

  get_cell_points(c,pts);

  return la::tri_area<double,3>
      (m_cell_pos[pts[0]],m_cell_pos[pts[1]],m_cell_pos[pts[2]]);
}

double tri_cc_geom_t::get_vert_area(cellid_t c) const
{
  ASSERT(get_cell_dim(c) == 0);

  double area = 0;

  cellid_t st[40];

  uint st_ct = get_vert_star(c,st);

  for(uint i = 1; i < st_ct; i++)
    area += la::tri_area<double,3>
        (m_cell_pos[st[i-1]],m_cell_pos[c],m_cell_pos[st[i]]);

  return area;
}

// a bunch of routines to ensure that tri_cc_t is consistent.

void check_tlist(const tri_cc_t::tri_idx_list_t &tlist,const uint & N)
{
  for ( int t = 0 ; t < tlist.size(); ++t )
  {
    for( int v = 0 ; v < 3 ; ++v)
    {
      ENSURES(is_in_range(tlist[t][v],0,N))
          << "tri vertex out of range"
          <<STRMVAR(t)
          <<STRMVAR(tlist[t][v])
          <<endl;
    }
  }
}

void check_verts(const tri_cc_t &te)
{
  for( int i = 0 ; i < te.vert_ct(); ++i)
  {
    ENSURES(is_in_range(te.m_verts[i],0,te.m_tris.size()))
        <<"Isolated vertex found"
        <<STRMVAR(i)
        <<endl;
  }
}

void check_tris(const tri_cc_t &te)
{
  vector<int> vert_tri_ct(te.vert_ct(),0);

  for( int i = 0 ; i < te.m_tris.size(); ++i)
  {
    vert_tri_ct[te.m_tris[i].v]++;
  }

  for( int i = 0 ; i < te.vert_ct(); ++i)
  {
    int tstart = te.m_verts[i],t = tstart,t_no = 0;

    do
    {
      t_no++;
      t = eprev(t);
      if (!te.has_fnext(t)) break;
      t = te.fnext(t);
    } while(t != tstart);

    ENSURES(t_no == vert_tri_ct[i])
        << "Non manifold vertex found"
        << STRMVAR(i)
        << endl;
  }
}

void check_edges(const tri_cc_t &te)
{
  edge_cmp cmp;

  set<edge_t,edge_cmp> eset(cmp);

  for( int i = 0 ; i < te.edge_ct(); ++i)
  {
    uint t = te.m_edges[i];

    uint u = te.m_tris[t].v;
    uint v = te.m_tris[enext(t)].v;

    edge_t uv = mk_edge(u,v);

    ENSURES(eset.count(uv) == 0)
        << "2 edges with same induced orientation found"
        <<STRMVAR(u)
        <<STRMVAR(v)
        <<endl;

    eset.insert(uv);
  }
}

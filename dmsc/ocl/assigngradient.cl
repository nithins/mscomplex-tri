////////////////////////////////////////////////////////////////////////////////

// Gradient assignment computation kernel
//

// #include <common_funcs.cl>


const sampler_t vert_fn_sampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;

cell_fn_t get_cell_fn(short2 c,__read_only image2d_t vert_fn_img)
{  
  int2 imgcrd;
  
  imgcrd.x = c.x/2;
  imgcrd.y = c.y/2;
  
  float4 col = read_imagef(vert_fn_img, vert_fn_sampler, imgcrd);
  
  return col.x;     
  
  //return sin(0.0f+0.125f*c.x + 0.5f)*sin(0.0f+0.125f*c.y + 0.5f);
}

int comparePoints(short2 c1,short2 c2,__read_only image2d_t vert_fn_img)
{    
  cell_fn_t f1 = get_cell_fn(c1, vert_fn_img);
  cell_fn_t f2 = get_cell_fn(c2, vert_fn_img);    
    
  if(f1 < f2 ) return 1;
  if(f2 < f1 ) return 0;
  
  if(c1.x < c2.x) return 1;
  if(c2.x < c1.x) return 0;
  
  if(c1.y < c2.y) return 1;   
  return 0;
}

int max_pt_idx(short2 *pts,int num_pts,int offset,__read_only image2d_t vert_fn_img)
{
  int pt_max_idx = offset;
  
  for(int i = offset+1 ; i < num_pts;i++)
    if(comparePoints(pts[pt_max_idx],pts[i], vert_fn_img) == 1)
      pt_max_idx = i;
    
  return pt_max_idx;  
}

int ins_sort_pts(short2 *pts,int num_pts,__read_only image2d_t vert_fn_img)
{  
    for(int i = 0 ; i < num_pts;i++)
    {      
      int swp_idx   = max_pt_idx(pts,num_pts,i,vert_fn_img);            
      short2 temp   = pts[i];
      pts[i]        = pts[swp_idx];
      pts[swp_idx]  = temp;
    }
}

int compareCells(short2 c1,short2 c2,__read_only image2d_t vert_fn_img)
{
  
  short2 pt1[4],pt2[4];
  
  int pt1_ct = get_cell_points(c1,pt1);
  int pt2_ct = get_cell_points(c2,pt2);      
  
  ins_sort_pts(pt1,pt1_ct, vert_fn_img);
  ins_sort_pts(pt2,pt2_ct, vert_fn_img);
  
  int num_lex_comp = min(pt1_ct,pt2_ct);
  
  for(int i = 0 ;i < num_lex_comp;++i)
  {
    if(comparePoints(pt1[i],pt2[i], vert_fn_img) == 1)
      return 1;
    
    if(comparePoints(pt2[i],pt1[i], vert_fn_img) == 1)
      return 0;
  }
  
  if (pt1_ct < pt2_ct) return 1;  
  
  return 0;  
}

__kernel void assign_gradient
( __read_only  image2d_t  vert_fn_img,
  __write_only image2d_t  cell_pr_img,
  __write_only image2d_t  cell_fg_img,
   const short2 int_bl,
   const short2 int_tr,
   const short2 ext_bl,
   const short2 ext_tr
)
   
{  
  short2 c,bb_ext_sz,bb_int_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  bb_int_sz.x = int_tr.x-int_bl.x;
  bb_int_sz.y = int_tr.y-int_bl.y;

  if(get_global_id(0) > bb_int_sz.x ||
     get_global_id(1) > bb_int_sz.y)
   return;

  c.x = get_global_id(0) + int_bl.x - ext_bl.x;
  c.y = get_global_id(1) + int_bl.y - ext_bl.y ;

  int cf_usable[4];

  short2 cf[4];

  int cf_ct = get_cell_cofacets(c,cf);

  int c_is_tb = is_cell_on_true_boundry(c,bb_ext_sz);

  for( int i = 0 ; i < 4;++i)
  {
    cf_usable[i] = 1;

    if(i >= cf_ct)
      cf_usable[i] &= 0;

    cf_usable[i] &= (is_cell_outside_true_boundry(cf[i],bb_ext_sz))?(0):(1);

    int cf_is_tb = is_cell_on_true_boundry(cf[i],bb_ext_sz);

    cf_usable[i] &= ((c_is_tb == 1) && (cf_is_tb == 0))?(0):(1);

    if(cf_usable[i] == 0 )
      continue;

    short2 f[4];

    int f_ct = get_cell_facets(cf[i],f);

    for( int j = 0 ; j < 4;++j)
    {
      if(j >= f_ct)
        continue;

      if(compareCells(c,f[j], vert_fn_img) == 1)
        cf_usable[i] = 0;
    }
  }

  short2 p;
  int is_paired = 0;

  for( int i = 0 ; i < 4;++i)
  {
    if(cf_usable[i] == 1)
    {
      if(is_paired == 0 )
      {
        p = cf[i];
        is_paired = 1;
      }
      else
      {
        if(compareCells(cf[i],p,vert_fn_img) == 1)
          p = cf[i];
      }
    }
  }

  if(is_paired == 1)
  {
    write_cell_pair(c,p,ext_bl,cell_pr_img);

    write_cell_flag(c,1,cell_fg_img);
  }
  else
  {
    write_cell_flag(c,2,cell_fg_img);
  }   
}
__kernel void complete_pairings(
__read_only  image2d_t  cell_pr_img,
__write_only image2d_t  cell_pr_img_out,
__read_only  image2d_t  cell_fg_img,
__write_only image2d_t  cell_fg_img_out,
const short2 ext_bl,
const short2 ext_tr
)
   
{  
  short2 c,bb_ext_sz;
  
  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;
  
  if(get_global_id(0) > bb_ext_sz.x || 
     get_global_id(1) > bb_ext_sz.y)
   return;    
  
  c.x = get_global_id(0);
  c.y = get_global_id(1);
  
  short2 f[4];

  int f_ct = get_cell_facets(c,f);

  short2 p;
  int is_paired = 0;

  for( int i = 0 ; i < f_ct;++i)
  {
   int2 imgcrd;
   imgcrd.y = f[i].x;
   imgcrd.x = f[i].y;

   uint   fflg = get_cell_flag(f[i],cell_fg_img);
   short2 fp   = get_cell_pair(f[i],ext_bl,cell_pr_img);

   if(fp.x == c.x&&
      fp.y == c.y&&
      is_cell_paired(fflg) == 1)
    {
     p = f[i];
     is_paired = 1;
     break;
    }
  }

  if(is_paired == 1)
  {
    write_cell_pair(c,p,ext_bl,cell_pr_img_out);

    write_cell_flag(c,1,cell_fg_img_out);
  }
}

__kernel void mark_boundrypairs_critical_1(
__read_only  image2d_t  cell_pr_img,
__read_only  image2d_t  cell_fg_img,
__write_only image2d_t  cell_fg_img_out,
const short2 int_bl,
const short2 int_tr,
const short2 ext_bl,
const short2 ext_tr
)

{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  if(get_global_id(0) > bb_ext_sz.x ||
     get_global_id(1) > bb_ext_sz.y)
   return;

  c.x = get_global_id(0);
  c.y = get_global_id(1);

  uint flag = get_cell_flag(c,cell_fg_img);

  short2 p  = get_cell_pair(c,ext_bl,cell_pr_img);

  int has_cof_in_ext = 0;

  short2 cf[4];

  int cf_ct = get_cell_cofacets(c,cf);

  short2 rel_int_bl,rel_int_tr;
  rel_int_bl.x = int_bl.x - ext_bl.x;
  rel_int_bl.y = int_bl.y - ext_bl.y;

  rel_int_tr.x = int_tr.x - ext_bl.x;
  rel_int_tr.y = int_tr.y - ext_bl.y;

  if(is_cell_paired(flag))
  {
    for(int i = 0 ;i < 4 ;++i)
    {
      if(i>=cf_ct)
        continue;

      if(is_cell_outside_true_boundry(cf[i],bb_ext_sz) ==1)
        continue;

      has_cof_in_ext |= (cf[i].x < rel_int_bl.x );
      has_cof_in_ext |= (cf[i].y < rel_int_bl.y );
      has_cof_in_ext |= (cf[i].x > rel_int_tr.x );
      has_cof_in_ext |= (cf[i].y > rel_int_tr.y );
    }

    if(has_cof_in_ext == 1)
    {
      write_cell_flag(c,3,cell_fg_img_out);
    }
  }
}

__kernel void mark_boundrypairs_critical_2(
__read_only  image2d_t  cell_pr_img,
__read_only  image2d_t  cell_fg_img,
__write_only image2d_t  cell_fg_img_out,
const short2 ext_bl,
const short2 ext_tr
)
{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  if(get_global_id(0) > bb_ext_sz.x ||
     get_global_id(1) > bb_ext_sz.y)
   return;

  c.x = get_global_id(0);
  c.y = get_global_id(1);

  uint flag = get_cell_flag(c,cell_fg_img);

  short2 p  = get_cell_pair(c,ext_bl,cell_pr_img);

  if(is_cell_paired(flag))
  {
    unsigned int pflag = get_cell_flag(p,cell_fg_img);

    if(is_cell_critical(pflag) && is_cell_paired(pflag))
    {
      write_cell_flag(c,3,cell_fg_img_out);
    }
  }
}
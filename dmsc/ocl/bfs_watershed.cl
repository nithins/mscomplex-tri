#define BLOCK_DIM_X 16
#define BLOCK_DIM_Y 16

// #include <common_funcs.cl>

__kernel void dobfs_markowner_extrema_init(
__read_only  image2d_t  cell_fg_img,
__read_only  image2d_t  cell_pr_img,
__write_only image2d_t  cell_own_image_out,
const short2 ext_bl,
const short2 ext_tr
)
{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  c.x = get_global_id(0);
  c.y = get_global_id(1);

  if(c.x > bb_ext_sz.x ||
     c.y > bb_ext_sz.y)
   return;

  short2 inf_cell;

  inf_cell.x = -1;
  inf_cell.y = -1;

  uint flag = get_cell_flag(c,cell_fg_img);

  if(is_cell_critical(flag) == 1)
    inf_cell  = c;
  
  if(is_cell_critical(flag) == 0 && is_cell_paired(flag) == 1)
  {
    short2 p = get_cell_pair(c,ext_bl,cell_pr_img);

    int d  = get_cell_dim(c);
    int pd = get_cell_dim(p);

    short2 cets[2];

    if(d == 0 || d == 2)
      inf_cell = p;
    else /*if(d ==1)*/
    {
      if(pd == 0)
        get_cell_facets(c,cets);
      else
        get_cell_cofacets(c,cets);

      int inf_idx = 0;

      if(cets[0].x == p.x && cets[0].y == p.y)
        inf_idx = 1;
      
      if(cets[inf_idx].x <= bb_ext_sz.x &&
         cets[inf_idx].y <= bb_ext_sz.y)
        inf_cell = cets[inf_idx];
    }          
  }
  write_to_owner_image(c,inf_cell,cell_own_image_out,ext_bl);
}

__kernel void dobfs_markowner_extrema(
__read_only  image2d_t  cell_own_image_in,
__write_only image2d_t  cell_own_image_out,
__global unsigned int * g_changed,
const short2 ext_bl,
const short2 ext_tr
)
{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  c.x = get_global_id(0);
  c.y = get_global_id(1);

  if(c.x > bb_ext_sz.x ||
     c.y > bb_ext_sz.y)
    return;

  short2 own = read_from_owner_image(c,cell_own_image_in,ext_bl);

  short2 own_own = own;

  if(is_cell_outside_true_boundry(own,bb_ext_sz) == 0) 
    own_own = read_from_owner_image(own,cell_own_image_in,ext_bl);
  
  write_to_owner_image(c,own_own,cell_own_image_out,ext_bl);

  if(own_own.x != own.x || own_own.y != own.y)
    g_changed[0] = 1;
}

#define cell_fn_t    float
#define cell_coord_t short
#define cell_flag_t  uchar

const sampler_t cell_fg_sampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
const sampler_t cell_pr_sampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;


short get_cell_dim(const short2 c)
{
  int dim = (c.x&0x01) + (c.y&0x01);
  return dim;
}

short get_cell_facets(const short2 c,short2 *f)
{
  if(get_cell_dim(c) == 0)
  {
    return 0;
  }
  
  if(get_cell_dim(c) == 1)
  {
    f[0].x = c.x + c.x%2;
    f[0].y = c.y + c.y%2;
    
    f[1].x = c.x - c.x%2;
    f[1].y = c.y - c.y%2;
    return 2;
  }
  
  if(get_cell_dim(c) == 2)
  {
    f[0].x = c.x + 1;
    f[0].y = c.y;
    
    f[1].x = c.x - 1;
    f[1].y = c.y;
    
    f[2].x = c.x;
    f[2].y = c.y + 1;
    
    f[3].x = c.x;
    f[3].y = c.y - 1;
    return 4;
  }
  
  return 0;
}

short get_cell_cofacets(const short2 c,short2 *cf)
{
  if(get_cell_dim(c) == 0)
  {
    cf[0].x = c.x + 1;
    cf[0].y = c.y;
    
    cf[1].x = c.x - 1;
    cf[1].y = c.y;
    
    cf[2].x = c.x;
    cf[2].y = c.y + 1;
    
    cf[3].x = c.x;
    cf[3].y = c.y - 1;
    return 4;
  }
  
  if(get_cell_dim(c) == 1)
  {
    cf[0].x = c.x + c.y%2;
    cf[0].y = c.y + c.x%2;
    
    cf[1].x = c.x - c.y%2;
    cf[1].y = c.y - c.x%2;
    return 2;
  }
  
  if(get_cell_dim(c) == 2)
  {
    return 0;
  }
  
  return 0;
}

short get_cell_points(const short2 c,short2 *p)
{
  if(get_cell_dim(c) == 0)
  {
    p[0] = c;
    return 1;
  }
  
  if(get_cell_dim(c) == 1)
  {
    p[0].x = c.x + c.x%2;
    p[0].y = c.y + c.y%2;
    
    p[1].x = c.x - c.x%2;
    p[1].y = c.y - c.y%2;
    return 2;
  }
  
  if(get_cell_dim(c) == 2)
  {
    p[0].x = c.x + 1;
    p[0].y = c.y + 1;
    
    p[1].x = c.x - 1;
    p[1].y = c.y + 1;
    
    p[2].x = c.x - 1;
    p[2].y = c.y - 1;
    
    p[3].x = c.x + 1;
    p[3].y = c.y - 1;
    return 4;
  }  
  
  return 0;
}

int is_cell_critical(unsigned int flag)
{  
   return (flag&0x02)?1:0;  
}

int is_cell_paired(unsigned int flag)
{
   return (flag&0x01)?1:0;
}

int is_cell_on_true_boundry(short2 c, short2 bb_ext_sz)
{
  int ret = 0;
  ret |= (c.x == 0)?(1):(0);
  ret |= (c.y == 0)?(1):(0);
  ret |= (c.x == bb_ext_sz.x)?(1):(0);
  ret |= (c.y == bb_ext_sz.y)?(1):(0);
  return ret;
}

int is_cell_outside_true_boundry(short2 c, short2 bb_ext_sz)
{
  int ret = 0;
  ret |= (c.x < 0)?(1):(0);
  ret |= (c.y < 0)?(1):(0);
  ret |= (c.x > bb_ext_sz.x)?(1):(0);
  ret |= (c.y > bb_ext_sz.y)?(1):(0);
  return ret;
}

unsigned int  get_cell_flag(short2 c, __read_only image2d_t cell_fg_img)
{  
  int2 imgcrd;
  
  imgcrd.x = c.y;
  imgcrd.y = c.x;
  
  uint4 val = read_imageui(cell_fg_img, cell_fg_sampler, imgcrd);
  
  return val.x;
}

void write_cell_flag(short2 c,unsigned int flag ,__write_only image2d_t cell_fg_img)
{
  int2 imgcrd;

  imgcrd.x = c.y;
  imgcrd.y = c.x;

  uint4 val;
  val.x = flag;
  val.y = 0;
  val.z = 0;
  val.w = 0;

  write_imageui(cell_fg_img, imgcrd,val);
}

void write_cell_pair(short2 c,short2 p,short2 ext_bl ,__write_only image2d_t cell_pr_img)
{
  int2 imgcrd;

  imgcrd.x = c.y;
  imgcrd.y = c.x;


  int4 pr;
  pr.x = p.x + ext_bl.x;
  pr.y = p.y + ext_bl.y;
  pr.z = 0;
  pr.w = 0;

  write_imagei(cell_pr_img,imgcrd,pr);
}

short2  get_cell_pair(short2 c,short2 ext_bl, __read_only image2d_t cell_pr_img)
{
  int2 imgcrd;

  imgcrd.x = c.y;
  imgcrd.y = c.x;

  int4 val = read_imagei(cell_pr_img, cell_pr_sampler, imgcrd);

  short2 pr;

  pr.x = val.x - ext_bl.x;
  pr.y = val.y - ext_bl.y;

  return pr;
}

const sampler_t cell_own_sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;

void write_to_owner_image(short2 c,short2 data, __write_only image2d_t cell_own_image,short2 ext_bl)
{
  int2 imgcrd;

  imgcrd.x = c.y;
  imgcrd.y = c.x;

  int4 data_val;

  data_val.x = data.x;
  data_val.y = data.y;
  data_val.z = 0;
  data_val.w = 0;

  if(data.x != -1 && data.y != -1)
  {
    data_val.x += ext_bl.x ;
    data_val.y += ext_bl.y ;
  }

  write_imagei(cell_own_image, imgcrd,data_val);
}

short2 read_from_owner_image(short2 c, __read_only image2d_t cell_own_image,short2 ext_bl)
{
  int2 imgcrd;

  imgcrd.x = c.y;
  imgcrd.y = c.x;

  int4 data_val = read_imagei(cell_own_image,cell_own_sampler,imgcrd);

  short2 own;

  own.x = data_val.x;
  own.y = data_val.y;

  if(data_val.x != -1 && data_val.y != -1)
  {
    own.x -= ext_bl.x ;
    own.y -= ext_bl.y ;
  }

  return own;
}
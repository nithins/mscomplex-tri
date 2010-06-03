#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable
#extension GL_ARB_texture_rectangle: enable

//HEADER_REPLACE_BEGIN

const int is_dual = 1;

//HEADER_REPLACE_END

const float line_sz   = 1.0;

uniform sampler2DRect rawdata_texture;

vec3[2] get_line(vec3 c)
{
  vec3[2] p; vec3  sz = vec3(0,0,0); 
  
  sz.x   = (((int(c.x)&1)^is_dual) == 1)?(line_sz):(0.0);
  sz.z   = (((int(c.z)&1)^is_dual) == 1)?(line_sz):(0.0);

  p[0]   = c - sz;
  p[1]   = c + sz;

  p[0].y = texture2DRect(rawdata_texture, (p[0].xz)/2).x;
  p[1].y = texture2DRect(rawdata_texture, (p[1].xz)/2).x;

  return p;
}

void draw_line(vec3[2] l)
{
  gl_Position = (gl_ModelViewProjectionMatrix*vec4(l[0],1.0)); EmitVertex(); 
  gl_Position = (gl_ModelViewProjectionMatrix*vec4(l[1],1.0)); EmitVertex(); 
  EndPrimitive();  
}

void main()
{
  for(int i=0; i< gl_VerticesIn; i++)
  {
    gl_FrontColor = gl_FrontColorIn[i];

    draw_line(get_line(gl_PositionIn[i].xyz));
  } 
}
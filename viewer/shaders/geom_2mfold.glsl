#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable
#extension GL_ARB_texture_rectangle: enable

//HEADER_REPLACE_BEGIN

const int is_dual = 1;

//HEADER_REPLACE_END

const float small_sz = 0.1;
const float big_sz   = 1.0;

uniform sampler2DRect rawdata_texture;

vec4 diffuse;
vec4 ambient;
vec4 in_color;

vec3 lightDir;
vec3 halfVector;

void set_front_color_xfm(vec4 p1,vec4 p2,vec4 p3)
{
  vec3 v1 = (p1.xyz-p2.xyz);
  vec3 v2 = (p1.xyz-p3.xyz);

  vec3  halfV,viewV,ldir;
  float NdotL,NdotHV;
  vec4  color  = ambient;

  vec3  n  = normalize(cross(v1,v2));
  NdotL = max(dot(n,lightDir),0.0);

  if (NdotL > 0.0)
  {
    halfV = normalize(halfVector);
    NdotHV = max(dot(n,halfV),0.0);
    if(NdotHV > 0.0)
    {
      color += gl_FrontMaterial.specular *gl_LightSource[0].specular *pow(NdotHV,gl_FrontMaterial.shininess);
      color += diffuse * NdotL;
    }
  }
  gl_FrontColor = in_color*color;
}

void set_light_constants()
{
  diffuse    = gl_LightSource[0].diffuse;
  ambient    = gl_LightSource[0].ambient + gl_LightModel.ambient;
  lightDir   = normalize(vec3(gl_LightSource[0].position));
  halfVector = normalize(gl_LightSource[0].halfVector.xyz);
}

vec3[4] get_quad(vec3 c)
{
  vec3[4] p;  int pos = 0;  ivec3 i = ivec3(0,0,0);ivec3 ex = ivec3(0,0,0); vec3  sz = vec3(0,0,0);
  
  ex.x     = (int(c.x)&1)^is_dual;
  ex.z     = (int(c.z)&1)^is_dual;
  
  sz.x     = (ex.x == 1)?(big_sz):(small_sz);
  sz.z     = (ex.z == 1)?(big_sz):(small_sz);

  for(i[2] = -1 ; i[2] <= 1 ;i[2]+=2)
      for(i[0] = -1 ; i[0] <= 1 ;i[0]+=2)
      {
        p[pos] 	 = c+i*sz;
        p[pos].y = texture2DRect(rawdata_texture, ((c+i*ex).xz)/2).x;
        pos++;
      }

  return p;
}

void draw_quad_post_xfm(vec4 p1,vec4 p2,vec4 p3,vec4 p4)
{
  set_front_color_xfm(p1,p2,p3);

  gl_Position = p1; EmitVertex();
  gl_Position = p2; EmitVertex();
  gl_Position = p3; EmitVertex();
  EndPrimitive();

  gl_Position = p3; EmitVertex();
  gl_Position = p2; EmitVertex();
  gl_Position = p4; EmitVertex();

  EndPrimitive();
}

void draw_quad(vec3[4] b)
{
  vec4[4] xb;

  for ( int i = 0 ; i <4 ; ++i)
    xb[i] = (gl_ModelViewProjectionMatrix*vec4(b[i],1.0));

  draw_quad_post_xfm(xb[1],xb[0],xb[3],xb[2]);
}

void main()
{
  set_light_constants();

  for(int i=0; i< gl_VerticesIn; i++)
  {
    in_color = gl_FrontColorIn[i];

    draw_quad(get_quad(gl_PositionIn[i].xyz));
  } 
}
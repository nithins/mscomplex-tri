#include <grid_dataset.h>

#include <vector>
#include <queue>

#include <timer.h>

#include <QFile>
#include <prefix_scan.h>
#include <bitonic_sort.h>
#include <logutil.h>

#include <grid_mscomplex.h>

std::string cl_error_code_to_string(cl_int err)
{
  switch(err)
  {
  case CL_DEVICE_NOT_FOUND               :return "CL_DEVICE_NOT_FOUND";
  case CL_DEVICE_NOT_AVAILABLE           :return "CL_DEVICE_NOT_AVAILABLE";
  case CL_COMPILER_NOT_AVAILABLE         :return "CL_COMPILER_NOT_AVAILABLE";
  case CL_MEM_OBJECT_ALLOCATION_FAILURE  :return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
  case CL_OUT_OF_RESOURCES               :return "CL_OUT_OF_RESOURCES";
  case CL_OUT_OF_HOST_MEMORY             :return "CL_OUT_OF_HOST_MEMORY";
  case CL_PROFILING_INFO_NOT_AVAILABLE   :return "CL_PROFILING_INFO_NOT_AVAILABLE";
  case CL_MEM_COPY_OVERLAP               :return "CL_MEM_COPY_OVERLAP";
  case CL_IMAGE_FORMAT_MISMATCH          :return "CL_IMAGE_FORMAT_MISMATCH";
  case CL_IMAGE_FORMAT_NOT_SUPPORTED     :return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
  case CL_BUILD_PROGRAM_FAILURE          :return "CL_BUILD_PROGRAM_FAILURE";
  case CL_MAP_FAILURE                    :return "CL_MAP_FAILURE";
  case CL_INVALID_VALUE                  :return "CL_INVALID_VALUE";
  case CL_INVALID_DEVICE_TYPE            :return "CL_INVALID_DEVICE_TYPE";
  case CL_INVALID_PLATFORM               :return "CL_INVALID_PLATFORM";
  case CL_INVALID_DEVICE                 :return "CL_INVALID_DEVICE";
  case CL_INVALID_CONTEXT                :return "CL_INVALID_CONTEXT";
  case CL_INVALID_QUEUE_PROPERTIES       :return "CL_INVALID_QUEUE_PROPERTIES";
  case CL_INVALID_COMMAND_QUEUE          :return "CL_INVALID_COMMAND_QUEUE";
  case CL_INVALID_HOST_PTR               :return "CL_INVALID_HOST_PTR";
  case CL_INVALID_MEM_OBJECT             :return "CL_INVALID_MEM_OBJECT";
  case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
  case CL_INVALID_IMAGE_SIZE             :return "CL_INVALID_IMAGE_SIZE";
  case CL_INVALID_SAMPLER                :return "CL_INVALID_SAMPLER";
  case CL_INVALID_BINARY                 :return "CL_INVALID_BINARY";
  case CL_INVALID_BUILD_OPTIONS          :return "CL_INVALID_BUILD_OPTIONS";
  case CL_INVALID_PROGRAM                :return "CL_INVALID_PROGRAM";
  case CL_INVALID_PROGRAM_EXECUTABLE     :return "CL_INVALID_PROGRAM_EXECUTABLE";
  case CL_INVALID_KERNEL_NAME            :return "CL_INVALID_KERNEL_NAME";
  case CL_INVALID_KERNEL_DEFINITION      :return "CL_INVALID_KERNEL_DEFINITION";
  case CL_INVALID_KERNEL                 :return "CL_INVALID_KERNEL";
  case CL_INVALID_ARG_INDEX              :return "CL_INVALID_ARG_INDEX";
  case CL_INVALID_ARG_VALUE              :return "CL_INVALID_ARG_VALUE";
  case CL_INVALID_ARG_SIZE               :return "CL_INVALID_ARG_SIZE";
  case CL_INVALID_KERNEL_ARGS            :return "CL_INVALID_KERNEL_ARGS";
  case CL_INVALID_WORK_DIMENSION         :return "CL_INVALID_WORK_DIMENSION";
  case CL_INVALID_WORK_GROUP_SIZE        :return "CL_INVALID_WORK_GROUP_SIZE";
  case CL_INVALID_WORK_ITEM_SIZE         :return "CL_INVALID_WORK_ITEM_SIZE";
  case CL_INVALID_GLOBAL_OFFSET          :return "CL_INVALID_GLOBAL_OFFSET";
  case CL_INVALID_EVENT_WAIT_LIST        :return "CL_INVALID_EVENT_WAIT_LIST";
  case CL_INVALID_EVENT                  :return "CL_INVALID_EVENT";
  case CL_INVALID_OPERATION              :return "CL_INVALID_OPERATION";
  case CL_INVALID_GL_OBJECT              :return "CL_INVALID_GL_OBJECT";
  case CL_INVALID_BUFFER_SIZE            :return "CL_INVALID_BUFFER_SIZE";
  case CL_INVALID_MIP_LEVEL              :return "CL_INVALID_MIP_LEVEL";
  case CL_INVALID_GLOBAL_WORK_SIZE       :return "CL_INVALID_GLOBAL_WORK_SIZE";
  }

  std::stringstream ss;

  ss<<"unknown cl error code = "<<err;

  return ss.str();
}

#define _CHECKCL_ERR_CODE(_ERROR,_MESSAGE)\
if(_ERROR != CL_SUCCESS) throw std::runtime_error(cl_error_code_to_string(_ERROR)+"::" + _MESSAGE);

#define _GET_GLOBAL(s,l) (((s)/(2*(l)))*(2*(l)) + ((((s)%(2*(l))) == 0)?(0):(2*l)))

#ifdef _LOG_TIMINGS

#define _START_TIMER(TIMERVAR) Timer TIMERVAR; TIMERVAR.start();

#define _LOG_TIMER(TIMERVAR,MESSAGE) \
_LOG("t = "<<TIMERVAR.getElapsedTimeInMilliSec()<<" ms "<<MESSAGE);

#define _END_TIMER(TIMERVAR) TIMERVAR.stop();

#else

#define _START_TIMER(TIMERVAR)
#define _LOG_TIMER(TIMERVAR,MESSAGE)
#define _END_TIMER(TIMERVAR)

#endif

cl_device_id s_device_id;             // compute device id
cl_context   s_context;               // compute context
PrefixScan   s_pre_scan;              // prefix scanning program
BitonicSortProgram  s_bi_sort;               // bitonic sorting program

struct oclProgInfo
{
  const char * sourcefile;
  const char * additional_include;
  const char * compilation_flags;
  cl_program   _handle;
};

enum eOclProgram
{
  OCLPROG_BEGIN=0,
  OCLPROG_GRADIENT_ASSIGN=0,
  OCLPROG_COLLATE_CRITPTS,
  OCLPROG_BFS_WATERSHED,
  OCLPROG_END,
};

oclProgInfo s_programs[OCLPROG_END-OCLPROG_BEGIN] = {
  {":/oclsources/assigngradient.cl",":/oclsources/common_funcs.cl","",NULL},
  {":/oclsources/collate_critpts.cl",":/oclsources/common_funcs.cl","",NULL},
  {":/oclsources/bfs_watershed.cl",":/oclsources/common_funcs.cl","",NULL},
};

struct oclKernelSourceInfo
{
  eOclProgram  program_idx;
  const char * kernel_name;
  cl_kernel    _handle;
};

enum eOclKernels
{
  OCLKERN_BEGIN=0,
  OCLKERN_ASSIGN_GRADIENT=0,
  OCLKERN_COMPLETE_PAIRINGS,
  OCLKERN_MARKBOUNDRY_PAIRS_CRITICAL_1,
  OCLKERN_MARKBOUNDRY_PAIRS_CRITICAL_2,
  OCLKERN_COLLATE_CPS_INITCOUNT,
  OCLKERN_COLLATE_CPS_WRITEIDS,
  OCLKERN_COUNT_CRITPT_INCIDENCES,
  OCLKERN_WRITE_CRITPT_INCIDENCES,
  OCLKERN_MARK_OWNER_EXTREMA_INIT,
  OCLKERN_MARK_OWNER_EXTREMA,
  OCLKERN_END,
};

oclKernelSourceInfo s_kernels[OCLKERN_END-OCLKERN_BEGIN] =
{
  {OCLPROG_GRADIENT_ASSIGN,"assign_gradient",NULL},
  {OCLPROG_GRADIENT_ASSIGN,"complete_pairings",NULL},
  {OCLPROG_GRADIENT_ASSIGN,"mark_boundrypairs_critical_1",NULL},
  {OCLPROG_GRADIENT_ASSIGN,"mark_boundrypairs_critical_2",NULL},
  {OCLPROG_COLLATE_CRITPTS,"collate_cps_initcount",NULL},
  {OCLPROG_COLLATE_CRITPTS,"collate_cps_writeids",NULL},
  {OCLPROG_COLLATE_CRITPTS,"count_critpt_incidences",NULL},
  {OCLPROG_COLLATE_CRITPTS,"write_critpt_incidences",NULL},
  {OCLPROG_BFS_WATERSHED,"dobfs_markowner_extrema_init",NULL},
  {OCLPROG_BFS_WATERSHED,"dobfs_markowner_extrema",NULL},
};

const int max_threads_1D   = 128;
const int max_threads_2D_x = 16;
const int max_threads_2D_y = 8;

void compile_cl_program(std::string prog_filename,std::string header_filename,
                        std::string compile_flags,cl_program &prog,cl_context & context,cl_device_id &device_id)
{
  std::string prog_src;

  if(header_filename.size() != 0 )
  {
    QFile head_src_qf ( header_filename.c_str() );
    head_src_qf.open(QIODevice::ReadOnly);

    prog_src = head_src_qf.readAll().constData();
    prog_src += "\n";
  }

  QFile prog_src_qf ( prog_filename.c_str() );
  prog_src_qf.open(QIODevice::ReadOnly);

  prog_src += prog_src_qf.readAll().constData();

  int error_code;               // error code returned from api calls

  const char * prog_src_cptr= prog_src.c_str();

  prog = clCreateProgramWithSource
         (context, 1, & prog_src_cptr,
          NULL, &error_code);

  _CHECKCL_ERR_CODE(error_code,"Failed to create program");

  const char * comp_flags_cptr= (compile_flags.size() >0)?compile_flags.c_str():NULL;

  error_code = clBuildProgram(prog, 0, NULL, comp_flags_cptr, NULL, NULL);

  if(error_code != CL_SUCCESS)
  {

    const size_t def_len = 2048;

    char *buffer = new char[def_len];

    size_t len;

    clGetProgramBuildInfo(prog, device_id, CL_PROGRAM_BUILD_LOG,
                          def_len, buffer, &len);

    if(len+1 > def_len)
    {
      buffer = new char[len+1];

      clGetProgramBuildInfo(prog, device_id, CL_PROGRAM_BUILD_LOG,
                            len+1, buffer, &len);
    }

    std::string buf_str(buffer);

    delete []buffer;

    // Log the binary generated

    size_t bin_size;
    error_code = clGetProgramInfo(prog,CL_PROGRAM_BINARY_SIZES,
                                  sizeof(size_t),&bin_size,NULL);

    if(error_code != CL_SUCCESS)
    {
      char * ptx_buffer = new char[bin_size];

      clGetProgramInfo(prog,CL_PROGRAM_BINARIES,
                       sizeof(ptx_buffer),&ptx_buffer,NULL);

      std::string ptx_filename(prog_src_qf.fileName().toStdString());
      ptx_filename += ".ptx";

//      _LOG_TO_FILE(std::string(ptx_buffer),ptx_filename.c_str());
      delete []ptx_buffer;
    }
    throw std::runtime_error(buf_str);
  }

}


void get_preferred_ocl_platform(cl_platform_id* clSelectedPlatformID)
{
  char chBuffer[1024];
  cl_uint num_platforms;
  cl_platform_id* clPlatformIDs;
  cl_int error_code;
  *clSelectedPlatformID = NULL;

  // Get OpenCL platform count
  error_code = clGetPlatformIDs (0, NULL, &num_platforms);

  _CHECKCL_ERR_CODE(error_code,"Failed clGetPlatformIDs");

  if(num_platforms == 0) error_code = CL_SUCCESS - 5000;

  _CHECKCL_ERR_CODE(error_code,"Failed clGetPlatformIDs");

  // if there's a platform or more, make space for ID's
  clPlatformIDs = new cl_platform_id[num_platforms];

  // get platform info for each platform and trap the NVIDIA platform if found
  error_code = clGetPlatformIDs (num_platforms, clPlatformIDs, NULL);

  _CHECKCL_ERR_CODE(error_code,"Failed clGetPlatformIDs");

  for(cl_uint i = 0; i < num_platforms; ++i)
  {
    error_code = clGetPlatformInfo (clPlatformIDs[i], CL_PLATFORM_NAME, 1024, &chBuffer, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed clGetPlatformIDs");

    if(strstr(chBuffer, "NVIDIA") != NULL)
    {
      *clSelectedPlatformID = clPlatformIDs[i];
      break;
    }
  }

  // default to zeroeth platform if NVIDIA not found
  if(*clSelectedPlatformID == NULL)
  {
    _LOG("WARNING: NVIDIA OpenCL platform not found - defaulting to first platform!\n\n");
    *clSelectedPlatformID = clPlatformIDs[0];
  }

  delete []clPlatformIDs;
}

namespace grid
{


  void dataset_t::init_opencl()
  {

    int error_code;

    cl_platform_id cpPlatform;

    get_preferred_ocl_platform(&cpPlatform);

    error_code = clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU , 1, &s_device_id, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to get device ids");

    s_context = clCreateContext(0, 1, &s_device_id, NULL, NULL, &error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create a compute context");


    for(uint pgm_idx = OCLPROG_BEGIN;pgm_idx != OCLPROG_END; pgm_idx++ )
    {

      compile_cl_program(s_programs[pgm_idx].sourcefile,
                         s_programs[pgm_idx].additional_include,
                         s_programs[pgm_idx].compilation_flags,
                         s_programs[pgm_idx]._handle,
                         s_context,s_device_id);
    }

    for(uint kern_idx = OCLKERN_BEGIN;kern_idx != OCLKERN_END; kern_idx++ )
    {
      s_kernels[kern_idx]._handle =
          clCreateKernel(s_programs[s_kernels[kern_idx].program_idx]._handle,
                         s_kernels[kern_idx].kernel_name, &error_code);

      _CHECKCL_ERR_CODE(error_code,"Failed to create kernel");
    }

    s_pre_scan.init(s_context,s_device_id);

    s_bi_sort.initBitonicSort(s_context,s_device_id);
  }

  void dataset_t::stop_opencl()
  {
    for(uint kern_idx = OCLKERN_BEGIN;kern_idx != OCLKERN_END; kern_idx++ )
    {
      clReleaseKernel(s_kernels[kern_idx]._handle);
    }

    for(uint pgm_idx = OCLPROG_BEGIN;pgm_idx != OCLPROG_END; pgm_idx++ )
    {
      clReleaseProgram(s_programs[pgm_idx]._handle);
    }

    clReleaseContext(s_context);

    s_pre_scan.cleanup();

    s_bi_sort.closeBitonicSort();

  }

  void  dataset_t::create_pair_flag_imgs_ocl()
  {

    rect_size_t sz = m_ext_rect.size();

    size_t cell_img_rgn[3] = {sz[1]+1,sz[0]+1,1};

    int error_code;                       // error code returned from api calls

    cl_image_format cell_pr_imgfmt,cell_fg_imgfmt;

    cell_pr_imgfmt.image_channel_data_type = CL_SIGNED_INT16;
    cell_pr_imgfmt.image_channel_order     = CL_RG;

    cell_fg_imgfmt.image_channel_data_type = CL_UNSIGNED_INT8;
    cell_fg_imgfmt.image_channel_order     = CL_R;

    m_cell_pair_img = clCreateImage2D
                      (s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                       &cell_pr_imgfmt,cell_img_rgn[0],cell_img_rgn[1],0,
                       (*m_cell_pairs).data(),&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create cell pair image");

    m_cell_flag_img = clCreateImage2D
                      (s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                       &cell_fg_imgfmt,cell_img_rgn[0],cell_img_rgn[1],0,
                       (*m_cell_flags).data(),&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create cell flag texture");

  }

  void  dataset_t::read_flag_img_ocl(cl_command_queue &commands)
  {

    rect_size_t sz = m_ext_rect.size();

    size_t cell_img_ogn[3] = {0,0,0};
    size_t cell_img_rgn[3] = {sz[1]+1,sz[0]+1,1};

    int error_code;

    error_code = clEnqueueReadImage( commands, m_cell_flag_img, CL_TRUE,
                                     cell_img_ogn,cell_img_rgn,0,0,
                                     (*m_cell_flags).data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read back cell flag image");
  }

  void  dataset_t::read_own_img_ocl(cl_command_queue &commands)
  {

    rect_size_t sz = m_ext_rect.size();

    size_t cell_img_ogn[3] = {0,0,0};
    size_t cell_img_rgn[3] = {sz[1]+1,sz[0]+1,1};

    int error_code;

    error_code = clEnqueueReadImage( commands, m_cell_own_img, CL_TRUE,
                                     cell_img_ogn,cell_img_rgn,0,0,
                                     (*m_cell_own).data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read back cell flag image");
  }

  void  dataset_t::read_pair_img_ocl(cl_command_queue &commands)
  {

    rect_size_t sz = m_ext_rect.size();

    size_t cell_img_ogn[3] = {0,0,0};
    size_t cell_img_rgn[3] = {sz[1]+1,sz[0]+1,1};

    int error_code;

    error_code = clEnqueueReadImage( commands, m_cell_pair_img, CL_TRUE,
                                     cell_img_ogn,cell_img_rgn,0,0,
                                     (*m_cell_pairs).data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read back cell pair image");

  }

  void  dataset_t::clear_buffers_ocl()
  {
    if(m_cell_pair_img != NULL)
      clReleaseMemObject(m_cell_pair_img);

    if(m_cell_flag_img != NULL)
      clReleaseMemObject(m_cell_flag_img);

    if(m_critical_cells_buf != NULL)
      clReleaseMemObject(m_critical_cells_buf);

    if(m_cell_own_img != NULL)
      clReleaseMemObject(m_cell_own_img);

    m_cell_pair_img = NULL;
    m_cell_flag_img = NULL;
    m_critical_cells_buf = NULL;
    m_cell_own_img = NULL;
  }

  void  dataset_t::work_ocl(bool collect_cps )
  {
    int error_code;                       // error code returned from api calls

    _START_TIMER(timer);

    cl_command_queue commands = clCreateCommandQueue
                                (s_context, s_device_id, 0, &error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create commands queue");

    create_pair_flag_imgs_ocl();

    _LOG_TIMER(timer,"Created data imgs");

    assignGradients_ocl(commands);

    _LOG_TIMER(timer,"Gradient Assignment Done");

    read_pair_img_ocl(commands);
    read_flag_img_ocl(commands);

    _LOG_TIMER(timer,"Read back Pair and flag results");

    collateCritcalPoints_ocl(commands);

    _LOG_TIMER(timer,"Collated crit pts");

    uint it_ct =   assignCellOwnerExtrema_ocl(commands);

    if(collect_cps)
    {
      _LOG_TIMER(timer,"Watershed Done in "<<it_ct<<" iterations");

      collect_saddle_conn_ocl(commands);
    }

    _LOG_TIMER(timer,"Saddle incidence collection Done");

    read_own_img_ocl(commands);

    clear_buffers_ocl();

    error_code = clReleaseCommandQueue(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to release commands queue");
  }


  void dataset_t::assignGradients_ocl(cl_command_queue &commands)
  {
    int error_code;                       // error code returned from api calls

    rect_size_t int_sz = m_rect.size();

    rect_size_t ext_sz = m_ext_rect.size();

    size_t vert_img_rgn[3]  = {(ext_sz[0]>>1)+1,(ext_sz[1]>>1)+1,1};

    size_t local[] = {max_threads_2D_x,max_threads_2D_y};

    size_t global[2];

    global[0] = _GET_GLOBAL(int_sz[0]+1,local[0]) ;
    global[1] = _GET_GLOBAL(int_sz[1]+1,local[0]) ;

    cl_image_format vert_fn_imgfmt;

    vert_fn_imgfmt.image_channel_data_type = CL_FLOAT;
    vert_fn_imgfmt.image_channel_order     = CL_R;

    cl_mem vfn_img_cl = clCreateImage2D
                        (s_context,CL_MEM_READ_ONLY| CL_MEM_COPY_HOST_PTR,
                         &vert_fn_imgfmt,vert_img_rgn[0],vert_img_rgn[1],0,
                         (*m_vert_fns_ref).data(),&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create cell fn texture");

    cl_short2 int_bl,int_tr,ext_bl,ext_tr;

    int_bl.x = m_rect.lower_corner()[0];
    int_bl.y = m_rect.lower_corner()[1];
    int_tr.x = m_rect.upper_corner()[0];
    int_tr.y = m_rect.upper_corner()[1];

    ext_bl.x = m_ext_rect.lower_corner()[0];
    ext_bl.y = m_ext_rect.lower_corner()[1];
    ext_tr.x = m_ext_rect.upper_corner()[0];
    ext_tr.y = m_ext_rect.upper_corner()[1];

    unsigned int a = 0;

    cl_kernel kernel= s_kernels[OCLKERN_ASSIGN_GRADIENT]._handle;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &vfn_img_cl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &int_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &int_tr);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _LOG_VAR(local[0]);
    _LOG_VAR(local[1]);
    _LOG_VAR(global[0]);
    _LOG_VAR(global[1]);

    _CHECKCL_ERR_CODE(error_code,"Failed to set assign grad kern args");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enque assign grad kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute assign grad kernel");

    error_code = clReleaseMemObject(vfn_img_cl);

    _CHECKCL_ERR_CODE(error_code,"Failed to release grad vfn_img_cl");

    a = 0;

    kernel = s_kernels[OCLKERN_COMPLETE_PAIRINGS]._handle;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set complete_pairings arguments");

    global[0] = _GET_GLOBAL(ext_sz[0]+1,local[0]) ;
    global[1] = _GET_GLOBAL(ext_sz[1]+1,local[0]) ;

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue complete_pairings kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute complete pairings kernel");

    a = 0;

    kernel = s_kernels[OCLKERN_MARKBOUNDRY_PAIRS_CRITICAL_1]._handle;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &int_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &int_tr);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set mark_boundrypairs_critical_1 arguments");

    global[0] = _GET_GLOBAL(ext_sz[0]+1,local[0]);
    global[1] = _GET_GLOBAL(ext_sz[1]+1,local[0]);

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue mark_boundrypairs_critical_1 kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute mark_boundrypairs_critical_1 kernel");

    a = 0;

    kernel = s_kernels[OCLKERN_MARKBOUNDRY_PAIRS_CRITICAL_2]._handle;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set mark_boundrypairs_critical_2 arguments");

    global[0] = _GET_GLOBAL(ext_sz[0]+1,local[0]);
    global[1] = _GET_GLOBAL(ext_sz[1]+1,local[0]);

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue mark_boundrypairs_critical_2 kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute mark_boundrypairs_critical_2 kernel");
  }

  void dataset_t::collateCritcalPoints_ocl(cl_command_queue &commands)
  {
    int error_code;

    rect_size_t ext_sz = m_ext_rect.size();

    int grid_size = (ext_sz[0]+1)*(ext_sz[1]+1);

    int critpt_idx_buf_sz = (grid_size+1)*sizeof(uint);

    cl_kernel kernel;

    size_t local[] = {max_threads_2D_x,max_threads_2D_y};

    size_t global[2];

    global[0] = _GET_GLOBAL(ext_sz[0]+1,local[0]) ;
    global[1] = _GET_GLOBAL(ext_sz[1]+1,local[0]) ;

    cl_short2 ext_bl,ext_tr;

    ext_bl.x = m_ext_rect.lower_corner()[0];
    ext_bl.y = m_ext_rect.lower_corner()[1];
    ext_tr.x = m_ext_rect.upper_corner()[0];
    ext_tr.y = m_ext_rect.upper_corner()[1];

    cl_mem critpt_idx_buf =
        clCreateBuffer(s_context,CL_MEM_READ_WRITE,critpt_idx_buf_sz,NULL,&error_code);

    _CHECKCL_ERR_CODE(error_code,"couldnt create prefixsum_buf");

    kernel = s_kernels[OCLKERN_COLLATE_CPS_INITCOUNT]._handle;

    uint a=0;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &critpt_idx_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set collate_critpts_initcount kernel arguments");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue collate_critpts_initcount kernel");

    error_code =  clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute collate_critpts_initcount kernel");

    uint crit_pt_ct = 0 ;

    s_pre_scan.CreatePartialSumBuffers(grid_size+1,s_context);
    s_pre_scan.PreScanBuffer(critpt_idx_buf,critpt_idx_buf,grid_size+1,commands);
    s_pre_scan.ReleasePartialSums();


    _CHECKCL_ERR_CODE(error_code,"Failed to execute kernel");

    error_code = clEnqueueReadBuffer
                 (commands,critpt_idx_buf,CL_TRUE,critpt_idx_buf_sz -sizeof(uint),
                  sizeof(uint),&crit_pt_ct,0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read crit_pt_ct");

    uint crit_pt_id_buf_sz = crit_pt_ct*sizeof(cell_coord_t)*2;

    m_critical_cells_buf = clCreateBuffer(s_context,CL_MEM_READ_WRITE,
                                          crit_pt_id_buf_sz,NULL,&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create crit_pt_id_buf");

    kernel = s_kernels[OCLKERN_COLLATE_CPS_WRITEIDS]._handle;

    a = 0;

    error_code = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &critpt_idx_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_critical_cells_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set collate_cps_writeids' arguments");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue collate_cps_writeids ");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute collate_cps_writeids ");

    m_critical_cells.resize(crit_pt_ct);

    error_code = clEnqueueReadBuffer(commands,m_critical_cells_buf,CL_TRUE,0,
                                     crit_pt_id_buf_sz,m_critical_cells.data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read critpt_id_buf");

    error_code = clReleaseMemObject(critpt_idx_buf);

    _CHECKCL_ERR_CODE(error_code,"Failed to reelase critpt_id_buf");
  }




  template <typename T> struct log_2D_array
  {
    void operator()(T* data,uint sz_x,uint sz_y)
    {
      for(int y = 0;y< sz_y;++y)
      {
        for(int x = 0;x< sz_x;++x)
          std::cout<<data[y*sz_x+x];
        std::cout<<std::endl;
      }
    }
  };


  template<typename T,typename array_logger_t = log_2D_array<T> >

  struct read_n_log_2D_img
  {
    void operator()(cl_mem img,cl_command_queue &commands,
                    uint sz_x,uint sz_y,array_logger_t& log_ftor)
    {
      int error_code;

      size_t cell_img_ogn[3] = {0,0,0};

      size_t cell_img_rgn[3] = {sz_x,sz_y,1};

      T* h_buf = new T[sz_x*sz_y];


      error_code = clEnqueueReadImage( commands, img, CL_TRUE,
                                       cell_img_ogn,cell_img_rgn,0,0,
                                       h_buf,0,NULL,NULL);

      _CHECKCL_ERR_CODE(error_code,"Failed to read cell own image");

      log_ftor(h_buf,sz_x,sz_y);

      delete []h_buf;
    }

    void operator()(cl_mem img,cl_command_queue &commands,
                    uint sz_x,uint sz_y)
    {
      array_logger_t array_logger;
      (*this)(img,commands,sz_x,sz_y,array_logger);
    }
  };

  template<typename T,typename array_logger_t = log_2D_array<T> >
  struct read_n_log_2D_buf
  {
    void operator()(cl_mem img,cl_command_queue &commands,
                    uint sz_x,uint sz_y,array_logger_t& log_ftor)
    {
      int error_code;

      T* h_buf = new T[sz_x*sz_y];

      error_code = clEnqueueReadBuffer
                   ( commands, img, CL_TRUE,0,sz_x*sz_y*sizeof(T),h_buf,0,NULL,NULL);

      _CHECKCL_ERR_CODE(error_code,"Failed to read cell own image");

      log_ftor(h_buf,sz_x,sz_y);

      delete []h_buf;
    }

    void operator()(cl_mem img,cl_command_queue &commands,
                    uint sz_x,uint sz_y)
    {
      array_logger_t array_logger;
      (*this)(img,commands,sz_x,sz_y,array_logger);
    }
  };

  struct check_n_log_array
  {
    cellid_t m_c;

    check_n_log_array(cellid_t c):m_c(c){}

    void operator()(cellid_t* data,uint sz_x,uint sz_y)
    {
      std::cout<<sz_x<<"x"<<sz_y<<std::endl;

      for(int x = 0;x< sz_x;++x)
      {
        for(int y = 0;y< sz_y;++y)
        {
          if(data[y*sz_x+x] == m_c)
            std::cout<<"1 ";
          else
            std::cout<<"0 ";
        }
        std::cout<<std::endl;
      }
    }
  };

  int  dataset_t::assignCellOwnerExtrema_ocl(cl_command_queue &commands)
  {

    cl_kernel kernel;

    rect_size_t ext_sz = m_ext_rect.size();

    size_t cell_img_rgn[3] = {ext_sz[1]+1,ext_sz[0]+1,1};

    int error_code;                       // error code returned from api calls

    size_t local[] = {max_threads_2D_x,max_threads_2D_y};

    size_t global[2];

    global[0] = _GET_GLOBAL(ext_sz[0]+1,local[0]) ;
    global[1] = _GET_GLOBAL(ext_sz[1]+1,local[0]) ;

    cl_short2 ext_bl,ext_tr;

    ext_bl.x = m_ext_rect.lower_corner()[0];
    ext_bl.y = m_ext_rect.lower_corner()[1];
    ext_tr.x = m_ext_rect.upper_corner()[0];
    ext_tr.y = m_ext_rect.upper_corner()[1];

    cl_image_format cell_own_imgfmt;

    cell_own_imgfmt.image_channel_data_type = CL_SIGNED_INT16;
    cell_own_imgfmt.image_channel_order     = CL_RG;

    cl_mem cell_own_img[2];

    for(int i = 0 ;i < 2; ++i)
    {
      cell_own_img[i] = clCreateImage2D
                        (s_context,CL_MEM_READ_WRITE,
                         &cell_own_imgfmt,cell_img_rgn[0],cell_img_rgn[1],0,
                         NULL,&error_code);
    }

    _CHECKCL_ERR_CODE(error_code,"Failed to create owner image");

    kernel = s_kernels[OCLKERN_MARK_OWNER_EXTREMA_INIT]._handle;

    uint a = 0 ;

    error_code  = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_flag_img);
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &cell_own_img[0]);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set args for dobfs_init kernel");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue dobfs_init kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute dobfs_init kernel");

    uint is_changed = 0;

    uint iteration_ct = 0;

    cl_mem is_changed_buf =
        clCreateBuffer(s_context,CL_MEM_READ_WRITE,sizeof(uint)*2,NULL,&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create is_changed_buf");

    kernel = s_kernels[OCLKERN_MARK_OWNER_EXTREMA]._handle;

    do
    {
      is_changed = 0;

      error_code = clEnqueueWriteBuffer(commands,is_changed_buf,CL_TRUE,
                                        0,sizeof(uint),&is_changed,0,NULL,NULL);

      _CHECKCL_ERR_CODE(error_code,"Failed to write to is_changed_buf");

      // Set the arguments to our compute kernel
      //
      a = 0 ;

      error_code  = 0;
      error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &cell_own_img[iteration_ct%2]);
      error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &cell_own_img[(1+iteration_ct)%2]);
      error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &is_changed_buf);
      error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
      error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);


      _CHECKCL_ERR_CODE(error_code,"Failed to set args for dobfs_markowner_extrema kernel");

      error_code = clEnqueueNDRangeKernel(commands, kernel, 2, NULL,
                                          global, local, 0, NULL, NULL);

      _CHECKCL_ERR_CODE(error_code,"Failed to enqueue dobfs_markowner_extrema kernel");

      error_code = clFinish(commands);

      _CHECKCL_ERR_CODE(error_code,"Failed to execute dobfs_markowner_extrema kernel");

      error_code = clEnqueueReadBuffer(commands,is_changed_buf,CL_TRUE,
                                       0,sizeof(uint),&is_changed,0,NULL,NULL);

      _CHECKCL_ERR_CODE(error_code,"Failed to read to is_changed_buf");

      iteration_ct++;
    }
    while(is_changed == 1);

    m_cell_own_img = cell_own_img[iteration_ct%2];

    error_code = clReleaseMemObject(is_changed_buf);

    _CHECKCL_ERR_CODE(error_code,"Failed to release ischanged buf ");

    error_code = clReleaseMemObject(cell_own_img[(iteration_ct+1)%2]);

    _CHECKCL_ERR_CODE(error_code,"Failed to release unsaved cell_own_img buf ");

    return iteration_ct;
  }
  template <typename T>
      void read_n_log_1D_buf(cl_command_queue &commands,cl_mem &buf,uint sz)
  {
    T* h_buf= new T[sz];

    int error_code;

    error_code = clEnqueueReadBuffer(commands,buf,CL_TRUE,0,sz*sizeof(T),h_buf,0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"failed to read buf");

    for(uint i =0 ; i < sz; ++i)
      std::cout<<h_buf[i]<<" ";
    std::cout<<std::endl;
  }

  void dataset_t::collect_saddle_conn_ocl(cl_command_queue &commands)
  {

    cl_kernel kernel;

    int error_code;                       // error code returned from api calls

    cl_short2 ext_bl,ext_tr;

    ext_bl.x = m_ext_rect.lower_corner()[0];
    ext_bl.y = m_ext_rect.lower_corner()[1];
    ext_tr.x = m_ext_rect.upper_corner()[0];
    ext_tr.y = m_ext_rect.upper_corner()[1];

    uint critpt_ct = m_critical_cells.size();

    uint incidence_ct_buf_sz = (critpt_ct+1)*sizeof(uint);

    size_t local[] = {max_threads_1D};

    size_t global[1];

    global[0] = _GET_GLOBAL(critpt_ct,local[0]) ;

    cl_mem incidence_ct_buf =
        clCreateBuffer(s_context,CL_MEM_READ_WRITE,incidence_ct_buf_sz,
                       NULL,&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create incidence_ct_buf");

    kernel = s_kernels[OCLKERN_COUNT_CRITPT_INCIDENCES]._handle;

    uint a = 0 ;

    error_code  = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_critical_cells_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_own_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &incidence_ct_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(uint), &critpt_ct);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set args for count_critpt_incidences kernel");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 1, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue count_critpt_incidences kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute count_critpt_incidences kernel");

    s_pre_scan.CreatePartialSumBuffers(critpt_ct+1,s_context);
    s_pre_scan.PreScanBuffer(incidence_ct_buf,incidence_ct_buf,critpt_ct+1,commands);
    s_pre_scan.ReleasePartialSums();

    uint incidence_buf_ct = 0;

    error_code = clEnqueueReadBuffer
                 (commands,incidence_ct_buf,CL_TRUE,incidence_ct_buf_sz-sizeof(uint),
                  sizeof(uint),&incidence_buf_ct,0,NULL,NULL);

    uint incidence_buf_sz = incidence_buf_ct*sizeof(uint);

    cl_mem incidence_buf =
        clCreateBuffer(s_context,CL_MEM_READ_WRITE,incidence_buf_sz,
                       NULL,&error_code);

    _CHECKCL_ERR_CODE(error_code,"Failed to create incidence_buf");

    kernel = s_kernels[OCLKERN_WRITE_CRITPT_INCIDENCES]._handle;

    a = 0 ;

    error_code  = 0;
    error_code  = clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_critical_cells_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_own_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &m_cell_pair_img);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &incidence_ct_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_mem), &incidence_buf);
    error_code |= clSetKernelArg(kernel, a++, sizeof(uint), &critpt_ct);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_bl);
    error_code |= clSetKernelArg(kernel, a++, sizeof(cl_short2), &ext_tr);

    _CHECKCL_ERR_CODE(error_code,"Failed to set args for write_critpt_incidences kernel");

    error_code = clEnqueueNDRangeKernel(commands, kernel, 1, NULL,
                                        global, local, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to enqueue write_critpt_incidences kernel");

    error_code = clFinish(commands);

    _CHECKCL_ERR_CODE(error_code,"Failed to execute write_critpt_incidences kernel");

    m_saddle_incidence_idx_offset.resize(critpt_ct+1);

    error_code =clEnqueueReadBuffer
                (commands,incidence_ct_buf,CL_TRUE,0,incidence_ct_buf_sz,
                 m_saddle_incidence_idx_offset.data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read back saddle inc idx offsets");

    m_saddle_incidence_idx.resize(incidence_buf_ct);

    error_code =clEnqueueReadBuffer
                (commands,incidence_buf,CL_TRUE,0,incidence_buf_sz,
                 m_saddle_incidence_idx.data(),0,NULL,NULL);

    _CHECKCL_ERR_CODE(error_code,"Failed to read back saddle inc idxs ");

    clReleaseMemObject(incidence_ct_buf);

    clReleaseMemObject(incidence_buf);
  }

  cellid_t get_cp_cellid(mscomplex_t *msgraph,uint idx)
  {
    return msgraph->m_cps[idx]->cellid;
  }

  static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
  {
    &dataset_t::getCellFacets,
    &dataset_t::getCellCofacets
  };

  inline uint   dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
  {
    inc[0] = cellid_t (c[0]  ,c[1]+1);
    inc[1] = cellid_t (c[0]  ,c[1]-1);
    inc[2] = cellid_t (c[0]-1,c[1]);
    inc[3] = cellid_t (c[0]+1,c[1]);
    return 4;
  }

  void compute_disc_bfs
      (dataset_t *dataset,
       critpt_disc_t *disc,
       cellid_t start_cellId,
       eGradientDirection gradient_dir
       )
  {
    typedef cellid_t id_type;

    std::queue<id_type> cell_queue;

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      id_type top_cell = cell_queue.front();

      cell_queue.pop();

      disc->push_back(top_cell);

      id_type cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( !dataset->isCellCritical ( cets[i] ) )
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            id_type next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  int dataset_t::postMergeFillDiscs(mscomplex_t *msgraph)
  {
    msgraph->add_disc_tracking_seed_cps();

    for(uint i = 0 ; i < msgraph->m_cps.size() ; ++i)
    {
      critpt_t * cp = msgraph->m_cps[i];

//      if(cp->index != 1) continue;

      for(uint dir = 0 ; dir < GRADDIR_COUNT;++dir)
      {
        if(cp->disc[dir].size() == 1)
        {
          cp->disc[dir].clear();
          compute_disc_bfs(this,&cp->disc[dir],cp->cellid,(eGradientDirection)dir);
        }
      }
    }

    return 0;
  }

  void connectCps (mscomplex_t *msgraph,cellid_t c1,cellid_t c2)
  {
    if(msgraph->m_rect.contains(c1) && msgraph->m_rect.contains(c2))
      msgraph->connect_cps(c1,c2);
  }

  void connectCps (mscomplex_t *msgraph,uint i1,uint i2)
  {
    connectCps(msgraph,msgraph->m_cps[i1]->cellid,msgraph->m_cps[i2]->cellid);
  }

  void dataset_t::writeout_connectivity_ocl(mscomplex_t *msgraph)
  {
    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c));
    }

    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      if(!isCellPaired(c))  continue;

      uint cp_idx = i;

      msgraph->m_cps[cp_idx]->is_paired = true;

      msgraph->m_cps[cp_idx]->pair_idx =
          msgraph->m_id_cp_map[getCellPairId(c)];
    }

    for(uint i = 1 ; i < m_saddle_incidence_idx_offset.size();++i)
    {
      for(uint j = m_saddle_incidence_idx_offset[i-1];
          j < m_saddle_incidence_idx_offset[i];++j)
      {
        connectCps(msgraph,m_saddle_incidence_idx[j],i-1);
      }
    }

  }

  inline bool lowestPairableCoFacet
      (dataset_t *dataset,cellid_t cellId,cellid_t& pairid
       )
  {
    cellid_t cofacets[20];
    bool    cofacet_usable[20];

    uint cofacet_count = dataset->getCellCofacets ( cellId,cofacets );

    bool isTrueBoundryCell = dataset->isTrueBoundryCell ( cellId ) ;

    // for each co facet
    for ( uint i = 0 ; i < cofacet_count ; i++ )
    {
      cellid_t facets[20];
      uint facet_count = dataset->getCellFacets ( cofacets[i],facets );

      cofacet_usable[i] = true;

      if ( isTrueBoundryCell &&
           !dataset->isTrueBoundryCell ( cofacets[i] ) )
      {
        cofacet_usable[i] = false;
        continue;
      }

      for ( uint j = 0 ; j < facet_count ; j++ )
      {
        if ( dataset->compareCells ( cellId,facets[j] ))
        {
          cofacet_usable[i] = false;
          break;
        }
      }
    }

    bool pairid_usable = false;

    for ( uint i =0 ; i < cofacet_count;i++ )
    {
      if ( cofacet_usable[i] == false )
        continue;

      if(pairid_usable == false)
      {
        pairid_usable = true;
        pairid = cofacets[i];
        continue;
      }

      if ( dataset->compareCells ( cofacets[i],pairid ) )
        pairid = cofacets[i];

    }
    return pairid_usable;
  }


  void track_gradient_tree_bfs
      (dataset_t *dataset,cellid_t start_cellId,eGradientDirection gradient_dir)
  {
    std::queue<cellid_t> cell_queue;

    // mark here that that cellid has no parent.

    cell_queue.push ( start_cellId );

    while ( !cell_queue.empty() )
    {
      cellid_t top_cell = cell_queue.front();

      cell_queue.pop();

      (*dataset->m_cell_own)(top_cell) = start_cellId;

      cellid_t      cets[20];

      uint cet_ct = ( dataset->*getcets[gradient_dir] ) ( top_cell,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( dataset->isCellCritical ( cets[i] ) )
        {
          //        connectCps(msgraph,start_cellId,cets[i]);
        }
        else
        {
          if ( !dataset->isCellExterior ( cets[i] ) )
          {
            cellid_t next_cell = dataset->getCellPairId ( cets[i] );

            if ( dataset->getCellDim ( top_cell ) ==
                 dataset->getCellDim ( next_cell ) &&
                 next_cell != top_cell )
            {
              (*dataset->m_cell_own)(cets[i]) = start_cellId;

              // mark here that the parent of next cell is top_cell
              cell_queue.push ( next_cell );
            }
          }
        }
      }
    }
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e) :
      m_rect (r),m_ext_rect (e),m_ptcomp(this),
      m_cell_pair_img(NULL),
      m_cell_flag_img(NULL),
      m_critical_cells_buf(NULL),
      m_cell_own_img(NULL),
      m_vert_fns_ref(NULL)
  {

    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
    //
  }

  dataset_t::dataset_t () :
      m_ptcomp(this),
      m_cell_pair_img(NULL),
      m_cell_flag_img(NULL),
      m_critical_cells_buf(NULL),
      m_cell_own_img(NULL)
  {
    m_vert_fns_ref = NULL;
    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_own     = NULL;
  }

  dataset_t::~dataset_t ()
  {
    clear();

    clear_fnref();
  }

  void dataset_t::init()
  {
    rect_size_t   s = m_ext_rect.size();

    m_cell_flags = new cellflag_array_t( (boost::extents[1+s[0]][1+s[1]]));
    m_cell_pairs = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));
    m_cell_own   = new cellpair_array_t( (boost::extents[1+s[0]][1+s[1]]));

    for (int y = 0 ; y<=s[1];++y)
      for (int x = 0 ; x<=s[0];++x)
        (*m_cell_flags)[x][y] = CELLFLAG_UNKNOWN;

    rect_point_t bl = m_ext_rect.lower_corner();

    (*m_cell_flags).reindex (bl);
    (*m_cell_pairs).reindex (bl);
    (*m_cell_own).reindex (bl);
  }

  void  dataset_t::clear()
  {
    if(m_cell_flags != NULL)
      delete m_cell_flags;

    if(m_cell_pairs != NULL)
      delete m_cell_pairs;

    if(m_cell_own != NULL)
      delete m_cell_own;

    m_critical_cells.clear();

    m_saddle_incidence_idx.clear();
    m_saddle_incidence_idx_offset.clear();

    m_cell_flags   = NULL;
    m_cell_pairs   = NULL;
    m_cell_own     = NULL;
  }

  void dataset_t::init_fnref(cell_fn_t * pData)
  {
    rect_size_t   s = m_ext_rect.size();

    if(pData != NULL)
      m_vert_fns_ref =
          new varray_ref_t(pData,boost::extents[1+s[0]/2][1+s[1]/2],
                           boost::fortran_storage_order());

    rect_point_t bl = m_ext_rect.lower_corner();

    if(pData != NULL)
      (*m_vert_fns_ref).reindex (bl/2);

  }

  void dataset_t::clear_fnref()
  {
    if(m_vert_fns_ref != NULL)
      delete m_vert_fns_ref;

    m_vert_fns_ref = NULL;
  }

  cellid_t   dataset_t::getCellPairId (cellid_t c) const
  {
    if ((*m_cell_flags) (c) &CELLFLAG_PAIRED == 0)
      throw std::logic_error ("invalid pair requested");

    return (*m_cell_pairs) (c);
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    if(getCellDim(c1) == 0)
      return ptLt(c1,c2);

    cellid_t pts1[20];
    cellid_t pts2[20];

    uint pts1_ct = getCellPoints ( c1,pts1);
    uint pts2_ct = getCellPoints ( c2,pts2);

    std::sort ( pts1,pts1+pts1_ct,m_ptcomp );
    std::sort ( pts2,pts2+pts2_ct,m_ptcomp);

    return std::lexicographical_compare
        ( pts1,pts1+pts1_ct,pts2,pts2+pts2_ct,
          m_ptcomp );
  }

  cell_fn_t dataset_t::get_cell_fn (cellid_t c) const
  {

    if(m_vert_fns_ref == NULL)
      return 0.0;

    cell_fn_t  fn = 0.0;

    cellid_t pts[20];

    uint pts_ct = getCellPoints (c,pts);

    for (int j = 0 ; j <pts_ct ;++j)
      fn += (*m_vert_fns_ref) (pts[j]/2);

    fn /= pts_ct;

    return fn;
  }

  void dataset_t::set_cell_fn (cellid_t c,cell_fn_t f)
  {
    if (getCellDim (c) != 0)
      throw std::logic_error ("values only for vertices are specified");

    c[0] /=2;

    c[1] /=2;

    (*m_vert_fns_ref) (c) = f;
  }

  uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
  {
    switch (getCellDim (c))
    {
    case 0:
      p[0] = c;
      return 1;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        p[0] = cellid_t (c[0]+d0,c[1]+d1);
        p[1] = cellid_t (c[0]-d0,c[1]-d1);
      }

      return 2;
    case 2:
      p[0] = cellid_t (c[0]+1,c[1]+1);
      p[1] = cellid_t (c[0]+1,c[1]-1);
      p[2] = cellid_t (c[0]-1,c[1]-1);
      p[3] = cellid_t (c[0]-1,c[1]+1);
      return 4;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }
  }

  uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
  {
    switch (getCellDim (c))
    {
    case 0:
      return 0;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        f[0] = cellid_t (c[0]+d0,c[1]+d1);
        f[1] = cellid_t (c[0]-d0,c[1]-d1);
      }

      return 2;
    case 2:
      f[0] = cellid_t (c[0]  ,c[1]+1);
      f[1] = cellid_t (c[0]  ,c[1]-1);
      f[2] = cellid_t (c[0]-1,c[1]);
      f[3] = cellid_t (c[0]+1,c[1]);
      return 4;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }
  }

  uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = 0;

    switch (getCellDim (c))
    {
    case 0:
      cf[0] = cellid_t (c[0]  ,c[1]+1);
      cf[1] = cellid_t (c[0]  ,c[1]-1);
      cf[2] = cellid_t (c[0]-1,c[1]);
      cf[3] = cellid_t (c[0]+1,c[1]);
      cf_ct =  4;
      break;
    case 1:
      {
        cell_coord_t d0 = c[0]&0x01,d1 = c[1]&0x01;
        cf[0] = cellid_t (c[0]+d1,c[1]+d0);
        cf[1] = cellid_t (c[0]-d1,c[1]-d0);
        cf_ct =  2;
      }

      break;
    case 2:
      return 0;
    default:
      throw std::logic_error ("impossible dim");
      return 0;
    }

    // position in cf[] where the next valid cf should be placed
    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;

  }

  uint dataset_t::getMaxCellDim() const
  {
    return 2;
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

  bool dataset_t::isCellMarked (cellid_t c) const
  {
    return ! ((*m_cell_flags) (c) == CELLFLAG_UNKNOWN);
  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_CRITCAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return ((*m_cell_flags) (c) & CELLFLAG_PAIRED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    (*m_cell_pairs) (c) = p;
    (*m_cell_pairs) (p) = c;

    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_PAIRED;
    (*m_cell_flags) (p) = (*m_cell_flags) (p) |CELLFLAG_PAIRED;
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    (*m_cell_flags) (c) = (*m_cell_flags) (c) |CELLFLAG_CRITCAL;
  }

  bool dataset_t::isTrueBoundryCell (cellid_t c) const
  {
    return (m_ext_rect.isOnBoundry (c));
  }

  bool dataset_t::isFakeBoundryCell (cellid_t c) const
  {
    return (m_rect.isOnBoundry (c) && (!m_ext_rect.isOnBoundry (c)));
  }

  bool dataset_t::isCellExterior (cellid_t c) const
  {
    return (!m_rect.contains (c) && m_ext_rect.contains (c));
  }

  std::string dataset_t::getCellFunctionDescription (cellid_t c) const
  {
    std::stringstream ss;

    ( (std::ostream &) ss) <<c;

    return ss.str();

  }

  std::string dataset_t::getCellDescription (cellid_t c) const
  {

    std::stringstream ss;

    ( (std::ostream &) ss) <<c;

    return ss.str();

  }

  void dataset_t::work()
  {
    assignGradients();

    collateCriticalPoints();

    assignCellOwnerExtrema();
  }

  void  dataset_t::assignGradients()
  {

    // determine all the pairings of all cells in m_rect
    for (cell_coord_t y = m_rect.lower_corner()[1]; y <= m_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y),p;

      if (isCellMarked (c))
        continue;

      if (lowestPairableCoFacet (this,c,p))
        pairCells (c,p);
    }

    for (cell_coord_t y = m_rect.lower_corner()[1]; y <= m_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y);

      if (!isCellMarked (c)) markCellCritical (c);
    }

    // mark artificial boundry as critical

    for (cell_coord_t x = m_rect.lower_corner()[0]; x <= m_rect.upper_corner()[0];x += 1)
    {
      cellid_t bcs[] = {cellid_t (x,m_rect.lower_corner()[1]),cellid_t (x,m_rect.upper_corner()[1]) };

      for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
      {
        cellid_t &c = bcs[i];

        if (isCellCritical (c)) continue;

        cellid_t cf[20];

        u_int cf_ct =  getCellCofacets (c,cf);

        for (u_int j = 0 ; j <cf_ct;++j)
        {
          if (isCellExterior (cf[j]))
          {
            markCellCritical (c);
            markCellCritical (getCellPairId (c));
            break;
          }
        }
      }
    }

    for (cell_coord_t y = m_rect.lower_corner()[1] +1; y < m_rect.upper_corner()[1];y += 1)
    {
      cellid_t bcs[] = {cellid_t (m_rect.lower_corner()[0],y),cellid_t (m_rect.upper_corner()[0],y) };

      for (uint i = 0 ; i <sizeof (bcs) /sizeof (cellid_t);++i)
      {
        cellid_t &c = bcs[i];

        if (isCellCritical (c)) continue;

        cellid_t cf[20];

        u_int cf_ct =  getCellCofacets (c,cf);

        for (u_int j = 0 ; j <cf_ct;++j)
        {
          if (isCellExterior (cf[j]))
          {
            markCellCritical (c);
            markCellCritical (getCellPairId (c));
            break;
          }
        }
      }
    }
  }

  void  dataset_t::collateCriticalPoints()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
      cellid_t c (x,y);

      if (isCellCritical (c))
        m_critical_cells.push_back(c);
    }
  }


  void  dataset_t::assignCellOwnerExtrema()
  {
    for (cellid_list_t::iterator it = m_critical_cells.begin() ;
    it != m_critical_cells.end();++it)
    {

      (*m_cell_own)(*it) = *it;

      switch (getCellDim (*it))
      {
      case 0:
        track_gradient_tree_bfs(this,*it,GRADDIR_ASCENDING);
        break;
      case 2:
        track_gradient_tree_bfs(this,*it,GRADDIR_DESCENDING);
        break;
      default:
        break;
      }
    }

  }

  void  dataset_t::writeout_connectivity(mscomplex_t *msgraph)
  {

    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      msgraph->add_critpt(c,getCellDim(c),get_cell_fn(c));
    }

    for (uint i = 0 ; i <m_critical_cells.size(); ++i)
    {
      cellid_t &c = m_critical_cells[i];

      if(!isCellPaired(c))  continue;

      uint cp_idx = i;

      msgraph->m_cps[cp_idx]->is_paired = true;

      msgraph->m_cps[cp_idx]->pair_idx =
          msgraph->m_id_cp_map[getCellPairId(c)];
    }

    for (cellid_list_t::iterator it = m_critical_cells.begin() ;
    it != m_critical_cells.end();++it)
    {
      cellid_t c = *it;

      if(getCellDim(c) == 1)
      {
        cellid_t f[4],cf[4];

        uint f_ct = getCellFacets(c,f);
        uint cf_ct = getCellCofacets(c,cf);

        for(uint i = 0 ; i < f_ct;++i)
        {
          cellid_t f_own_cp = (*m_cell_own)(f[i]);

          if(f_own_cp != cellid_t(-1,-1))
            connectCps(msgraph,c,f_own_cp);
        }

        for(uint i = 0 ; i < cf_ct;++i)
        {
          cellid_t cf_own_cp = (*m_cell_own)(cf[i]);

          if(cf_own_cp != cellid_t(-1,-1))
            connectCps(msgraph,c,cf_own_cp);
        }
      }
    }
  }

  void dataset_t::getCellCoord (cellid_t c,double &x,double &y,double &z)
  {
    x = c[0];
    y = 0;
    z = c[1];

    cellid_t pts[20];

    if(m_ext_rect.contains(c))
    {
      y= get_cell_fn(c);

    }
  }

  void dataset_t::log_flags()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
    {
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
        cellid_t c(x,y);

        int val = (*m_cell_flags)(c);

        std::cout<<val<<" ";
      }
      std::cout<<std::endl;
    }
  }

  void dataset_t::log_pairs()
  {
    for (cell_coord_t y = m_ext_rect.lower_corner()[1]; y <= m_ext_rect.upper_corner()[1];y += 1)
    {
      for (cell_coord_t x = m_ext_rect.lower_corner()[0]; x <= m_ext_rect.upper_corner()[0];x += 1)
      {
        cellid_t c(x,y);
        std::cout<<(*m_cell_pairs)(c)<<" ";
      }
      std::cout<<std::endl;
    }
  }
}

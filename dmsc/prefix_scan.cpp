#include <stdio.h>
#include <prefix_scan.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

#include <QFile>

#define NUM_BANKS       (16)
size_t  GROUP_SIZE      = 256;

PrefixScan::PrefixScan():
    ScanPartialSums(NULL),
    ElementsAllocated(0),
    LevelsAllocated(0)
{
}


char *PrefixScan::LoadProgramSourceFromFile(const char *filename)
{
  struct stat statbuf;
  FILE        *fh;
  char        *source;

  fh = fopen(filename, "r");
  if (fh == 0)
    return 0;

  stat(filename, &statbuf);
  source = (char *) malloc(statbuf.st_size + 1);
  fread(source, statbuf.st_size, 1, fh);
  source[statbuf.st_size] = '\0';

  return source;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

bool IsPowerOfTwo(int n)
{
  return ((n&(n-1))==0) ;
}

int floorPow2(int n)
{
  int exp;
  frexp((float)n, &exp);
  return 1 << (exp - 1);
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int PrefixScan::CreatePartialSumBuffers(unsigned int count,cl_context& ComputeContext)
{
  ElementsAllocated = count;

  unsigned int group_size = GROUP_SIZE;
  unsigned int element_count = count;

  int level = 0;

  do
  {
    unsigned int group_count = (int)fmax(1, (int)ceil((float)element_count / (2.0f * group_size)));
    if (group_count > 1)
    {
      level++;
    }
    element_count = group_count;

  } while (element_count > 1);

  ScanPartialSums = (cl_mem*) malloc(level * sizeof(cl_mem));
  LevelsAllocated = level;
  memset(ScanPartialSums, 0, sizeof(cl_mem) * level);

  element_count = count;
  level = 0;

  do
  {
    unsigned int group_count = (int)fmax(1, (int)ceil((float)element_count / (2.0f * group_size)));
    if (group_count > 1)
    {
      size_t buffer_size = group_count * sizeof(SCAN_T);
      ScanPartialSums[level++] =
          clCreateBuffer(ComputeContext, CL_MEM_READ_WRITE, buffer_size, NULL, NULL);
    }

    element_count = group_count;

  } while (element_count > 1);

  return CL_SUCCESS;
}

void
    PrefixScan::ReleasePartialSums(void)
{
  unsigned int i;
  for (i = 0; i < LevelsAllocated; i++)
  {
    clReleaseMemObject(ScanPartialSums[i]);
  }

  free(ScanPartialSums);
  ScanPartialSums = 0;
  ElementsAllocated = 0;
  LevelsAllocated = 0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

int
    PrefixScan::PreScan(
        size_t *global,
        size_t *local,
        size_t shared,
        cl_mem output_data,
        cl_mem input_data,
        unsigned int n,
        int group_index,
        int base_index,
        cl_command_queue &ComputeCommands)
{
#if DEBUG_INFO
  printf("PreScan: Global[%4d] Local[%4d] Shared[%4d] BlockIndex[%4d] BaseIndex[%4d] Entries[%d]\n",
         (int)global[0], (int)local[0], (int)shared, group_index, base_index, n);
#endif

  unsigned int k = PRESCAN;
  unsigned int a = 0;

  int err = CL_SUCCESS;
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &output_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &input_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, shared,         0);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &group_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &base_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &n);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to set kernel arguments!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  err = CL_SUCCESS;
  err |= clEnqueueNDRangeKernel(ComputeCommands, ComputeKernels[k], 1, NULL, global, local, 0, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to execute kernel!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  return CL_SUCCESS;
}

int
    PrefixScan::PreScanStoreSum(
        size_t *global,
        size_t *local,
        size_t shared,
        cl_mem output_data,
        cl_mem input_data,
        cl_mem partial_sums,
        unsigned int n,
        int group_index,
        int base_index,
        cl_command_queue &ComputeCommands)
{
#if DEBUG_INFO
  printf("PreScan: Global[%4d] Local[%4d] Shared[%4d] BlockIndex[%4d] BaseIndex[%4d] Entries[%d]\n",
         (int)global[0], (int)local[0], (int)shared, group_index, base_index, n);
#endif

  unsigned int k = PRESCAN_STORE_SUM;
  unsigned int a = 0;

  int err = CL_SUCCESS;
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &output_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &input_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &partial_sums);
  err |= clSetKernelArg(ComputeKernels[k],  a++, shared,         0);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &group_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &base_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &n);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to set kernel arguments!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  err = CL_SUCCESS;
  err |= clEnqueueNDRangeKernel(ComputeCommands, ComputeKernels[k], 1, NULL, global, local, 0, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to execute kernel!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  return CL_SUCCESS;
}

int PrefixScan::PreScanStoreSumNonPowerOfTwo(
    size_t *global,
    size_t *local,
    size_t shared,
    cl_mem output_data,
    cl_mem input_data,
    cl_mem partial_sums,
    unsigned int n,
    int group_index,
    int base_index,
    cl_command_queue &ComputeCommands)
{
#if DEBUG_INFO
  printf("PreScanStoreSumNonPowerOfTwo: Global[%4d] Local[%4d] BlockIndex[%4d] BaseIndex[%4d] Entries[%d]\n",
         (int)global[0], (int)local[0], group_index, base_index, n);
#endif

  unsigned int k = PRESCAN_STORE_SUM_NON_POWER_OF_TWO;
  unsigned int a = 0;

  int err = CL_SUCCESS;
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &output_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &input_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &partial_sums);
  err |= clSetKernelArg(ComputeKernels[k],  a++, shared,         0);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &group_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &base_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &n);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to set kernel arguments!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  err = CL_SUCCESS;
  err |= clEnqueueNDRangeKernel(ComputeCommands, ComputeKernels[k], 1, NULL, global, local, 0, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to execute kernel!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  return CL_SUCCESS;
}

int    PrefixScan::PreScanNonPowerOfTwo(
    size_t *global,
    size_t *local,
    size_t shared,
    cl_mem output_data,
    cl_mem input_data,
    unsigned int n,
    int group_index,
    int base_index,
    cl_command_queue &ComputeCommands)
{
#if DEBUG_INFO
  printf("PreScanNonPowerOfTwo: Global[%4d] Local[%4d] BlockIndex[%4d] BaseIndex[%4d] Entries[%d]\n",
         (int)global[0], (int)local[0], group_index, base_index, n);
#endif

  unsigned int k = PRESCAN_NON_POWER_OF_TWO;
  unsigned int a = 0;

  int err = CL_SUCCESS;
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &output_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &input_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, shared,         0);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &group_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &base_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &n);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to set kernel arguments!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  err = CL_SUCCESS;
  err |= clEnqueueNDRangeKernel(ComputeCommands, ComputeKernels[k], 1, NULL, global, local, 0, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to execute kernel!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }
  return CL_SUCCESS;
}

int PrefixScan::UniformAdd(
    size_t *global,
    size_t *local,
    cl_mem output_data,
    cl_mem partial_sums,
    unsigned int n,
    unsigned int group_offset,
    unsigned int base_index,
    cl_command_queue &ComputeCommands)
{
#if DEBUG_INFO
  printf("UniformAdd: Global[%4d] Local[%4d] BlockOffset[%4d] BaseIndex[%4d] Entries[%d]\n",
         (int)global[0], (int)local[0], group_offset, base_index, n);
#endif

  unsigned int k = UNIFORM_ADD;
  unsigned int a = 0;

  int err = CL_SUCCESS;
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &output_data);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_mem), &partial_sums);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(SCAN_T),  0);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &group_offset);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &base_index);
  err |= clSetKernelArg(ComputeKernels[k],  a++, sizeof(cl_int), &n);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to set kernel arguments!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  err = CL_SUCCESS;
  err |= clEnqueueNDRangeKernel(ComputeCommands, ComputeKernels[k], 1, NULL, global, local, 0, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    printf("Error: %s: Failed to execute kernel!\n", KernelNames[k]);
    return EXIT_FAILURE;
  }

  return CL_SUCCESS;
}

int  PrefixScan::PreScanBufferRecursive(
    cl_mem output_data,
    cl_mem input_data,
    int element_count,
    int level,
    cl_command_queue        &ComputeCommands)
{

  int max_group_size      = GROUP_SIZE;
  int max_work_item_count = GROUP_SIZE;

  unsigned int group_size = max_group_size;
  unsigned int group_count = (int)fmax(1.0f, (int)ceil((float)element_count / (2.0f * group_size)));
  unsigned int work_item_count = 0;

  if (group_count > 1)
    work_item_count = group_size;
  else if (IsPowerOfTwo(element_count))
    work_item_count = element_count / 2;
  else
    work_item_count = floorPow2(element_count);

  work_item_count = (work_item_count > max_work_item_count) ? max_work_item_count : work_item_count;

  unsigned int element_count_per_group = work_item_count * 2;
  unsigned int last_group_element_count = element_count - (group_count-1) * element_count_per_group;
  unsigned int remaining_work_item_count = (int)fmax(1.0f, last_group_element_count / 2);
  remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
  unsigned int remainder = 0;
  size_t last_shared = 0;


  if (last_group_element_count != element_count_per_group)
  {
    remainder = 1;

    if(!IsPowerOfTwo(last_group_element_count))
      remaining_work_item_count = floorPow2(last_group_element_count);

    remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
    unsigned int padding = (2 * remaining_work_item_count) / NUM_BANKS;
    last_shared = sizeof(SCAN_T) * (2 * remaining_work_item_count + padding);
  }

  remaining_work_item_count = (remaining_work_item_count > max_work_item_count) ? max_work_item_count : remaining_work_item_count;
  size_t global[] = { (int)fmax(1, group_count - remainder) * work_item_count, 1 };
  size_t local[]  = { work_item_count, 1 };

  unsigned int padding = element_count_per_group / NUM_BANKS;
  size_t shared = sizeof(SCAN_T) * (element_count_per_group + padding);

  cl_mem partial_sums = ScanPartialSums[level];
  int err = CL_SUCCESS;

  if (group_count > 1)
  {
    err = PreScanStoreSum(global, local, shared, output_data, input_data, partial_sums, work_item_count * 2, 0, 0,ComputeCommands);
    if(err != CL_SUCCESS)
      return err;

    if (remainder)
    {
      size_t last_global[] = { 1 * remaining_work_item_count, 1 };
      size_t last_local[]  = { remaining_work_item_count, 1 };

      err = PreScanStoreSumNonPowerOfTwo(
          last_global, last_local, last_shared,
          output_data, input_data, partial_sums,
          last_group_element_count,
          group_count - 1,
          element_count - last_group_element_count,
          ComputeCommands);

      if(err != CL_SUCCESS)
        return err;

    }

    err = PreScanBufferRecursive(partial_sums, partial_sums, group_count, level + 1,ComputeCommands);
    if(err != CL_SUCCESS)
      return err;

    err = UniformAdd(global, local, output_data, partial_sums,  element_count - last_group_element_count, 0, 0,ComputeCommands);
    if(err != CL_SUCCESS)
      return err;

    if (remainder)
    {
      size_t last_global[] = { 1 * remaining_work_item_count, 1 };
      size_t last_local[]  = { remaining_work_item_count, 1 };

      err = UniformAdd(
          last_global, last_local,
          output_data, partial_sums,
          last_group_element_count,
          group_count - 1,
          element_count - last_group_element_count,
          ComputeCommands);

      if(err != CL_SUCCESS)
        return err;
    }
  }
  else if (IsPowerOfTwo(element_count))
  {
    err = PreScan(global, local, shared, output_data, input_data, work_item_count * 2, 0, 0,ComputeCommands);
    if(err != CL_SUCCESS)
      return err;
  }
  else
  {
    err = PreScanNonPowerOfTwo(global, local, shared, output_data, input_data, element_count, 0, 0,ComputeCommands);
    if(err != CL_SUCCESS)
      return err;
  }

  return CL_SUCCESS;
}

void PrefixScan::PreScanBuffer(
    cl_mem output_data,
    cl_mem input_data,
    unsigned int element_count,
    cl_command_queue        &ComputeCommands)
{
  PreScanBufferRecursive(output_data, input_data, element_count, 0,ComputeCommands);
}

void PrefixScan::init(cl_context &ComputeContext,cl_device_id &ComputeDeviceId,
                      int groupsize)
{

  if(groupsize >0)
    GROUP_SIZE = groupsize;

  int err;

  // Load the compute program from disk into a cstring buffer
  //

  QFile prog_src_file ( ":/oclsources/scan_kernel.cl" );
  prog_src_file.open(QIODevice::ReadOnly);

  QByteArray prog_src_ba = prog_src_file.readAll().constData();

  const char * source = prog_src_ba.constData();

  // Create the compute program from the source buffer
  //
  ComputeProgram = clCreateProgramWithSource(ComputeContext, 1, & source, NULL, &err);
  if (!ComputeProgram || err != CL_SUCCESS)
  {
    printf("%s\n", source);
    printf("Error: Failed to create compute program!\n");
    return ;
  }

  // Build the program executable
  //
  err = clBuildProgram(ComputeProgram, 0, NULL, NULL, NULL, NULL);
  if (err != CL_SUCCESS)
  {
    size_t length;
    char build_log[2048];
    printf("%s\n", source);
    printf("Error: Failed to build program executable!\n");
    clGetProgramBuildInfo(ComputeProgram, ComputeDeviceId, CL_PROGRAM_BUILD_LOG, sizeof(build_log), build_log, &length);
    printf("%s\n", build_log);
    return ;
  }

  ComputeKernels = (cl_kernel*) malloc(KernelCount * sizeof(cl_kernel));
  for(int i = 0; i < KernelCount; i++)
  {
    // Create each compute kernel from within the program
    //
    ComputeKernels[i] = clCreateKernel(ComputeProgram, KernelNames[i], &err);
    if (!ComputeKernels[i] || err != CL_SUCCESS)
    {
      printf("Error: Failed to create compute kernel!\n");
      return ;
    }

    size_t wgSize;
    err = clGetKernelWorkGroupInfo(ComputeKernels[i], ComputeDeviceId, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wgSize, NULL);
    if(err)
    {
      printf("Error: Failed to get kernel work group size\n");
    }
    GROUP_SIZE = std::min( GROUP_SIZE, wgSize );
  }

}
////////////////////////////////////////////////////////////////////////////////////////////////////

void PrefixScan::cleanup()
{
  for(int i = 0; i < KernelCount; i++)
    clReleaseKernel(ComputeKernels[i]);
  clReleaseProgram(ComputeProgram);

  free(ComputeKernels);
}


const char* PrefixScan::KernelNames[] =
{
  "PreScanKernel",
  "PreScanStoreSumKernel",
  "PreScanStoreSumNonPowerOfTwoKernel",
  "PreScanNonPowerOfTwoKernel",
  "UniformAddKernel"
};

const unsigned int PrefixScan::KernelCount
    = sizeof(PrefixScan::KernelNames) / sizeof(char *);

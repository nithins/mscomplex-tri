
#include <CL/cl.h>


typedef unsigned int SCAN_T;

class PrefixScan
{

public:

  PrefixScan();

  cl_program              ComputeProgram;
  cl_kernel*              ComputeKernels;

  cl_mem*                 ScanPartialSums;
  unsigned int            ElementsAllocated;
  unsigned int            LevelsAllocated;

  enum KernelMethods
  {
    PRESCAN                             = 0,
    PRESCAN_STORE_SUM                   = 1,
    PRESCAN_STORE_SUM_NON_POWER_OF_TWO  = 2,
    PRESCAN_NON_POWER_OF_TWO            = 3,
    UNIFORM_ADD                         = 4,
  };

  static const char* KernelNames[];

  static const unsigned int KernelCount;

  char * LoadProgramSourceFromFile(const char *filename);

  int CreatePartialSumBuffers(unsigned int count,cl_context& ComputeContext);

  void ReleasePartialSums(void);

  ////////////////////////////////////////////////////////////////////////////////////////////////////

  int PreScan
      (size_t *global,
       size_t *local,
       size_t shared,
       cl_mem output_data,
       cl_mem input_data,
       unsigned int n,
       int group_index,
       int base_index,
       cl_command_queue &ComputeCommands);

  int PreScanStoreSum
      (size_t *global,
       size_t *local,
       size_t shared,
       cl_mem output_data,
       cl_mem input_data,
       cl_mem partial_sums,
       unsigned int n,
       int group_index,
       int base_index,
       cl_command_queue &ComputeCommands);

  int PreScanStoreSumNonPowerOfTwo
      (size_t *global,
       size_t *local,
       size_t shared,
       cl_mem output_data,
       cl_mem input_data,
       cl_mem partial_sums,
       unsigned int n,
       int group_index,
       int base_index,
       cl_command_queue &ComputeCommands);

  int PreScanNonPowerOfTwo
      (size_t *global,
       size_t *local,
       size_t shared,
       cl_mem output_data,
       cl_mem input_data,
       unsigned int n,
       int group_index,
       int base_index,
       cl_command_queue &ComputeCommands);

  int UniformAdd
      (size_t *global,
       size_t *local,
       cl_mem output_data,
       cl_mem partial_sums,
       unsigned int n,
       unsigned int group_offset,
       unsigned int base_index,
       cl_command_queue &ComputeCommands);

  int PreScanBufferRecursive
      (cl_mem output_data,
       cl_mem input_data,
       int element_count,
       int level,
       cl_command_queue &ComputeCommands);

  void PreScanBuffer
      (cl_mem output_data,
       cl_mem input_data,
       unsigned int element_count,
       cl_command_queue &ComputeCommands);

  void init
      (cl_context &ComputeContext,
       cl_device_id &ComputeDeviceId,
       int groupsize = 256);

  void cleanup();
};

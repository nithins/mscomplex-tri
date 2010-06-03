/*
 * Copyright 1993-2009 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */

#ifndef BITONIC_SORT_H_INCLUDED
#define BITONIC_SORT_H_INCLUDED

#include <CL/cl.h>

class BitonicSortProgram
{
  typedef unsigned int uint;

public:
  void initBitonicSort(cl_context cxGPUContext, cl_device_id cdDevice);

  void closeBitonicSort();

  size_t bitonicSort(
      cl_command_queue cqCommandQueue,
      cl_mem d_DstKey,
      cl_mem d_DstVal,
      cl_mem d_SrcKey,
      cl_mem d_SrcVal,
      uint batch,
      uint arrayLength,
      uint dir
  );

private:
  //OpenCL bitonic sort program
  cl_program cpBitonicSort;

  //OpenCL bitonic sort kernels
  cl_kernel ckBitonicSortLocal,ckBitonicSortLocal1;
  cl_kernel ckBitonicMergeGlobal,ckBitonicMergeLocal;

};

#endif

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


#include <bitonic_sort.h>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <fstream>
#include <iostream>
#include <stdexcept>

#define _CHECKCL_ERR_CODE(_ERROR,_MESSAGE)\
if(_ERROR != CL_SUCCESS) throw std::runtime_error(_MESSAGE);

#define _GET_GLOBAL(s,l) (((s)/(2*(l)))*(2*(l)) + ((((s)%(2*(l))) == 0)?(0):(2*l)))

using namespace std;

void compile_cl_program(std::string prog_filename,std::string header_filename,
                        std::string compile_flags,cl_program &prog,cl_context & context,cl_device_id &device_id);


static const uint LOCAL_SIZE_LIMIT = 32U;
static const char  *compileOptions = "-D LOCAL_SIZE_LIMIT=32";

uint nextPow2(uint x)
{
  if (x == 0)
    return 0;

  x--;
  x = (x >> 1) | x;
  x = (x >> 2) | x;
  x = (x >> 4) | x;
  x = (x >> 8) | x;
  x = (x >> 16) | x;
  x++;
  return x;
};


void BitonicSortProgram::initBitonicSort(cl_context cxGPUContext,cl_device_id cdDevice){
  cl_int ciErrNum;

  compile_cl_program(":/oclsources/bitonic_sort.cl","",compileOptions,cpBitonicSort,cxGPUContext,cdDevice);


  ckBitonicSortLocal = clCreateKernel(cpBitonicSort, "bitonicSortLocal", &ciErrNum);
  _CHECKCL_ERR_CODE(ciErrNum,"failed to create bitonicSortLocal Kernel");

  ckBitonicSortLocal1 = clCreateKernel(cpBitonicSort, "bitonicSortLocal1", &ciErrNum);
  _CHECKCL_ERR_CODE(ciErrNum,"failed to create bitonicSortLocal1 Kernel");

  ckBitonicMergeGlobal = clCreateKernel(cpBitonicSort, "bitonicMergeGlobal", &ciErrNum);
  _CHECKCL_ERR_CODE(ciErrNum,"failed to create bitonicMergeGlobal Kernel");

  ckBitonicMergeLocal = clCreateKernel(cpBitonicSort, "bitonicMergeLocal", &ciErrNum);
  _CHECKCL_ERR_CODE(ciErrNum,"failed to create bitonicMergeLocal Kernel");


  //Check for work group size
  size_t szBitonicSortLocal, szBitonicSortLocal1, szBitonicMergeLocal;

  ciErrNum |= clGetKernelWorkGroupInfo
              (ckBitonicSortLocal,cdDevice, CL_KERNEL_WORK_GROUP_SIZE,
               sizeof(size_t),&szBitonicSortLocal, NULL);

  ciErrNum |= clGetKernelWorkGroupInfo
              (ckBitonicSortLocal1, cdDevice, CL_KERNEL_WORK_GROUP_SIZE,
               sizeof(size_t), &szBitonicSortLocal1, NULL);

  ciErrNum |= clGetKernelWorkGroupInfo
              (ckBitonicMergeLocal, cdDevice, CL_KERNEL_WORK_GROUP_SIZE,
               sizeof(size_t), &szBitonicMergeLocal, NULL);

  _CHECKCL_ERR_CODE(ciErrNum,"failed to get work group info");

  if( (szBitonicSortLocal < (LOCAL_SIZE_LIMIT / 2)) ||
      (szBitonicSortLocal1 < (LOCAL_SIZE_LIMIT / 2)) ||
      (szBitonicMergeLocal < (LOCAL_SIZE_LIMIT / 2)) )
  {
    throw std::runtime_error("Minimum work-group size required by this "\
                             "application is not supported on this device");
  }
}

void BitonicSortProgram::closeBitonicSort(void)
{
  cl_int ciErrNum;
  ciErrNum  = clReleaseKernel(ckBitonicMergeLocal);
  ciErrNum |= clReleaseKernel(ckBitonicMergeGlobal);
  ciErrNum |= clReleaseKernel(ckBitonicSortLocal1);
  ciErrNum |= clReleaseKernel(ckBitonicSortLocal);
  ciErrNum |= clReleaseProgram(cpBitonicSort);
  _CHECKCL_ERR_CODE(ciErrNum,"Failed to release kernels");
}

size_t BitonicSortProgram::bitonicSort(
    cl_command_queue cqCommandQueue,
    cl_mem d_DstKey,
    cl_mem d_DstVal,
    cl_mem d_SrcKey,
    cl_mem d_SrcVal,
    uint batch,
    uint arrayLength,
    uint dir
    )
{
  if(arrayLength < 2)
    return 0;

  dir = (dir != 0);

  cl_int ciErrNum;
  size_t localWorkSize;
  size_t globalWorkSize;

  if(arrayLength <= LOCAL_SIZE_LIMIT)
  {
    // oclCheckError( (batch * arrayLength) % LOCAL_SIZE_LIMIT == 0, shrTRUE );
    //Launch bitonicSortLocal
    ciErrNum  = clSetKernelArg(ckBitonicSortLocal, 0,   sizeof(cl_mem), (void *)&d_DstKey);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal, 1,   sizeof(cl_mem), (void *)&d_DstVal);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal, 2,   sizeof(cl_mem), (void *)&d_SrcKey);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal, 3,   sizeof(cl_mem), (void *)&d_SrcVal);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal, 4,  sizeof(cl_uint), (void *)&arrayLength);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal, 5,  sizeof(cl_uint), (void *)&dir);

    _CHECKCL_ERR_CODE(ciErrNum,"failed to set args for ckBitonicSortLocal");

    localWorkSize  = LOCAL_SIZE_LIMIT / 2;
    globalWorkSize = _GET_GLOBAL(batch * arrayLength / 2,localWorkSize);

    ciErrNum = clEnqueueNDRangeKernel
               (cqCommandQueue, ckBitonicSortLocal, 1, NULL,
                &globalWorkSize, &localWorkSize, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(ciErrNum,"failed to enqueue ckBitonicSortLocal");
  }
  else
  {
    //Launch bitonicSortLocal1
    ciErrNum  = clSetKernelArg(ckBitonicSortLocal1, 0,  sizeof(cl_mem), (void *)&d_DstKey);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal1, 1,  sizeof(cl_mem), (void *)&d_DstVal);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal1, 2,  sizeof(cl_mem), (void *)&d_SrcKey);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal1, 3,  sizeof(cl_mem), (void *)&d_SrcVal);
    ciErrNum |= clSetKernelArg(ckBitonicSortLocal1, 4,  sizeof(cl_uint), (void *)&arrayLength);

    _CHECKCL_ERR_CODE(ciErrNum,"failed to set args for ckBitonicSortLocal1");

    localWorkSize = LOCAL_SIZE_LIMIT / 2;
    globalWorkSize = _GET_GLOBAL(batch * arrayLength / 2,localWorkSize);
    ciErrNum = clEnqueueNDRangeKernel
               (cqCommandQueue, ckBitonicSortLocal1, 1, NULL,
                &globalWorkSize, &localWorkSize, 0, NULL, NULL);

    _CHECKCL_ERR_CODE(ciErrNum,"failed to set enqueue ckBitonicSortLocal1");

    uint arrayLength_nextpow2 =  nextPow2( arrayLength);

    for(uint size = 2 * LOCAL_SIZE_LIMIT; size <= arrayLength_nextpow2; size <<= 1)
    {
      printf("size = %u\n",size);

      for(unsigned stride = size / 2; stride > 0; stride >>= 1)
      {
        printf("stride = %u\n",stride);

        if(stride >= LOCAL_SIZE_LIMIT)
        {
          //Launch bitonicMergeGlobal
          ciErrNum  = clSetKernelArg(ckBitonicMergeGlobal, 0,  sizeof(cl_mem), (void *)&d_DstKey);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 1,  sizeof(cl_mem), (void *)&d_DstVal);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 2,  sizeof(cl_mem), (void *)&d_DstKey);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 3,  sizeof(cl_mem), (void *)&d_DstVal);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 4, sizeof(cl_uint), (void *)&arrayLength);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 5, sizeof(cl_uint), (void *)&size);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 6, sizeof(cl_uint), (void *)&stride);
          ciErrNum |= clSetKernelArg(ckBitonicMergeGlobal, 7, sizeof(cl_uint), (void *)&dir);
          _CHECKCL_ERR_CODE(ciErrNum,"failed to set args ckBitonicMergeGlobal");

          globalWorkSize = _GET_GLOBAL(batch * arrayLength / 2,localWorkSize);
          ciErrNum = clEnqueueNDRangeKernel
                     (cqCommandQueue, ckBitonicMergeGlobal, 1, NULL,
                      &globalWorkSize, NULL, 0, NULL, NULL);

          _CHECKCL_ERR_CODE(ciErrNum,"failed to enqueue ckBitonicMergeGlobal");
        }
        else
        {
          //          //Launch bitonicMergeLocal
          ciErrNum  = clSetKernelArg(ckBitonicMergeLocal, 0,  sizeof(cl_mem), (void *)&d_DstKey);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 1,  sizeof(cl_mem), (void *)&d_DstVal);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 2,  sizeof(cl_mem), (void *)&d_DstKey);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 3,  sizeof(cl_mem), (void *)&d_DstVal);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 4, sizeof(cl_uint), (void *)&arrayLength);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 5, sizeof(cl_uint), (void *)&stride);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 6, sizeof(cl_uint), (void *)&size);
          ciErrNum |= clSetKernelArg(ckBitonicMergeLocal, 7, sizeof(cl_uint), (void *)&dir);

          _CHECKCL_ERR_CODE(ciErrNum,"failed to set args bitonicMergeLocal");

          localWorkSize  = LOCAL_SIZE_LIMIT / 2;
          globalWorkSize = _GET_GLOBAL(batch * arrayLength / 2,localWorkSize);

          ciErrNum = clEnqueueNDRangeKernel(cqCommandQueue, ckBitonicMergeLocal, 1, NULL, &globalWorkSize, &localWorkSize, 0, NULL, NULL);
          _CHECKCL_ERR_CODE(ciErrNum,"failed to enqueue bitonicMergeLocal");
          break;
        }
      }
    }
  }

  return localWorkSize;
}

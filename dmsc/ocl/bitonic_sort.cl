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



//Passed down by clBuildProgram
//#define LOCAL_SIZE_LIMIT 1024


#define TYPE uint

inline int compare(TYPE a,TYPE b)
{
  return  (a > b)?(1):(0);
}

inline void ComparatorPrivate(
    TYPE *keyA,
    TYPE *valA,
    TYPE *keyB,
    TYPE *valB,
    uint arrowDir
)
{
    if( compare(*keyA,*keyB) == arrowDir )
    {
        TYPE t;
        t = *keyA; *keyA = *keyB; *keyB = t;
        t = *valA; *valA = *valB; *valB = t;
    }
}

inline void ComparatorLocal(
    __local TYPE *keyA,
    __local TYPE *valA,
    __local TYPE *keyB,
    __local TYPE *valB,
    uint arrowDir
){
    if( compare(*keyA,*keyB) == arrowDir )
    {
        TYPE t;
        t = *keyA; *keyA = *keyB; *keyB = t;
        t = *valA; *valA = *valB; *valB = t;
    }
}

inline uint nextPow2(uint x)
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
}

////////////////////////////////////////////////////////////////////////////////
// Monolithic bitonic sort kernel for short arrays fitting into local memory
////////////////////////////////////////////////////////////////////////////////
__kernel __attribute__((reqd_work_group_size(LOCAL_SIZE_LIMIT / 2, 1, 1)))
void bitonicSortLocal(
    __global TYPE *d_DstKey,
    __global TYPE *d_DstVal,
    __global TYPE *d_SrcKey,
    __global TYPE *d_SrcVal,
    uint arrayLength,
    uint sortDir
){
    __local  TYPE l_key[LOCAL_SIZE_LIMIT];
    __local  TYPE l_val[LOCAL_SIZE_LIMIT];

    uint idx = get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    //Offset to the beginning of subarray and load data
    d_SrcKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_SrcVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    if(idx < arrayLength)
    {
      l_key[get_local_id(0) +                      0] = d_SrcKey[                     0];
      l_val[get_local_id(0) +                      0] = d_SrcVal[                     0];
    }
    if(idx + (LOCAL_SIZE_LIMIT / 2) < arrayLength)
    {
      l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcKey[(LOCAL_SIZE_LIMIT / 2)];
      l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcVal[(LOCAL_SIZE_LIMIT / 2)];
    }

    uint arrayLength_nextPow2 = nextPow2(arrayLength);

    for(uint size = 2; size < arrayLength_nextPow2; size <<= 1){
        //Bitonic merge
        uint dir = ( (get_local_id(0) & (size / 2)) != 0 );
        for(uint stride = size / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            if(pos + stride + get_group_id(0) * LOCAL_SIZE_LIMIT < arrayLength)
            {
              ComparatorLocal(
                  &l_key[pos +      0], &l_val[pos +      0],
                  &l_key[pos + stride], &l_val[pos + stride],
                  dir
              );
            }
        }
    }

    //dir == sortDir for the last bitonic merge step
    {
        for(uint stride = arrayLength_nextPow2 / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            if(pos + stride + get_group_id(0) * LOCAL_SIZE_LIMIT < arrayLength)
            {
              ComparatorLocal(
                  &l_key[pos +      0], &l_val[pos +      0],
                  &l_key[pos + stride], &l_val[pos + stride],
                  sortDir
              );
            }
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if(idx < arrayLength)
    {
      d_DstKey[                     0] = l_key[get_local_id(0) +                      0];
      d_DstVal[                     0] = l_val[get_local_id(0) +                      0];
    }
    if(idx + LOCAL_SIZE_LIMIT / 2< arrayLength)
    {
      d_DstKey[(LOCAL_SIZE_LIMIT / 2)] = l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
      d_DstVal[(LOCAL_SIZE_LIMIT / 2)] = l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
    }

}

////////////////////////////////////////////////////////////////////////////////
// Bitonic sort kernel for large arrays (not fitting into local memory)
////////////////////////////////////////////////////////////////////////////////
//Bottom-level bitonic sort
//Almost the same as bitonicSortLocal with the only exception
//of even / odd subarrays (of LOCAL_SIZE_LIMIT points) being
//sorted in opposite directions
__kernel __attribute__((reqd_work_group_size(LOCAL_SIZE_LIMIT / 2, 1, 1)))
void bitonicSortLocal1(
    __global TYPE *d_DstKey,
    __global TYPE *d_DstVal,
    __global TYPE *d_SrcKey,
    __global TYPE *d_SrcVal,
    uint arrayLength
){
    __local TYPE l_key[LOCAL_SIZE_LIMIT];
    __local TYPE l_val[LOCAL_SIZE_LIMIT];

    uint idx = get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    //Offset to the beginning of subarray and load data
    d_SrcKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_SrcVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    if(idx < arrayLength)
    {
      l_key[get_local_id(0) +                      0] = d_SrcKey[                     0];
      l_val[get_local_id(0) +                      0] = d_SrcVal[                     0];
    }
    if(idx + (LOCAL_SIZE_LIMIT / 2) < arrayLength)
    {
      l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcKey[(LOCAL_SIZE_LIMIT / 2)];
      l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcVal[(LOCAL_SIZE_LIMIT / 2)];
    }

    uint comparatorI = get_global_id(0) & ((LOCAL_SIZE_LIMIT / 2) - 1);

    for(uint size = 2; size < LOCAL_SIZE_LIMIT; size <<= 1){
        //Bitonic merge
        uint dir = (comparatorI & (size / 2)) != 0;
        for(uint stride = size / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));

            if(pos + stride + get_group_id(0) * LOCAL_SIZE_LIMIT < arrayLength)
            {
              ComparatorLocal(
                  &l_key[pos +      0], &l_val[pos +      0],
                  &l_key[pos + stride], &l_val[pos + stride],
                  dir
              );
            }
        }
    }
//
    //Odd / even arrays of LOCAL_SIZE_LIMIT elements
    //sorted in opposite directions
    {
        uint dir = (get_group_id(0) & 1);
        for(uint stride = LOCAL_SIZE_LIMIT / 2; stride > 0; stride >>= 1){
            barrier(CLK_LOCAL_MEM_FENCE);
            uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));
            if(pos + stride + get_group_id(0) * LOCAL_SIZE_LIMIT < arrayLength)
            {
              ComparatorLocal(
                  &l_key[pos +      0], &l_val[pos +      0],
                  &l_key[pos + stride], &l_val[pos + stride],
                dir
              );
            }
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if(idx < arrayLength)
    {
      d_DstKey[                     0] = l_key[get_local_id(0) +                      0];
      d_DstVal[                     0] = l_val[get_local_id(0) +                      0];
    }
    if(idx + LOCAL_SIZE_LIMIT / 2< arrayLength)
    {
      d_DstKey[(LOCAL_SIZE_LIMIT / 2)] = l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
      d_DstVal[(LOCAL_SIZE_LIMIT / 2)] = l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
    }
}

//Bitonic merge iteration for 'stride' >= LOCAL_SIZE_LIMIT
__kernel void bitonicMergeGlobal(
    __global TYPE *d_DstKey,
    __global TYPE *d_DstVal,
    __global TYPE *d_SrcKey,
    __global TYPE *d_SrcVal,
    uint arrayLength,
    uint size,
    uint stride,
    uint sortDir
){
    uint global_comparatorI = get_global_id(0);
    uint        comparatorI = global_comparatorI & (arrayLength / 2 - 1);

    //Bitonic merge
    uint dir = sortDir ^ ( (comparatorI & (size / 2)) != 0 );
    uint pos = 2 * global_comparatorI - (global_comparatorI & (stride - 1));

    if(pos + stride < arrayLength)
    {
      TYPE keyA = d_SrcKey[pos +      0];
      TYPE valA = d_SrcVal[pos +      0];
      TYPE keyB = d_SrcKey[pos + stride];
      TYPE valB = d_SrcVal[pos + stride];

      ComparatorPrivate(
          &keyA, &valA,
          &keyB, &valB,
          dir
      );

      d_DstKey[pos +      0] = keyA;
      d_DstVal[pos +      0] = valA;
      d_DstKey[pos + stride] = keyB;
      d_DstVal[pos + stride] = valB;
    }
}

//Combined bitonic merge steps for
//'size' > LOCAL_SIZE_LIMIT and 'stride' = [1 .. LOCAL_SIZE_LIMIT / 2]
__kernel __attribute__((reqd_work_group_size(LOCAL_SIZE_LIMIT / 2, 1, 1)))
void bitonicMergeLocal(
    __global TYPE *d_DstKey,
    __global TYPE *d_DstVal,
    __global TYPE *d_SrcKey,
    __global TYPE *d_SrcVal,
    uint arrayLength,
    uint stride,
    uint size,
    uint sortDir
){
    __local TYPE l_key[LOCAL_SIZE_LIMIT];
    __local TYPE l_val[LOCAL_SIZE_LIMIT];

    uint idx = get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    d_SrcKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_SrcVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstKey += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);
    d_DstVal += get_group_id(0) * LOCAL_SIZE_LIMIT + get_local_id(0);

    if(idx < arrayLength)
    {
      l_key[get_local_id(0) +                      0] = d_SrcKey[                     0];
      l_val[get_local_id(0) +                      0] = d_SrcVal[                     0];
    }

    if(idx + (LOCAL_SIZE_LIMIT / 2) < arrayLength)
    {
      l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcKey[(LOCAL_SIZE_LIMIT / 2)];
      l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)] = d_SrcVal[(LOCAL_SIZE_LIMIT / 2)];
    }

    //Bitonic merge
    uint comparatorI = get_global_id(0) & ((arrayLength / 2) - 1);
    uint         dir = sortDir ^ ( (comparatorI & (size / 2)) != 0 );
    for(; stride > 0; stride >>= 1)
    {
        barrier(CLK_LOCAL_MEM_FENCE);
        uint pos = 2 * get_local_id(0) - (get_local_id(0) & (stride - 1));

        if(pos + stride + get_group_id(0) * LOCAL_SIZE_LIMIT < arrayLength)
        {
          ComparatorLocal(
              &l_key[pos +      0], &l_val[pos +      0],
              &l_key[pos + stride], &l_val[pos + stride],
              dir
          );
        }
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    if(idx < arrayLength)
    {
      d_DstKey[                     0] = l_key[get_local_id(0) +                      0];
      d_DstVal[                     0] = l_val[get_local_id(0) +                      0];
    }
    if(idx + (LOCAL_SIZE_LIMIT / 2) < arrayLength)
    {
      d_DstKey[(LOCAL_SIZE_LIMIT / 2)] = l_key[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
      d_DstVal[(LOCAL_SIZE_LIMIT / 2)] = l_val[get_local_id(0) + (LOCAL_SIZE_LIMIT / 2)];
    }
}

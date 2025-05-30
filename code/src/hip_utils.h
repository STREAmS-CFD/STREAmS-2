#include <hip/hip_runtime.h>
#include "hip/hip_runtime_api.h"
#include "hipcub/hipcub.hpp"
#include <float.h>
#include <limits>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <vector>

//Precision
#ifdef SINGLE_PRECISION
  typedef float real;
#else
  typedef double real;
#endif

//Loop bounds template
template <typename I, typename E, typename S> __device__ __forceinline__ bool loop_cond(I idx,E end,S stride) {
  return idx <= end;     
}

// Thread specification for flux evaluation
#define EULERCENTRAL_THREADS_X 256
#define EULERCENTRAL_THREADS_Y 3
#define EULERWENO_THREADS_X 128
#define EULERWENO_THREADS_Y 2

// One loop parallel thread definition
#define ONE_X 128

// Two loop parallel thread definition
#define TWO_X 128
#define TWO_Y 2

// Three loop parallel thread definition
#define THREE_X 64
#define THREE_Y 4
#define THREE_Z 1

// Access the thread index and offset if required
#define __GIDX(idx,off) 1+(threadIdx.idx + blockIdx.idx * blockDim.idx)+((off)-1)

// Round up the grid dimensions
#define divideAndRoundUp(x, y) ((x) / (y) + ((x) % (y) != 0))

// Max function
#define SIGN(a,b) ( ( (b) < 0 ) ? (-(abs(a))) : (abs(a)) )

// HIP C++ version of HIP_CHECK used for reductions 
#define HIP_CHECK(condition)         \
  {                                  \
    hipError_t error = condition;    \
    if(error != hipSuccess){         \
        std::cout << "HIP error: " << hipGetErrorString(error) << " line: " << __LINE__ << std::endl; \
        exit(error); \
    } \
  }

#ifdef INCLUDE_KERNELS
#define ___I3_REDN_3D(i,j,k) (((i)-(1))+((nx)-(1)+1)*((j)-(1))+((nx)-(1)+1)*((ny)-(1)+1)*((k)-(1)))
__global__ void  reduce_init_kernel(int nx,int ny,int nz,real *redn_3d_gpu){
//Kernel for initialization of reductions
int i;int j;int k;
i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);
if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[___I3_REDN_3D(i,j,k)]=0.0;
}
}
#endif

// Hipcub reduction implementation
// Taken directly from GPUFORT
// https://github.com/ROCmSoftwarePlatform/gpufort/blob/main/include/gpufort_reductions.h
namespace {
  struct reduce_op_mult {
    template <typename T> 
    static __host__ __device__ __forceinline__ T ival() { return (T)1; }
    template <typename T> 
    static __host__ __device__ __forceinline__ void init(T &a) { a = ival<T>(); }
    template <typename T>
    __device__ __forceinline__ T operator()(const T &a, const T &b) const {
      return a * b;
    }
  };
  
  struct reduce_op_add {
    template <typename T> 
    static __host__ __device__ __forceinline__ T ival() { return (T)0; }
    template <typename T>
    static __host__ __device__ __forceinline__ void init(T &a) { a = ival<T>(); }
    template <typename T>
    __device__ __forceinline__ T operator()(const T &a, const T &b) const {
      return a + b;
    }
  };
  
  struct reduce_op_max {
    template <typename T> 
    static __host__ __device__ __forceinline__ T ival() {
      return -std::numeric_limits<T>::max(); // has negative sign
    }
    template <typename T> 
    static __host__ __device__ __forceinline__ void init(T &a) { a = ival<T>(); }
    template <typename T>
    __device__ __forceinline__ T operator()(const T &a, const T &b) const {
      return std::max(a, b);
    }
  };
  
  struct reduce_op_min {
    template <typename T>
    static __host__ __device__ __forceinline__ T ival() {
      return std::numeric_limits<T>::max();
    }
    template <typename T>
    static __host__ __device__ __forceinline__ void init(T &a) { a = ival<T>(); }
    template <typename T>
    __device__ __forceinline__ T operator()(const T &a, const T &b) const {
      return std::min(a, b);
    }
  };
  
  template <typename T, typename ReduceOpT>
  void reduce(const T *const d_in, const int &NX, const T *h_out) {
    T *d_out = nullptr;
    HIP_CHECK(hipMalloc((void **)&d_out, sizeof(T)));
    //Determine temporary device storage requirements
    void *temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    ReduceOpT reduceOp;
    hipcub::DeviceReduce::Reduce(temp_storage, temp_storage_bytes, d_in, d_out, NX,
                                 ReduceOpT(), ReduceOpT::template ival<T>());
    //Allocate temporary storage
    HIP_CHECK(hipMalloc(&temp_storage, temp_storage_bytes));
    //Run reduction
    hipcub::DeviceReduce::Reduce(temp_storage, temp_storage_bytes, d_in, d_out, NX,
                                 ReduceOpT(), ReduceOpT::template ival<T>());
    HIP_CHECK(hipMemcpy((void *)h_out, d_out, sizeof(T), hipMemcpyDeviceToHost));
    //Clean up
    HIP_CHECK(hipFree(d_out));
    HIP_CHECK(hipFree(temp_storage));
  }
}

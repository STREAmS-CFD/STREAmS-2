#include <hip/hip_runtime.h>
#include "hip_utils.h"
#include "amd_arrays.h"



__global__ void  bcte_step_1_kernel1(int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftus_gpu,real *wbuftes_gpu){
//Kernel for bcte_step_1_kernel1
int i;int j;int k;
int m;int iercuda;

k = __GIDX(x,1);


if(loop_cond(k,nz,1)){
wbuftes_gpu[__I2_WBUFTES(k,1)] = w_gpu[__I4_W(ite_l,1,k,1)];
if(nv > 1) wbuftes_gpu[__I2_WBUFTES(k,2)] = w_gpu[__I4_W(ite_l,1,k,5)];

}
}


extern "C"{
void bcte_step_1_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftus_gpu,real *wbuftes_gpu){
dim3 block(ONE_X);
dim3 grid(divideAndRoundUp((nz)-(1)+1,block.x));

hipLaunchKernelGGL((bcte_step_1_kernel1),grid,block,0,stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu);
}
}


__global__ void  bcte_step_1_kernel2(int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftus_gpu,real *wbuftes_gpu){
//Kernel for bcte_step_1_kernel2
int i;int j;int k;
int m;int iercuda;

k = __GIDX(x,1);


if(loop_cond(k,nz,1)){
wbuftus_gpu[__I2_WBUFTUS(k,1)] = w_gpu[__I4_W(itu_l,1,k,1)];
if(nv > 1) wbuftus_gpu[__I2_WBUFTUS(k,2)] = w_gpu[__I4_W(itu_l,1,k,5)];

}
}


extern "C"{
void bcte_step_1_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftus_gpu,real *wbuftes_gpu){
dim3 block(ONE_X);
dim3 grid(divideAndRoundUp((nz)-(1)+1,block.x));

hipLaunchKernelGGL((bcte_step_1_kernel2),grid,block,0,stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftus_gpu,wbuftes_gpu);
}
}



__global__ void  bcte_step_3_kernel1(int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftur_gpu,real *wbufter_gpu){
//Kernel for bcte_step_3_kernel1
int k;int iercuda;

k = __GIDX(x,1);


if(loop_cond(k,nz,1)){
w_gpu[__I4_W(ite_l,1,k,1)] = 0.50*(w_gpu[__I4_W(ite_l,1,k,1)]+wbuftur_gpu[__I2_WBUFTUR(k,1)]);
if(nv > 1) w_gpu[__I4_W(ite_l,1,k,5)] = 0.50*(w_gpu[__I4_W(ite_l,1,k,5)]+wbuftur_gpu[__I2_WBUFTUR(k,2)]);

}
}


extern "C"{
void bcte_step_3_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftur_gpu,real *wbufter_gpu){
dim3 block(ONE_X);
dim3 grid(divideAndRoundUp((nz)-(1)+1,block.x));

hipLaunchKernelGGL((bcte_step_3_kernel1),grid,block,0,stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu);
}
}


__global__ void  bcte_step_3_kernel2(int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftur_gpu,real *wbufter_gpu){
//Kernel for bcte_step_3_kernel2
int k;int iercuda;

k = __GIDX(x,1);


if(loop_cond(k,nz,1)){
w_gpu[__I4_W(itu_l,1,k,1)] = 0.50*(w_gpu[__I4_W(itu_l,1,k,1)]+wbufter_gpu[__I2_WBUFTER(k,1)]);
if(nv > 1) w_gpu[__I4_W(itu_l,1,k,5)] = 0.50*(w_gpu[__I4_W(itu_l,1,k,5)]+wbufter_gpu[__I2_WBUFTER(k,2)]);

}
}


extern "C"{
void bcte_step_3_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int nrank_x,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,real *w_gpu,real *wbuftur_gpu,real *wbufter_gpu){
dim3 block(ONE_X);
dim3 grid(divideAndRoundUp((nz)-(1)+1,block.x));

hipLaunchKernelGGL((bcte_step_3_kernel2),grid,block,0,stream,nx,ny,nz,nv,ng,nrank_x,ite_rank_x,itu_rank_x,ite_l,itu_l,w_gpu,wbuftur_gpu,wbufter_gpu);
}
}



__global__ void  bcswap_wake_step_1_kernel(int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf4s_gpu){
//Kernel for bcswap_wake_step_1_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf4s_gpu[__I4_WBUF4S(i,j,k,m)] = w_gpu[__I4_W(nx-i+1,1+j,k,m)];
}

}
}


extern "C"{
void bcswap_wake_step_1_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf4s_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_wake_step_1_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu,wbuf4s_gpu);
}
}



__global__ void  bcswap_wake_step_3_kernel(int nx,int ny,int nz,int ng,int nv,int *wall_tag_gpu,real *w_gpu,real *wbuf3r_gpu){
//Kernel for bcswap_wake_step_3_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
if (wall_tag_gpu[__I1_WALL_TAG(i)] > 0) {
w_gpu[__I4_W(i,1-j,k,m)] = wbuf3r_gpu[__I4_WBUF3R(i,j,k,m)];
}
}

}
}


extern "C"{
void bcswap_wake_step_3_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int *wall_tag_gpu,real *w_gpu,real *wbuf3r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_wake_step_3_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,wall_tag_gpu,w_gpu,wbuf3r_gpu);
}
}



__global__ void  bcswap_step_1_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
//Kernel for bcswap_step_1_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf1s_gpu[__I4_WBUF1S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,m)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_1_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_1_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu);
}
}


__global__ void  bcswap_step_1_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
//Kernel for bcswap_step_1_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf3s_gpu[__I4_WBUF3S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf4s_gpu[__I4_WBUF4S(i,j,k,m)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
}

}
}


extern "C"{
void bcswap_step_1_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_1_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu);
}
}


__global__ void  bcswap_step_1_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
//Kernel for bcswap_step_1_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
wbuf5s_gpu[__I4_WBUF5S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf6s_gpu[__I4_WBUF6S(i,j,k,m)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
}

}
}


extern "C"{
void bcswap_step_1_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_1_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu);
}
}



__global__ void  bcswap_c2_step_1_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2_step_1_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
wbuf1s_gpu[__I4_WBUF1S(i,j,k,1)] = w_gpu[__I4_W(i,j,k,1)];
wbuf1s_gpu[__I4_WBUF1S(i,j,k,2)] = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
wbuf1s_gpu[__I4_WBUF1S(i,j,k,3)] = w_gpu[__I4_W(i,j,k,2)]*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+w_gpu[__I4_W(i,j,k,3)]*detadync2_gpu[__I2_DETADYNC2(i,j)];
wbuf1s_gpu[__I4_WBUF1S(i,j,k,4)] = w_gpu[__I4_W(i,j,k,4)];
wbuf1s_gpu[__I4_WBUF1S(i,j,k,5)] = w_gpu[__I4_W(i,j,k,5)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,1)] = w_gpu[__I4_W(nx-ng+i,j,k,1)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,2)] = w_gpu[__I4_W(nx-ng+i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(nx-ng+i,j)]+w_gpu[__I4_W(nx-ng+i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(nx-ng+i,j)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,3)] = w_gpu[__I4_W(nx-ng+i,j,k,2)]*detadxnc2_gpu[__I2_DETADXNC2(nx-ng+i,j)]+w_gpu[__I4_W(nx-ng+i,j,k,3)]*detadync2_gpu[__I2_DETADYNC2(nx-ng+i,j)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,4)] = w_gpu[__I4_W(nx-ng+i,j,k,4)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,5)] = w_gpu[__I4_W(nx-ng+i,j,k,5)];

}
}


extern "C"{
void bcswap_c2_step_1_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_1_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2_step_1_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2_step_1_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf1s_gpu[__I4_WBUF1S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf2s_gpu[__I4_WBUF2S(i,j,k,m)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_1_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_1_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2_step_1_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2_step_1_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf3s_gpu[__I4_WBUF3S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf4s_gpu[__I4_WBUF4S(i,j,k,m)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_1_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_1_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2_step_1_kernel4(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2_step_1_kernel4
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
wbuf5s_gpu[__I4_WBUF5S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf6s_gpu[__I4_WBUF6S(i,j,k,m)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_1_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1s_gpu,real *wbuf2s_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_1_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1s_gpu,wbuf2s_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}



__global__ void  bcswap_c2xybc_step_1_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2xybc_step_1_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,1)] = w_gpu[__I4_W(i,j,k,1)];
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,2)] = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,3)] = w_gpu[__I4_W(i,j,k,2)]*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+w_gpu[__I4_W(i,j,k,3)]*detadync2_gpu[__I2_DETADYNC2(i,j)];
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,4)] = w_gpu[__I4_W(i,j,k,4)];
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,5)] = w_gpu[__I4_W(i,j,k,5)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,1)] = w_gpu[__I4_W(nx-ng+i,j,k,1)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,2)] = w_gpu[__I4_W(nx-ng+i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(nx-ng+i,j)]+w_gpu[__I4_W(nx-ng+i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(nx-ng+i,j)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,3)] = w_gpu[__I4_W(nx-ng+i,j,k,2)]*detadxnc2_gpu[__I2_DETADXNC2(nx-ng+i,j)]+w_gpu[__I4_W(nx-ng+i,j,k,3)]*detadync2_gpu[__I2_DETADYNC2(nx-ng+i,j)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,4)] = w_gpu[__I4_W(nx-ng+i,j,k,4)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,5)] = w_gpu[__I4_W(nx-ng+i,j,k,5)];

}
}


extern "C"{
void bcswap_c2xybc_step_1_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_1_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_1_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2xybc_step_1_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf1xybcs_gpu[__I4_WBUF1XYBCS(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf2xybcs_gpu[__I4_WBUF2XYBCS(i,j,k,m)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_1_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_1_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_1_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2xybc_step_1_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf3s_gpu[__I4_WBUF3S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf4s_gpu[__I4_WBUF4S(i,j,k,m)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_1_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_1_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_1_kernel4(int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
//Kernel for bcswap_c2xybc_step_1_kernel4
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
wbuf5s_gpu[__I4_WBUF5S(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)];
wbuf6s_gpu[__I4_WBUF6S(i,j,k,m)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_1_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,real *w_gpu,real *wbuf1xybcs_gpu,real *wbuf2xybcs_gpu,real *wbuf3s_gpu,real *wbuf4s_gpu,real *wbuf5s_gpu,real *wbuf6s_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_1_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,w_gpu,wbuf1xybcs_gpu,wbuf2xybcs_gpu,wbuf3s_gpu,wbuf4s_gpu,wbuf5s_gpu,wbuf6s_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu);
}
}



__global__ void  bcswap_corner_step_1_kernel(int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf1s_c_gpu,real *wbuf2s_c_gpu,real *wbuf3s_c_gpu,real *wbuf4s_c_gpu){
//Kernel for bcswap_corner_step_1_kernel
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(j,ny,1)){
for(int k=1; k<ng+1; k++){
for(int i=1; i<ng+1; i++){
wbuf1s_c_gpu[__I4_WBUF1S_C(m,i,j,k)] = w_gpu[__I4_W(i,j,k,m)];
wbuf2s_c_gpu[__I4_WBUF2S_C(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,j,nz-ng+k,m)];
wbuf3s_c_gpu[__I4_WBUF3S_C(m,i,j,k)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
wbuf4s_c_gpu[__I4_WBUF4S_C(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
}
}

}
}


extern "C"{
void bcswap_corner_step_1_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf1s_c_gpu,real *wbuf2s_c_gpu,real *wbuf3s_c_gpu,real *wbuf4s_c_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_1_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu,wbuf1s_c_gpu,wbuf2s_c_gpu,wbuf3s_c_gpu,wbuf4s_c_gpu);
}
}



__global__ void  bcswap_corner_step_1b_kernel1(int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
//Kernel for bcswap_corner_step_1b_kernel1
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,nz,1)){
for(int j=1; j<ng+1; j++){
for(int i=1; i<ng+1; i++){
wbuf_lxly_s_gpu[__I4_WBUF_LXLY_S(m,i,j,k)] = w_gpu[__I4_W(i,j,k,m)];
wbuf_lxry_s_gpu[__I4_WBUF_LXRY_S(m,i,j,k)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
wbuf_rxly_s_gpu[__I4_WBUF_RXLY_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
wbuf_rxry_s_gpu[__I4_WBUF_RXRY_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,ny-ng+j,k,m)];
}
}

}
}


extern "C"{
void bcswap_corner_step_1b_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_1b_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu);
}
}


__global__ void  bcswap_corner_step_1b_kernel2(int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
//Kernel for bcswap_corner_step_1b_kernel2
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
i = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(i,nx,1)){
for(int k=1; k<ng+1; k++){
for(int j=1; j<ng+1; j++){
wbuf_lylz_s_gpu[__I4_WBUF_LYLZ_S(m,i,j,k)] = w_gpu[__I4_W(i,j,k,m)];
wbuf_lyrz_s_gpu[__I4_WBUF_LYRZ_S(m,i,j,k)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
wbuf_rylz_s_gpu[__I4_WBUF_RYLZ_S(m,i,j,k)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
wbuf_ryrz_s_gpu[__I4_WBUF_RYRZ_S(m,i,j,k)] = w_gpu[__I4_W(i,ny-ng+j,nz-ng+k,m)];
}
}

}
}


extern "C"{
void bcswap_corner_step_1b_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_1b_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu);
}
}


__global__ void  bcswap_corner_step_1b_kernel3(int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
//Kernel for bcswap_corner_step_1b_kernel3
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int j=1; j<ng+1; j++){
for(int i=1; i<ng+1; i++){
wbuf_lxlylz_s_gpu[__I4_WBUF_LXLYLZ_S(m,i,j,k)] = w_gpu[__I4_W(i,j,k,m)];
wbuf_lxlyrz_s_gpu[__I4_WBUF_LXLYRZ_S(m,i,j,k)] = w_gpu[__I4_W(i,j,nz-ng+k,m)];
wbuf_lxrylz_s_gpu[__I4_WBUF_LXRYLZ_S(m,i,j,k)] = w_gpu[__I4_W(i,ny-ng+j,k,m)];
wbuf_lxryrz_s_gpu[__I4_WBUF_LXRYRZ_S(m,i,j,k)] = w_gpu[__I4_W(i,ny-ng+j,nz-ng+k,m)];
wbuf_rxlylz_s_gpu[__I4_WBUF_RXLYLZ_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,j,k,m)];
wbuf_rxlyrz_s_gpu[__I4_WBUF_RXLYRZ_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,j,nz-ng+k,m)];
wbuf_rxrylz_s_gpu[__I4_WBUF_RXRYLZ_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,ny-ng+j,k,m)];
wbuf_rxryrz_s_gpu[__I4_WBUF_RXRYRZ_S(m,i,j,k)] = w_gpu[__I4_W(nx-ng+i,ny-ng+j,nz-ng+k,m)];
}
}

}
}


extern "C"{
void bcswap_corner_step_1b_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf_lxly_s_gpu,real *wbuf_lxry_s_gpu,real *wbuf_rxly_s_gpu,real *wbuf_rxry_s_gpu,real *wbuf_lylz_s_gpu,real *wbuf_lyrz_s_gpu,real *wbuf_rylz_s_gpu,real *wbuf_ryrz_s_gpu,real *wbuf_lxlylz_s_gpu,real *wbuf_lxlyrz_s_gpu,real *wbuf_lxrylz_s_gpu,real *wbuf_lxryrz_s_gpu,real *wbuf_rxlylz_s_gpu,real *wbuf_rxlyrz_s_gpu,real *wbuf_rxrylz_s_gpu,real *wbuf_rxryrz_s_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_1b_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu,wbuf_lxly_s_gpu,wbuf_lxry_s_gpu,wbuf_rxly_s_gpu,wbuf_rxry_s_gpu,wbuf_lylz_s_gpu,wbuf_lyrz_s_gpu,wbuf_rylz_s_gpu,wbuf_ryrz_s_gpu,wbuf_lxlylz_s_gpu,wbuf_lxlyrz_s_gpu,wbuf_lxrylz_s_gpu,wbuf_lxryrz_s_gpu,wbuf_rxlylz_s_gpu,wbuf_rxlyrz_s_gpu,wbuf_rxrylz_s_gpu,wbuf_rxryrz_s_gpu);
}
}



__global__ void  bcswap_step_3_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i-ng,j,k,m)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}


__global__ void  bcswap_step_3_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+i,j,k,m)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}


__global__ void  bcswap_step_3_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j-ng,k,m)] = wbuf3r_gpu[__I4_WBUF3R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}


__global__ void  bcswap_step_3_kernel4(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel4
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+j,k,m)] = wbuf4r_gpu[__I4_WBUF4R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}


__global__ void  bcswap_step_3_kernel5(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel5
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,k-ng,m)] = wbuf5r_gpu[__I4_WBUF5R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel5),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}


__global__ void  bcswap_step_3_kernel6(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
//Kernel for bcswap_step_3_kernel6
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,nz+k,m)] = wbuf6r_gpu[__I4_WBUF6R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_step_3_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_step_3_kernel6),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu);
}
}



__global__ void  bcswap_c2_step_3_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
w_gpu[__I4_W(i-ng,j,k,1)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,1)];
w_gpu[__I4_W(i-ng,j,k,2)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,2)]*dxdcsinc2_gpu[__I2_DXDCSINC2(i-ng,j)]+wbuf1r_gpu[__I4_WBUF1R(i,j,k,3)]*dxdetanc2_gpu[__I2_DXDETANC2(i-ng,j)];
w_gpu[__I4_W(i-ng,j,k,3)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,2)]*dydcsinc2_gpu[__I2_DYDCSINC2(i-ng,j)]+wbuf1r_gpu[__I4_WBUF1R(i,j,k,3)]*dydetanc2_gpu[__I2_DYDETANC2(i-ng,j)];
w_gpu[__I4_W(i-ng,j,k,4)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,4)];
w_gpu[__I4_W(i-ng,j,k,5)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,5)];

}
}


extern "C"{
void bcswap_c2_step_3_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
w_gpu[__I4_W(nx+i,j,k,1)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,1)];
w_gpu[__I4_W(nx+i,j,k,2)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,2)]*dxdcsinc2_gpu[__I2_DXDCSINC2(nx+i,j)]+wbuf2r_gpu[__I4_WBUF2R(i,j,k,3)]*dxdetanc2_gpu[__I2_DXDETANC2(nx+i,j)];
w_gpu[__I4_W(nx+i,j,k,3)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,2)]*dydcsinc2_gpu[__I2_DYDCSINC2(nx+i,j)]+wbuf2r_gpu[__I4_WBUF2R(i,j,k,3)]*dydetanc2_gpu[__I2_DYDETANC2(nx+i,j)];
w_gpu[__I4_W(nx+i,j,k,4)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,4)];
w_gpu[__I4_W(nx+i,j,k,5)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,5)];

}
}


extern "C"{
void bcswap_c2_step_3_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i-ng,j,k,m)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel4(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel4
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+i,j,k,m)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel5(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel5
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j-ng,k,m)] = wbuf3r_gpu[__I4_WBUF3R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel5),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel6(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel6
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+j,k,m)] = wbuf4r_gpu[__I4_WBUF4R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel6),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel7(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel7
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,k-ng,m)] = wbuf5r_gpu[__I4_WBUF5R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel7_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel7),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2_step_3_kernel8(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2_step_3_kernel8
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,nz+k,m)] = wbuf6r_gpu[__I4_WBUF6R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2_step_3_kernel8_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1r_gpu,real *wbuf2r_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2_step_3_kernel8),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1r_gpu,wbuf2r_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}



__global__ void  bcswap_c2xybc_step_3_kernel1(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
w_gpu[__I4_W(i-ng,j,k,1)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,1)];
w_gpu[__I4_W(i-ng,j,k,2)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,2)]*dxdcsinc2_gpu[__I2_DXDCSINC2(i-ng,j)]+wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,3)]*dxdetanc2_gpu[__I2_DXDETANC2(i-ng,j)];
w_gpu[__I4_W(i-ng,j,k,3)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,2)]*dydcsinc2_gpu[__I2_DYDCSINC2(i-ng,j)]+wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,3)]*dydetanc2_gpu[__I2_DYDETANC2(i-ng,j)];
w_gpu[__I4_W(i-ng,j,k,4)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,4)];
w_gpu[__I4_W(i-ng,j,k,5)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,5)];

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel2(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
w_gpu[__I4_W(nx+i,j,k,1)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,1)];
w_gpu[__I4_W(nx+i,j,k,2)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,2)]*dxdcsinc2_gpu[__I2_DXDCSINC2(nx+i,j)]+wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,3)]*dxdetanc2_gpu[__I2_DXDETANC2(nx+i,j)];
w_gpu[__I4_W(nx+i,j,k,3)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,2)]*dydcsinc2_gpu[__I2_DYDCSINC2(nx+i,j)]+wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,3)]*dydetanc2_gpu[__I2_DYDETANC2(nx+i,j)];
w_gpu[__I4_W(nx+i,j,k,4)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,4)];
w_gpu[__I4_W(nx+i,j,k,5)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,5)];

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel3(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i-ng,j,k,m)] = wbuf1xybcr_gpu[__I4_WBUF1XYBCR(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel4(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel4
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1-ng);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny+ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+i,j,k,m)] = wbuf2xybcr_gpu[__I4_WBUF2XYBCR(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny+ng)-(1-ng)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel5(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel5
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j-ng,k,m)] = wbuf3r_gpu[__I4_WBUF3R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel5),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel6(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel6
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+j,k,m)] = wbuf4r_gpu[__I4_WBUF4R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel6),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel7(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel7
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,k-ng,m)] = wbuf5r_gpu[__I4_WBUF5R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel7_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel7),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}


__global__ void  bcswap_c2xybc_step_3_kernel8(int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcswap_c2xybc_step_3_kernel8
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,ng,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,nz+k,m)] = wbuf6r_gpu[__I4_WBUF6R(i,j,k,m)];
}

}
}


extern "C"{
void bcswap_c2xybc_step_3_kernel8_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ndim,int ileftx,int ilefty,int ileftz,int irightx,int irighty,int irightz,real *w_gpu,real *wbuf1xybcr_gpu,real *wbuf2xybcr_gpu,real *wbuf3r_gpu,real *wbuf4r_gpu,real *wbuf5r_gpu,real *wbuf6r_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcswap_c2xybc_step_3_kernel8),grid,block,0,stream,nx,ny,nz,ng,nv,ndim,ileftx,ilefty,ileftz,irightx,irighty,irightz,w_gpu,wbuf1xybcr_gpu,wbuf2xybcr_gpu,wbuf3r_gpu,wbuf4r_gpu,wbuf5r_gpu,wbuf6r_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}



__global__ void  bcswap_corner_step_3_kernel1(int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
//Kernel for bcswap_corner_step_3_kernel1
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(j,ny,1)){
for(int k=1; k<ng+1; k++){
for(int i=1; i<ng+1; i++){
w_gpu[__I4_W(i-ng,j,k-ng,m)] = wbuf1r_c_gpu[__I4_WBUF1R_C(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu);
}
}


__global__ void  bcswap_corner_step_3_kernel2(int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
//Kernel for bcswap_corner_step_3_kernel2
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(j,ny,1)){
for(int k=1; k<ng+1; k++){
for(int i=1; i<ng+1; i++){
w_gpu[__I4_W(i+nx,j,k+nz,m)] = wbuf2r_c_gpu[__I4_WBUF2R_C(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu);
}
}


__global__ void  bcswap_corner_step_3_kernel3(int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
//Kernel for bcswap_corner_step_3_kernel3
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(j,ny,1)){
for(int k=1; k<ng+1; k++){
for(int i=1; i<ng+1; i++){
w_gpu[__I4_W(i-ng,j,k+nz,m)] = wbuf3r_c_gpu[__I4_WBUF3R_C(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu);
}
}


__global__ void  bcswap_corner_step_3_kernel4(int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
//Kernel for bcswap_corner_step_3_kernel4
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(j,ny,1)){
for(int k=1; k<ng+1; k++){
for(int i=1; i<ng+1; i++){
w_gpu[__I4_W(i+nx,j,k-ng,m)] = wbuf4r_c_gpu[__I4_WBUF4R_C(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ileftbottom,int ilefttop,int irightbottom,int irighttop,real *w_gpu,real *wbuf1r_c_gpu,real *wbuf2r_c_gpu,real *wbuf3r_c_gpu,real *wbuf4r_c_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ileftbottom,ilefttop,irightbottom,irighttop,w_gpu,wbuf1r_c_gpu,wbuf2r_c_gpu,wbuf3r_c_gpu,wbuf4r_c_gpu);
}
}



__global__ void  bcswap_corner_step_3b_kernel1(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel1
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,nz,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j-ng,k,m)] = wbuf_lxly_r_gpu[__I4_WBUF_LXLY_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel2(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel2
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,nz,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j+ny,k,m)] = wbuf_rxry_r_gpu[__I4_WBUF_RXRY_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel3(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel3
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,nz,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j+ny,k,m)] = wbuf_lxry_r_gpu[__I4_WBUF_LXRY_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel4(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel4
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,nz,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j-ng,k,m)] = wbuf_rxly_r_gpu[__I4_WBUF_RXLY_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel5(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel5
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
i = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(i,nx,1)){
for(int k=1; k<ng+1; k++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i,j-ng,k-ng,m)] = wbuf_lylz_r_gpu[__I4_WBUF_LYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel5),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel6(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel6
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
i = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(i,nx,1)){
for(int k=1; k<ng+1; k++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i,j+ny,k+nz,m)] = wbuf_ryrz_r_gpu[__I4_WBUF_RYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel6),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel7(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel7
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
i = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(i,nx,1)){
for(int k=1; k<ng+1; k++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i,j-ng,k+nz,m)] = wbuf_lyrz_r_gpu[__I4_WBUF_LYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel7_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel7),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel8(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel8
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
i = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(i,nx,1)){
for(int k=1; k<ng+1; k++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i,j+ny,k-ng,m)] = wbuf_rylz_r_gpu[__I4_WBUF_RYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel8_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel8),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel9(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel9
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j-ng,k-ng,m)] = wbuf_lxlylz_r_gpu[__I4_WBUF_LXLYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel9_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel9),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel10(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel10
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j-ng,k+nz,m)] = wbuf_lxlyrz_r_gpu[__I4_WBUF_LXLYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel10_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel10),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel11(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel11
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j+ny,k-ng,m)] = wbuf_lxrylz_r_gpu[__I4_WBUF_LXRYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel11_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel11),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel12(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel12
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i-ng,j+ny,k+nz,m)] = wbuf_lxryrz_r_gpu[__I4_WBUF_LXRYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel12_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel12),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel13(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel13
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j-ng,k-ng,m)] = wbuf_rxlylz_r_gpu[__I4_WBUF_RXLYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel13_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel13),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel14(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel14
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j-ng,k+nz,m)] = wbuf_rxlyrz_r_gpu[__I4_WBUF_RXLYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel14_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel14),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel15(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel15
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j+ny,k-ng,m)] = wbuf_rxrylz_r_gpu[__I4_WBUF_RXRYLZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel15_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel15),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}


__global__ void  bcswap_corner_step_3b_kernel16(int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
//Kernel for bcswap_corner_step_3b_kernel16
int i;int j;int k;
int m;int iercuda;

m = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(m,nv,1)&&loop_cond(k,ng,1)){
for(int i=1; i<ng+1; i++){
for(int j=1; j<ng+1; j++){
w_gpu[__I4_W(i+nx,j+ny,k+nz,m)] = wbuf_rxryrz_r_gpu[__I4_WBUF_RXRYRZ_R(m,i,j,k)];
}
}

}
}


extern "C"{
void bcswap_corner_step_3b_kernel16_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int lxly,int lxry,int rxly,int rxry,int lylz,int lyrz,int rylz,int ryrz,int lxlylz,int lxlyrz,int lxrylz,int lxryrz,int rxlylz,int rxlyrz,int rxrylz,int rxryrz,real *w_gpu,real *wbuf_lxly_r_gpu,real *wbuf_lxry_r_gpu,real *wbuf_rxly_r_gpu,real *wbuf_rxry_r_gpu,real *wbuf_lylz_r_gpu,real *wbuf_lyrz_r_gpu,real *wbuf_rylz_r_gpu,real *wbuf_ryrz_r_gpu,real *wbuf_lxlylz_r_gpu,real *wbuf_lxlyrz_r_gpu,real *wbuf_lxrylz_r_gpu,real *wbuf_lxryrz_r_gpu,real *wbuf_rxlylz_r_gpu,real *wbuf_rxlyrz_r_gpu,real *wbuf_rxrylz_r_gpu,real *wbuf_rxryrz_r_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nv)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y));

hipLaunchKernelGGL((bcswap_corner_step_3b_kernel16),grid,block,0,stream,nx,ny,nz,ng,nv,lxly,lxry,rxly,rxry,lylz,lyrz,rylz,ryrz,lxlylz,lxlyrz,lxrylz,lxryrz,rxlylz,rxlyrz,rxrylz,rxryrz,w_gpu,wbuf_lxly_r_gpu,wbuf_lxry_r_gpu,wbuf_rxly_r_gpu,wbuf_rxry_r_gpu,wbuf_lylz_r_gpu,wbuf_lyrz_r_gpu,wbuf_rylz_r_gpu,wbuf_ryrz_r_gpu,wbuf_lxlylz_r_gpu,wbuf_lxlyrz_r_gpu,wbuf_lxrylz_r_gpu,wbuf_lxryrz_r_gpu,wbuf_rxlylz_r_gpu,wbuf_rxlyrz_r_gpu,wbuf_rxrylz_r_gpu,wbuf_rxryrz_r_gpu);
}
}



__global__ void  extr_corner_ymin_kernel1(int nx,int ny,int nz,int ng,int nv,real *w_gpu){
//Kernel for extr_corner_ymin_kernel1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(1-i,1-j,k,m)] = w_gpu[__I4_W(1-i,1,k,m)];
}

}
}


extern "C"{
void extr_corner_ymin_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((extr_corner_ymin_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu);
}
}


__global__ void  extr_corner_ymin_kernel2(int nx,int ny,int nz,int ng,int nv,real *w_gpu){
//Kernel for extr_corner_ymin_kernel2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ng,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+i,1-j,k,m)] = w_gpu[__I4_W(nx+i,1,k,m)];
}

}
}


extern "C"{
void extr_corner_ymin_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((extr_corner_ymin_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,w_gpu);
}
}


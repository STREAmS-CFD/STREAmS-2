#include <hip/hip_runtime.h>
#define INCLUDE_KERNELS
#include "../hip_utils.h"  
#include "../amd_arrays.h"



__global__ void  zero_flux_kernel(int nx,int ny,int nz,int nv,real *fl_gpu){
//Kernel for zero_flux_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
fl_gpu[__I4_FL(i,j,k,m)] = 0.0;
}

}
}


extern "C"{
void zero_flux_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,real *fl_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((zero_flux_kernel),grid,block,0,stream,nx,ny,nz,nv,fl_gpu);
}
}



__global__ void  init_flux_kernel(int nx,int ny,int nz,int nv,real rhodt,real *fl_gpu,real *fln_gpu){
//Kernel for init_flux_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
fln_gpu[__I4_FLN(i,j,k,m)] = - rhodt * fl_gpu[__I4_FL(i,j,k,m)];
fl_gpu[__I4_FL(i,j,k,m)] = 0.0;
}

}
}


extern "C"{
void init_flux_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,real rhodt,real *fl_gpu,real *fln_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((init_flux_kernel),grid,block,0,stream,nx,ny,nz,nv,rhodt,fl_gpu,fln_gpu);
}
}




__global__ void  count_weno_kernel1_count_weno_x(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_x,real *redn_3d_gpu){
//Kernel for count_weno_kernel1_count_weno_x
int i;int j;int k;
int ishk;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,eul_imin);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,eul_imax,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
ishk = 0;
for(int ii=i-weno_scheme+1; ii<i+weno_scheme+1; ii++){
if (w_aux_gpu[__I4_W_AUX(ii,j,k,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_kernel1_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_x,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((eul_imax)-(eul_imin)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_kernel1_count_weno_x),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_x,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_x);


}
}



__global__ void  count_weno_kernel2_count_weno_y(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_y,real *redn_3d_gpu){
//Kernel for count_weno_kernel2_count_weno_y
int i;int j;int k;
int ishk;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,eul_jmin);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,eul_jmax,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
ishk = 0;
for(int jj=j-weno_scheme+1; jj<j+weno_scheme+1; jj++){
if (w_aux_gpu[__I4_W_AUX(i,jj,k,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_kernel2_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_y,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((eul_jmax)-(eul_jmin)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_kernel2_count_weno_y),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_y,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_y);


}
}



__global__ void  count_weno_kernel3_count_weno_z(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_z,real *redn_3d_gpu){
//Kernel for count_weno_kernel3_count_weno_z
int i;int j;int k;
int ishk;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,eul_kmin);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,eul_kmax,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
ishk = 0;
for(int kk=k-weno_scheme+1; kk<k+weno_scheme+1; kk++){
if (w_aux_gpu[__I4_W_AUX(i,j,kk,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_kernel3_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,real *w_aux_gpu,real *count_weno_z,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((eul_kmax)-(eul_kmin)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_kernel3_count_weno_z),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,weno_scheme,sensor_threshold,ep_ord_change_gpu,w_aux_gpu,count_weno_z,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_z);


}
}




__global__ void  count_weno_c2_kernel1_count_weno_x(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_x,real *redn_3d_gpu){
//Kernel for count_weno_c2_kernel1_count_weno_x
int i;int j;int k;
int ishk;int weno_scheme_i;int weno_scheme_j;
int weno_scheme_k;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,eul_imin);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,eul_imax,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
weno_scheme_i = weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)];
if(weno_scheme_i < 1) {
weno_scheme_i = 1;
}
if(j == 1) {
weno_scheme_i = lmax_tag_gpu[__I1_LMAX_TAG(i)]+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)];
if(weno_scheme_i < 1) {
weno_scheme_i = 1;
}
}
ishk = 0;
for(int ii=i-weno_scheme_i+1; ii<i+weno_scheme_i+1; ii++){
if (w_aux_gpu[__I4_W_AUX(ii,j,k,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_c2_kernel1_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_x,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((eul_imax)-(eul_imin)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_c2_kernel1_count_weno_x),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_x,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_x);


}
}



__global__ void  count_weno_c2_kernel2_count_weno_y(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_y,real *redn_3d_gpu){
//Kernel for count_weno_c2_kernel2_count_weno_y
int i;int j;int k;
int ishk;int weno_scheme_i;int weno_scheme_j;
int weno_scheme_k;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,eul_jmin);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,eul_jmax,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
weno_scheme_j = weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)];
if(weno_scheme_j < 1) {
weno_scheme_j = 1;
}
if (j <= l_base && wall_tag_gpu[__I1_WALL_TAG(i)] > 0) {
weno_scheme_j = weno_scheme;
}
ishk = 0;
for(int jj=j-weno_scheme_j+1; jj<j+weno_scheme_j+1; jj++){
if (w_aux_gpu[__I4_W_AUX(i,jj,k,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_c2_kernel2_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_y,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((eul_jmax)-(eul_jmin)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_c2_kernel2_count_weno_y),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_y,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_y);


}
}



__global__ void  count_weno_c2_kernel3_count_weno_z(int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_z,real *redn_3d_gpu){
//Kernel for count_weno_c2_kernel3_count_weno_z
int i;int j;int k;
int ishk;int weno_scheme_i;int weno_scheme_j;
int weno_scheme_k;int ii;int jj;
int kk;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,eul_kmin);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,eul_kmax,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
ishk = 0;
for(int kk=k-weno_scheme+1; kk<k+weno_scheme+1; kk++){
if (w_aux_gpu[__I4_W_AUX(i,j,kk,8)] > sensor_threshold) {
ishk = 1;
}
}
if (ishk > 0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1;
}

}
}



extern "C"{
void count_weno_c2_kernel3_wrapper(hipStream_t stream,int nv,int nv_aux,int nx,int ny,int nz,int ng,int eul_imin,int eul_imax,int eul_jmin,int eul_jmax,int eul_kmin,int eul_kmax,int l_base,int weno_scheme,real sensor_threshold,int *ep_ord_change_gpu,int *lmax_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *count_weno_z,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((eul_kmax)-(eul_kmin)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((count_weno_c2_kernel3_count_weno_z),grid,block,0,stream,nv,nv_aux,nx,ny,nz,ng,eul_imin,eul_imax,eul_jmin,eul_jmax,eul_kmin,eul_kmax,l_base,weno_scheme,sensor_threshold,ep_ord_change_gpu,lmax_tag_gpu,wall_tag_gpu,w_aux_gpu,count_weno_z,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, count_weno_z);


}
}



__global__ void  euler_x_update_kernel(int nx,int ny,int nz,int ng,int nv,int eul_imin,int eul_imax,real *fhat_gpu,real *fl_gpu,real *dcsidx_gpu){
//Kernel for euler_x_update_kernel
int i;int j;int k;
int m;int iv;int iercuda;

i = __GIDX(x,eul_imin);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,eul_imax,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int iv=1; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + (fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i-1,j,k,iv)])*dcsidx_gpu[__I1_DCSIDX(i)];
}

}
}


extern "C"{
void euler_x_update_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int eul_imin,int eul_imax,real *fhat_gpu,real *fl_gpu,real *dcsidx_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((eul_imax)-(eul_imin)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_x_update_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,dcsidx_gpu);
}
}



__global__ void  euler_x_update_c2_kernel(int nx,int ny,int nz,int ng,int nv,int eul_imin,int eul_imax,real *fhat_gpu,real *fl_gpu,real *jac_gpu){
//Kernel for euler_x_update_c2_kernel
int i;int j;int k;
int m;int iv;int iercuda;

i = __GIDX(x,eul_imin);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,eul_imax,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int iv=1; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + (fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i-1,j,k,iv)])*jac_gpu[__I2_JAC(i,j)];
}

}
}


extern "C"{
void euler_x_update_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int eul_imin,int eul_imax,real *fhat_gpu,real *fl_gpu,real *jac_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((eul_imax)-(eul_imin)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_x_update_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,eul_imin,eul_imax,fhat_gpu,fl_gpu,jac_gpu);
}
}



__global__ void  euler_y_update_kernel(int nx,int ny,int nz,int ng,int nv,int eul_jmin,int eul_jmax,real *fhat_gpu,real *fl_gpu,real *detady_gpu){
//Kernel for euler_y_update_kernel
int i;int j;int k;
int m;int iv;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,eul_jmin);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,eul_jmax,1)&&loop_cond(k,nz,1)){
for(int iv=1; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + (fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i,j-1,k,iv)])*detady_gpu[__I1_DETADY(j)];
}

}
}


extern "C"{
void euler_y_update_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int eul_jmin,int eul_jmax,real *fhat_gpu,real *fl_gpu,real *detady_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((eul_jmax)-(eul_jmin)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_y_update_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,fhat_gpu,fl_gpu,detady_gpu);
}
}



__global__ void  euler_y_update_c2_kernel(int nx,int ny,int nz,int ng,int nv,int eul_jmin,int eul_jmax,int *wall_tag_gpu,real *fhat_gpu,real *fl_gpu,real *jac_gpu){
//Kernel for euler_y_update_c2_kernel
int i;int j;int k;
int m;int iv;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,eul_jmax,1)&&loop_cond(k,nz,1)){
if(j > 1 || eul_jmin == 1 || wall_tag_gpu[__I1_WALL_TAG(i)] > 0) {
for(int iv=1; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + (fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i,j-1,k,iv)])*jac_gpu[__I2_JAC(i,j)];
}
}

}
}


extern "C"{
void euler_y_update_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int eul_jmin,int eul_jmax,int *wall_tag_gpu,real *fhat_gpu,real *fl_gpu,real *jac_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((eul_jmax)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_y_update_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,eul_jmin,eul_jmax,wall_tag_gpu,fhat_gpu,fl_gpu,jac_gpu);
}
}



__global__ void  euler_z_update_kernel(int nx,int ny,int nz,int ng,int nv,int eul_kmin,int eul_kmax,real *fhat_gpu,real *fl_gpu,real *dzitdz_gpu){
//Kernel for euler_z_update_kernel
int i;int j;int k;
int m;int iv;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,eul_kmin);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,eul_kmax,1)){
for(int iv=1; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + (fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i,j,k-1,iv)])*dzitdz_gpu[__I1_DZITDZ(k)];
}

}
}


extern "C"{
void euler_z_update_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int eul_kmin,int eul_kmax,real *fhat_gpu,real *fl_gpu,real *dzitdz_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((eul_kmax)-(eul_kmin)+1,block.z));

hipLaunchKernelGGL((euler_z_update_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,eul_kmin,eul_kmax,fhat_gpu,fl_gpu,dzitdz_gpu);
}
}



__global__ void  force_rhs_2_kernel(int nx,int ny,int nz,int ng,real bulk_1,real bulk_2,int *fluid_mask_gpu,real *fln_gpu,real *w_aux_gpu){
//Kernel for force_rhs_2_kernel
real uu;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
fln_gpu[__I4_FLN(i,j,k,1)] = fln_gpu[__I4_FLN(i,j,k,1)] - bulk_1;
fln_gpu[__I4_FLN(i,j,k,2)] = fln_gpu[__I4_FLN(i,j,k,2)] - bulk_2;
fln_gpu[__I4_FLN(i,j,k,5)] = fln_gpu[__I4_FLN(i,j,k,5)] - uu*bulk_2;
}

}
}


extern "C"{
void force_rhs_2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real bulk_1,real bulk_2,int *fluid_mask_gpu,real *fln_gpu,real *w_aux_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((force_rhs_2_kernel),grid,block,0,stream,nx,ny,nz,ng,bulk_1,bulk_2,fluid_mask_gpu,fln_gpu,w_aux_gpu);
}
}



__global__ void  force_rhs_2_c2_kernel(int nx,int ny,int nz,int ng,real r_curv,real bulk_1,real bulk_5,int *fluid_mask_gpu,real *fln_gpu,real *yn_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *w_aux_gpu){
//Kernel for force_rhs_2_c2_kernel
real uu;real vv;real r;
real dpdcsi;real dpdx;real dpdy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
r = r_curv+0.50*(yn_gpu[__I1_YN(j)]+yn_gpu[__I1_YN(j+1)]);
dpdcsi = bulk_5/r;
dpdx = dpdcsi*dxdcsinc2_gpu[__I2_DXDCSINC2(i,j)];
dpdy = dpdcsi*dydcsinc2_gpu[__I2_DYDCSINC2(i,j)];
fln_gpu[__I4_FLN(i,j,k,1)] = fln_gpu[__I4_FLN(i,j,k,1)] - bulk_1;
fln_gpu[__I4_FLN(i,j,k,2)] = fln_gpu[__I4_FLN(i,j,k,2)] - dpdx;
fln_gpu[__I4_FLN(i,j,k,3)] = fln_gpu[__I4_FLN(i,j,k,3)] - dpdy;
fln_gpu[__I4_FLN(i,j,k,5)] = fln_gpu[__I4_FLN(i,j,k,5)] - uu*dpdx - vv*dpdy;
}

}
}


extern "C"{
void force_rhs_2_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real r_curv,real bulk_1,real bulk_5,int *fluid_mask_gpu,real *fln_gpu,real *yn_gpu,real *dxdcsinc2_gpu,real *dydcsinc2_gpu,real *w_aux_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((force_rhs_2_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,r_curv,bulk_1,bulk_5,fluid_mask_gpu,fln_gpu,yn_gpu,dxdcsinc2_gpu,dydcsinc2_gpu,w_aux_gpu);
}
}




__global__ void  force_rhs_1_kernel_bulk_1(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_kernel_bulk_1
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +fln_gpu[__I4_FLN(i,j,k,1)]*dy;
}

}
}


__global__ void  force_rhs_1_kernel_bulk_2(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_kernel_bulk_2
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +fln_gpu[__I4_FLN(i,j,k,2)]*dy;
}

}
}


__global__ void  force_rhs_1_kernel_bulk_3(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_kernel_bulk_3
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +w_gpu[__I4_W(i,j,k,1)]*dy;
}

}
}


__global__ void  force_rhs_1_kernel_bulk_4(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_kernel_bulk_4
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +w_gpu[__I4_W(i,j,k,2)]*dy;
}

}
}


__global__ void  force_rhs_1_kernel_bulk_5(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_kernel_bulk_5
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +w_gpu[__I4_W(i,j,k,2)]*dy*w_aux_gpu[__I4_W_AUX(i,j,k,6)];
}

}
}



extern "C"{
void force_rhs_1_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((force_rhs_1_kernel_bulk_1),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_1);

hipLaunchKernelGGL((force_rhs_1_kernel_bulk_2),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_2);

hipLaunchKernelGGL((force_rhs_1_kernel_bulk_3),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_3);

hipLaunchKernelGGL((force_rhs_1_kernel_bulk_4),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_4);

hipLaunchKernelGGL((force_rhs_1_kernel_bulk_5),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_5);


}
}




__global__ void  force_rhs_1_c2_kernel_bulk_1(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_c2_kernel_bulk_1
real dv;real rhou;real f_rhou;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rhou = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+ w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
f_rhou = fln_gpu[__I4_FLN(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+fln_gpu[__I4_FLN(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +fln_gpu[__I4_FLN(i,j,k,1)]*dv;
}

}
}


__global__ void  force_rhs_1_c2_kernel_bulk_2(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_c2_kernel_bulk_2
real dv;real rhou;real f_rhou;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rhou = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+ w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
f_rhou = fln_gpu[__I4_FLN(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+fln_gpu[__I4_FLN(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +w_gpu[__I4_W(i,j,k,1)]*dv;
}

}
}


__global__ void  force_rhs_1_c2_kernel_bulk_3(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_c2_kernel_bulk_3
real dv;real rhou;real f_rhou;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rhou = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+ w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
f_rhou = fln_gpu[__I4_FLN(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+fln_gpu[__I4_FLN(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +rhou*dv;
}

}
}


__global__ void  force_rhs_1_c2_kernel_bulk_4(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_c2_kernel_bulk_4
real dv;real rhou;real f_rhou;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rhou = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+ w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
f_rhou = fln_gpu[__I4_FLN(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+fln_gpu[__I4_FLN(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +rhou*w_aux_gpu[__I4_W_AUX(i,j,k,6)]*dv;
}

}
}


__global__ void  force_rhs_1_c2_kernel_bulk_5(int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_rhs_1_c2_kernel_bulk_5
real dv;real rhou;real f_rhou;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rhou = w_gpu[__I4_W(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+ w_gpu[__I4_W(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
f_rhou = fln_gpu[__I4_FLN(i,j,k,2)]*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+fln_gpu[__I4_FLN(i,j,k,3)]*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +f_rhou*dv;
}

}
}



extern "C"{
void force_rhs_1_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *w_gpu,real *w_aux_gpu,real *bulk_1,real *bulk_2,real *bulk_3,real *bulk_4,real *bulk_5,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((force_rhs_1_c2_kernel_bulk_1),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_1);

hipLaunchKernelGGL((force_rhs_1_c2_kernel_bulk_2),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_2);

hipLaunchKernelGGL((force_rhs_1_c2_kernel_bulk_3),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_3);

hipLaunchKernelGGL((force_rhs_1_c2_kernel_bulk_4),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_4);

hipLaunchKernelGGL((force_rhs_1_c2_kernel_bulk_5),grid,block,0,stream,nx,ny,nz,ng,fluid_mask_gpu,yn_gpu,fln_gpu,dcsidxnc2_gpu,dcsidync2_gpu,jac_gpu,w_gpu,w_aux_gpu,bulk_1,bulk_2,bulk_3,bulk_4,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_5);


}
}

__device__ real get_temperature_from_e_dev_force_var_1_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_force_var_1_kernel_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}





__global__ void  force_var_1_kernel_bulk_5(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tol_iter_nr,real bulkt,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_var_1_kernel_bulk_5
real rho;real rhou;real rhov;
real rhow;real rhoe;real ri;
real uu;real vv;real ww;
real qq;real ee;real tt;
real dy;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dy = yn_gpu[__I1_YN(j+1)]-yn_gpu[__I1_YN(j)];
rho = w_gpu[__I4_W(i,j,k,1)];
rhou = w_gpu[__I4_W(i,j,k,2)];
rhov = w_gpu[__I4_W(i,j,k,3)];
rhow = w_gpu[__I4_W(i,j,k,4)];
rhoe = w_gpu[__I4_W(i,j,k,5)];
ri = 1.0/rho;
uu = rhou*ri;
vv = rhov*ri;
ww = rhow*ri;
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_force_var_1_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
w_aux_gpu[__I4_W_AUX(i,j,k,6)] = tt;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +rhou*tt*dy;
}

}
}



extern "C"{
void force_var_1_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tol_iter_nr,real bulkt,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *bulk_5,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((force_var_1_kernel_bulk_5),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,t0,tol_iter_nr,bulkt,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_5);


}
}

__device__ real get_e_from_temperature_dev_force_var_2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_force_var_2_kernel_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  force_var_2_kernel(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real tbdiff,real t0,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
//Kernel for force_var_2_kernel
real rho;real rhou;real rhov;
real rhow;real tt;real ttnew;
real ee;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
rhou = w_gpu[__I4_W(i,j,k,2)];
rhov = w_gpu[__I4_W(i,j,k,3)];
rhow = w_gpu[__I4_W(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
ttnew = tt+tbdiff;
ee = get_e_from_temperature_dev_force_var_2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,ttnew,cv_coeff_gpu);
w_gpu[__I4_W(i,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;
}

}
}


extern "C"{
void force_var_2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real tbdiff,real t0,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((force_var_2_kernel),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,tbdiff,t0,fluid_mask_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu);
}
}

__device__ real get_temperature_from_e_dev_force_var_1_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_force_var_1_c2_kernel_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}





__global__ void  force_var_1_c2_kernel_bulk_5(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tol_iter_nr,real bulkt,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *jac_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *bulk_5,real *redn_3d_gpu){
//Kernel for force_var_1_c2_kernel_bulk_5
real rho;real rhou;real rhov;
real rhow;real rhoe;real ri;
real uu;real vv;real ww;
real qq;real ee;real tt;
real dy;real dv;real rhout;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rho = w_gpu[__I4_W(i,j,k,1)];
rhou = w_gpu[__I4_W(i,j,k,2)];
rhov = w_gpu[__I4_W(i,j,k,3)];
rhow = w_gpu[__I4_W(i,j,k,4)];
rhoe = w_gpu[__I4_W(i,j,k,5)];
rhout = rhou*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+rhov*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
ri = 1.0/rho;
uu = rhou*ri;
vv = rhov*ri;
ww = rhow*ri;
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_force_var_1_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
w_aux_gpu[__I4_W_AUX(i,j,k,6)] = tt;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +rhou*tt*dv;
}

}
}



extern "C"{
void force_var_1_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tol_iter_nr,real bulkt,int *fluid_mask_gpu,real *yn_gpu,real *fln_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *jac_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *bulk_5,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((force_var_1_c2_kernel_bulk_5),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,t0,tol_iter_nr,bulkt,fluid_mask_gpu,yn_gpu,fln_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu,dcsidxnc2_gpu,dcsidync2_gpu,bulk_5,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, bulk_5);


}
}

__device__ real get_e_from_temperature_dev_force_var_2_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_force_var_2_c2_kernel_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  force_var_2_c2_kernel(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real tbdiff,real t0,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *jac_gpu){
//Kernel for force_var_2_c2_kernel
real rho;real rhou;real rhov;
real rhow;real tt;real ttnew;
real ee;real dv;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
dv = 1.0/jac_gpu[__I2_JAC(i,j)];
rho = w_gpu[__I4_W(i,j,k,1)];
rhou = w_gpu[__I4_W(i,j,k,2)];
rhov = w_gpu[__I4_W(i,j,k,3)];
rhow = w_gpu[__I4_W(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
ttnew = tt+tbdiff;
ee = get_e_from_temperature_dev_force_var_2_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,ttnew,cv_coeff_gpu);
w_gpu[__I4_W(i,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;
}

}
}


extern "C"{
void force_var_2_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real tbdiff,real t0,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *jac_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((force_var_2_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,tbdiff,t0,fluid_mask_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,jac_gpu);
}
}



__global__ void  update_flux_kernel(int nx,int ny,int nz,int nv,real gamdt,real *fl_gpu,real *fln_gpu){
//Kernel for update_flux_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
fln_gpu[__I4_FLN(i,j,k,m)] = fln_gpu[__I4_FLN(i,j,k,m)]-gamdt*fl_gpu[__I4_FL(i,j,k,m)];
}

}
}


extern "C"{
void update_flux_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,real gamdt,real *fl_gpu,real *fln_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((update_flux_kernel),grid,block,0,stream,nx,ny,nz,nv,gamdt,fl_gpu,fln_gpu);
}
}



__global__ void  update_field_kernel(int nx,int ny,int nz,int nv,int ng,int *fluid_mask_gpu,real *w_gpu,real *fln_gpu){
//Kernel for update_field_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,k,m)] = w_gpu[__I4_W(i,j,k,m)]+fln_gpu[__I4_FLN(i,j,k,m)];
}
}

}
}


extern "C"{
void update_field_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int *fluid_mask_gpu,real *w_gpu,real *fln_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((update_field_kernel),grid,block,0,stream,nx,ny,nz,nv,ng,fluid_mask_gpu,w_gpu,fln_gpu);
}
}



__global__ void  visflx_div_ord2_kernel(int nx,int ny,int nz,int ng,real *w_aux_gpu,real *fl_gpu,real *x_gpu,real *y_gpu,real *z_gpu){
//Kernel for visflx_div_ord2_kernel
real dxl;real dyl;real dzl;
real uu;real vv;real ww;
real mu;real sigq;real sigx;
real sigy;real sigz;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
dxl = mu/(x_gpu[__I1_X(i+1)]-x_gpu[__I1_X(i-1)]);
dyl = mu/(y_gpu[__I1_Y(j+1)]-y_gpu[__I1_Y(j-1)]);
dzl = mu/(z_gpu[__I1_Z(k+1)]-z_gpu[__I1_Z(k-1)]);
sigx = dxl*(w_aux_gpu[__I4_W_AUX(i+1,j,k,10)]-w_aux_gpu[__I4_W_AUX(i-1,j,k,10)]);
sigy = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,10)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,10)]);
sigz = dzl*(w_aux_gpu[__I4_W_AUX(i,j,k+1,10)]-w_aux_gpu[__I4_W_AUX(i,j,k-1,10)]);
sigq = sigx*uu+sigy*vv+sigz*ww;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_div_ord2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_aux_gpu,real *fl_gpu,real *x_gpu,real *y_gpu,real *z_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_div_ord2_kernel),grid,block,0,stream,nx,ny,nz,ng,w_aux_gpu,fl_gpu,x_gpu,y_gpu,z_gpu);
}
}



__global__ void  visflx_div_kernel(int nx,int ny,int nz,int ng,int visc_order,int lmax,real *w_aux_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu){
//Kernel for visflx_div_kernel
real ccl;real uu;real vv;
real ww;real tt;real mu;
real sigq;real sigx;real sigy;
real sigz;real divx3l;real divy3l;
real divz3l;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
divx3l = 0.0;
divy3l = 0.0;
divz3l = 0.0;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
divx3l = divx3l+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,10)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,10)]);
divy3l = divy3l+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,10)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,10)]);
divz3l = divz3l+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,10)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,10)]);
}
divx3l = divx3l*dcsidx_gpu[__I1_DCSIDX(i)];
divy3l = divy3l*detady_gpu[__I1_DETADY(j)];
divz3l = divz3l*dzitdz_gpu[__I1_DZITDZ(k)];
sigx = mu*divx3l;
sigy = mu*divy3l;
sigz = mu*divz3l;
sigq = sigx*uu+sigy*vv+sigz*ww;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_div_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int lmax,real *w_aux_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_div_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,lmax,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu);
}
}



__global__ void  visflx_div_c2_kernel(int nx,int ny,int nz,int ng,int visc_order,int lmax,int *vis_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *dcsidyc2_gpu,real *detadxc2_gpu,real *dzitdz_gpu){
//Kernel for visflx_div_c2_kernel
real cli;real clj;real clk;
real uu;real vv;real ww;
real tt;real mu;real sigq;
real sigx;real sigy;real sigz;
real divcsi;real diveta;real divzit;
real divx3l;real divy3l;real divz3l;
int i;int j;int k;
int l;int iercuda;int lmaxi;
int lmaxj;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
lmaxi = lmax;
lmaxj = lmax;
if (j==1) lmaxi = vis_tag_gpu[__I1_VIS_TAG(i)];
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) lmaxj = min(j,lmax);
divcsi = 0.0;
diveta = 0.0;
divzit = 0.0;
for(int l=1; l<lmax+1; l++){
cli = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxi)];
clj = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxj)];
clk = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax )];
divcsi = divcsi+cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,10)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,10)]);
diveta = diveta+clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,10)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,10)]);
divzit = divzit+clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,10)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,10)]);
}
divx3l = divcsi*dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + diveta*detadxc2_gpu[__I2_DETADXC2(i,j)];
divy3l = divcsi*dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + diveta*detadyc2_gpu[__I2_DETADYC2(i,j)];
divz3l = divzit*dzitdz_gpu[__I1_DZITDZ(k)];
sigx = mu*divx3l;
sigy = mu*divy3l;
sigz = mu*divz3l;
sigq = sigx*uu+sigy*vv+sigz*ww;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_div_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int lmax,int *vis_tag_gpu,int *wall_tag_gpu,real *w_aux_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *dcsidyc2_gpu,real *detadxc2_gpu,real *dzitdz_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_div_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,lmax,vis_tag_gpu,wall_tag_gpu,w_aux_gpu,fl_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,detadyc2_gpu,dcsidyc2_gpu,detadxc2_gpu,dzitdz_gpu);
}
}



__global__ void  visflx_kernel(int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int lmax,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dcsidx2_gpu,real *detady2_gpu,real *dzitdz2_gpu,real *cp_coeff_gpu){
//Kernel for visflx_kernel
real ccl;real clapl;real sig11;
real sig12;real sig13;real sig22;
real sig23;real sig33;real uu;
real vv;real ww;real tt;
real mu;real ux;real uy;
real uz;real vx;real vy;
real vz;real wx;real wy;
real wz;real tx;real ty;
real tz;real mux;real muy;
real muz;real ulap;real ulapx;
real ulapy;real ulapz;real vlap;
real vlapx;real vlapy;real vlapz;
real wlap;real wlapx;real wlapy;
real wlapz;real tlap;real tlapx;
real tlapy;real tlapz;real sigq;
real sigx;real sigy;real sigz;
real sigqt;real sigah;real div;
real div3l;real omegax;real omegay;
real omegaz;real omod2;real cploc;
int i;int j;int k;
int l;int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
}
ux = 0.0;
vx = 0.0;
wx = 0.0;
tx = 0.0;
mux = 0.0;
uy = 0.0;
vy = 0.0;
wy = 0.0;
ty = 0.0;
muy = 0.0;
uz = 0.0;
vz = 0.0;
wz = 0.0;
tz = 0.0;
muz = 0.0;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
ux = ux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vx = vx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wx = wx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
tx = tx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
mux = mux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,7)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,7)]);
uy = uy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
vy = vy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
wy = wy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
ty = ty+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
muy = muy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,7)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,7)]);
uz = uz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vz = vz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wz = wz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tz = tz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
muz = muz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,7)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,7)]);
}
ux = ux*dcsidx_gpu[__I1_DCSIDX(i)];
vx = vx*dcsidx_gpu[__I1_DCSIDX(i)];
wx = wx*dcsidx_gpu[__I1_DCSIDX(i)];
tx = tx*dcsidx_gpu[__I1_DCSIDX(i)];
mux = mux*dcsidx_gpu[__I1_DCSIDX(i)];
uy = uy*detady_gpu[__I1_DETADY(j)];
vy = vy*detady_gpu[__I1_DETADY(j)];
wy = wy*detady_gpu[__I1_DETADY(j)];
ty = ty*detady_gpu[__I1_DETADY(j)];
muy = muy*detady_gpu[__I1_DETADY(j)];
uz = uz*dzitdz_gpu[__I1_DZITDZ(k)];
vz = vz*dzitdz_gpu[__I1_DZITDZ(k)];
wz = wz*dzitdz_gpu[__I1_DZITDZ(k)];
tz = tz*dzitdz_gpu[__I1_DZITDZ(k)];
muz = muz*dzitdz_gpu[__I1_DZITDZ(k)];
if (j==1) {
wallprop_gpu[__I3_WALLPROP(i,k,2)] = mu*uy;
wallprop_gpu[__I3_WALLPROP(i,k,3)] = mu*wy;
wallprop_gpu[__I3_WALLPROP(i,k,4)] = mu*ty*cploc/prandtl;
}
ulapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*uu;
ulapy = ulapx;
ulapz = ulapx;
vlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*vv;
vlapy = vlapx;
vlapz = vlapx;
wlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*ww;
wlapy = wlapx;
wlapz = wlapx;
tlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*tt;
tlapy = tlapx;
tlapz = tlapx;
for(int l=1; l<lmax+1; l++){
clapl = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax)];
ulapx = ulapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
ulapy = ulapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
ulapz = ulapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vlapx = vlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
vlapy = vlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
vlapz = vlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wlapx = wlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
wlapy = wlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
wlapz = wlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tlapx = tlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
tlapy = tlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
tlapz = tlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
}
ulapx = ulapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+ux*dcsidx2_gpu[__I1_DCSIDX2(i)];
vlapx = vlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+vx*dcsidx2_gpu[__I1_DCSIDX2(i)];
wlapx = wlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+wx*dcsidx2_gpu[__I1_DCSIDX2(i)];
tlapx = tlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+tx*dcsidx2_gpu[__I1_DCSIDX2(i)];
ulapy = ulapy*detadys_gpu[__I1_DETADYS(j)]+uy*detady2_gpu[__I1_DETADY2(j)];
vlapy = vlapy*detadys_gpu[__I1_DETADYS(j)]+vy*detady2_gpu[__I1_DETADY2(j)];
wlapy = wlapy*detadys_gpu[__I1_DETADYS(j)]+wy*detady2_gpu[__I1_DETADY2(j)];
tlapy = tlapy*detadys_gpu[__I1_DETADYS(j)]+ty*detady2_gpu[__I1_DETADY2(j)];
ulapz = ulapz*dzitdzs_gpu[__I1_DZITDZS(k)]+uz*dzitdz2_gpu[__I1_DZITDZ2(k)];
vlapz = vlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+vz*dzitdz2_gpu[__I1_DZITDZ2(k)];
wlapz = wlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+wz*dzitdz2_gpu[__I1_DZITDZ2(k)];
tlapz = tlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+tz*dzitdz2_gpu[__I1_DZITDZ2(k)];
ulap = ulapx+ulapy+ulapz;
vlap = vlapx+vlapy+vlapz;
wlap = wlapx+wlapy+wlapz;
tlap = tlapx+tlapy+tlapz;
div = ux+vy+wz;
div3l = div/3.0;
w_aux_gpu[__I4_W_AUX(i,j,k,10)] = div3l;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
w_aux_gpu[__I4_W_AUX(i,j,k,9)] = sqrt(omod2);
w_aux_gpu[__I4_W_AUX(i,j,k,8)]=(((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0)))*((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0))));
sig11 = 2.0*(ux-div3l);
sig12 = uy+vx;
sig13 = uz+wx;
sig22 = 2.0*(vy-div3l);
sig23 = vz+wy;
sig33 = 2.0*(wz-div3l);
sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap;
sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap;
sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap;
sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl;
sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu;
sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int lmax,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dcsidx2_gpu,real *detady2_gpu,real *dzitdz2_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,lmax,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,cp_coeff_gpu);
}
}



__global__ void  visflx_c2_kernel(int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int lmax,real prandtl,real u0,real l0,real t0,real teshk,int *vis_tag_gpu,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dzitdz2_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *g1_gpu,real *g2_gpu,real *g12_gpu,real *jac_gpu,real *cp_coeff_gpu){
//Kernel for visflx_c2_kernel
real ccl;real clapl;real cli;
real clj;real clk;real clapi;
real clapj;real clapk;real sig11;
real sig12;real sig13;real sig22;
real sig23;real sig33;real uu;
real vv;real ww;real tt;
real mu;real ux;real uy;
real uz;real vx;real vy;
real vz;real wx;real wy;
real wz;real tx;real ty;
real tz;real mux;real muy;
real muz;real ucsi;real ueta;
real uzit;real vcsi;real veta;
real vzit;real wcsi;real weta;
real wzit;real tcsi;real teta;
real tzit;real mucsi;real mueta;
real muzit;real ulapcsi;real ulapeta;
real ulapzit;real vlapcsi;real vlapeta;
real vlapzit;real wlapcsi;real wlapeta;
real wlapzit;real ulapcsieta;real vlapcsieta;
real wlapcsieta;real tlapcsieta;real dg1;
real dg2;real dg12csi;real dg12eta;
real tlapcsi;real tlapeta;real tlapzit;
real ulap;real ulapx;real ulapy;
real ulapz;real vlap;real vlapx;
real vlapy;real vlapz;real wlap;
real wlapx;real wlapy;real wlapz;
real tlap;real tlapx;real tlapy;
real tlapz;real sigq;real sigx;
real sigy;real sigz;real sigqt;
real sigah;real div;real div3l;
real omegax;real omegay;real omegaz;
real omod2;real cploc;
int i;int j;int k;
int l;int ll;int iercuda;
int lmaxi;int lmaxj;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
lmaxi = lmax;
lmaxj = lmax;
if (j == 1) lmaxi = vis_tag_gpu[__I1_VIS_TAG(i)];
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) lmaxj = min(j,lmax);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
}
ucsi = 0.0;
vcsi = 0.0;
wcsi = 0.0;
tcsi = 0.0;
mucsi = 0.0;
ueta = 0.0;
veta = 0.0;
weta = 0.0;
teta = 0.0;
mueta = 0.0;
uzit = 0.0;
vzit = 0.0;
wzit = 0.0;
tzit = 0.0;
muzit = 0.0;
dg1 = 0.0;
dg2 = 0.0;
dg12csi = 0.0;
dg12eta = 0.0;
for(int l=1; l<lmax+1; l++){
cli = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxi)];
clj = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxj)];
clk = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax )];
ucsi = ucsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vcsi = vcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wcsi = wcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
tcsi = tcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
mucsi = mucsi+cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,7)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,7)]);
ueta = ueta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
veta = veta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
weta = weta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
teta = teta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
mueta = mueta+clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,7)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,7)]);
uzit = uzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vzit = vzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wzit = wzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tzit = tzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
muzit = muzit+clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,7)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,7)]);
dg1 = dg1 +cli*(g1_gpu[__I2_G1(i+l,j)]-g1_gpu[__I2_G1(i-l,j)]);
dg12csi=dg12csi+cli*(g12_gpu[__I2_G12(i+l,j)]-g12_gpu[__I2_G12(i-l,j)]);
dg2 = dg2 +clj*(g2_gpu[__I2_G2(i,j+l)]-g2_gpu[__I2_G2(i,j-l)]);
dg12eta=dg12eta+clj*(g12_gpu[__I2_G12(i,j+l)]-g12_gpu[__I2_G12(i,j-l)]);
}
ux = ucsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + ueta *detadxc2_gpu[__I2_DETADXC2(i,j)];
vx = vcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + veta *detadxc2_gpu[__I2_DETADXC2(i,j)];
wx = wcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + weta *detadxc2_gpu[__I2_DETADXC2(i,j)];
tx = tcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + teta *detadxc2_gpu[__I2_DETADXC2(i,j)];
mux = mucsi*dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + mueta*detadxc2_gpu[__I2_DETADXC2(i,j)];
uy = ucsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + ueta *detadyc2_gpu[__I2_DETADYC2(i,j)];
vy = vcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + veta *detadyc2_gpu[__I2_DETADYC2(i,j)];
wy = wcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + weta *detadyc2_gpu[__I2_DETADYC2(i,j)];
ty = tcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + teta *detadyc2_gpu[__I2_DETADYC2(i,j)];
muy = mucsi*dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + mueta*detadyc2_gpu[__I2_DETADYC2(i,j)];
uz = uzit *dzitdz_gpu[__I1_DZITDZ(k)];
vz = vzit *dzitdz_gpu[__I1_DZITDZ(k)];
wz = wzit *dzitdz_gpu[__I1_DZITDZ(k)];
tz = tzit *dzitdz_gpu[__I1_DZITDZ(k)];
muz = muzit*dzitdz_gpu[__I1_DZITDZ(k)];
if (j==1) {
wallprop_gpu[__I3_WALLPROP(i,k,2)] = mu*ueta;
wallprop_gpu[__I3_WALLPROP(i,k,3)] = mu*weta;
wallprop_gpu[__I3_WALLPROP(i,k,4)] = mu*teta*cploc/prandtl;
}
ulapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*uu;
ulapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*uu;
ulapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*uu;
vlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*vv;
vlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*vv;
vlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*vv;
wlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*ww;
wlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*ww;
wlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*ww;
tlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*tt;
tlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*tt;
tlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*tt;
for(int l=1; l<lmax+1; l++){
clapi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmaxi)];
clapj = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmaxj)];
clapk = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax )];
ulapcsi = ulapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
ulapeta = ulapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
ulapzit = ulapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vlapcsi = vlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
vlapeta = vlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
vlapzit = vlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wlapcsi = wlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
wlapeta = wlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
wlapzit = wlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tlapcsi = tlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
tlapeta = tlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
tlapzit = tlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
}
ulap = g1_gpu[__I2_G1(i,j)]*ulapcsi+ucsi*dg1+g2_gpu[__I2_G2(i,j)]*ulapeta+ueta*dg2;
vlap = g1_gpu[__I2_G1(i,j)]*vlapcsi+vcsi*dg1+g2_gpu[__I2_G2(i,j)]*vlapeta+veta*dg2;
wlap = g1_gpu[__I2_G1(i,j)]*wlapcsi+wcsi*dg1+g2_gpu[__I2_G2(i,j)]*wlapeta+weta*dg2;
tlap = g1_gpu[__I2_G1(i,j)]*tlapcsi+tcsi*dg1+g2_gpu[__I2_G2(i,j)]*tlapeta+teta*dg2;
ulap = ulap * jac_gpu[__I2_JAC(i,j)];
vlap = vlap * jac_gpu[__I2_JAC(i,j)];
wlap = wlap * jac_gpu[__I2_JAC(i,j)];
tlap = tlap * jac_gpu[__I2_JAC(i,j)];
ulapz = ulapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+uz*dzitdz2_gpu[__I1_DZITDZ2(k)];
vlapz = vlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+vz*dzitdz2_gpu[__I1_DZITDZ2(k)];
wlapz = wlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+wz*dzitdz2_gpu[__I1_DZITDZ2(k)];
tlapz = tlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+tz*dzitdz2_gpu[__I1_DZITDZ2(k)];
ulap = ulap+ulapz;
vlap = vlap+vlapz;
wlap = wlap+wlapz;
tlap = tlap+tlapz;
div = ux+vy+wz;
div3l = div/3.0;
w_aux_gpu[__I4_W_AUX(i,j,k,10)] = div3l;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
w_aux_gpu[__I4_W_AUX(i,j,k,9)] = sqrt(omod2);
w_aux_gpu[__I4_W_AUX(i,j,k,8)]=max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0);
if( j == 1) {
if(iblock == ite_rank_x && i == ite_l) w_aux_gpu[__I4_W_AUX(i,j,k,8)] = teshk;
if(iblock == itu_rank_x && i == itu_l) w_aux_gpu[__I4_W_AUX(i,j,k,8)] = teshk;
}
sig11 = 2.0*(ux-div3l);
sig12 = uy+vx;
sig13 = uz+wx;
sig22 = 2.0*(vy-div3l);
sig23 = vz+wy;
sig33 = 2.0*(wz-div3l);
sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap;
sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap;
sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap;
sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl;
sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu;
sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int lmax,real prandtl,real u0,real l0,real t0,real teshk,int *vis_tag_gpu,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dzitdz2_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *g1_gpu,real *g2_gpu,real *g12_gpu,real *jac_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,teshk,vis_tag_gpu,wall_tag_gpu,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu);
}
}



__global__ void  visflx_nosensor_kernel(int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int lmax,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dcsidx2_gpu,real *detady2_gpu,real *dzitdz2_gpu,real *cp_coeff_gpu){
//Kernel for visflx_nosensor_kernel
real ccl;real clapl;real sig11;
real sig12;real sig13;real sig22;
real sig23;real sig33;real uu;
real vv;real ww;real tt;
real mu;real ux;real uy;
real uz;real vx;real vy;
real vz;real wx;real wy;
real wz;real tx;real ty;
real tz;real mux;real muy;
real muz;real ulap;real ulapx;
real ulapy;real ulapz;real vlap;
real vlapx;real vlapy;real vlapz;
real wlap;real wlapx;real wlapy;
real wlapz;real tlap;real tlapx;
real tlapy;real tlapz;real sigq;
real sigx;real sigy;real sigz;
real sigqt;real sigah;real div;
real div3l;real omegax;real omegay;
real omegaz;real omod2;real cploc;
int i;int j;int k;
int l;int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
ux = 0.0;
vx = 0.0;
wx = 0.0;
tx = 0.0;
mux = 0.0;
uy = 0.0;
vy = 0.0;
wy = 0.0;
ty = 0.0;
muy = 0.0;
uz = 0.0;
vz = 0.0;
wz = 0.0;
tz = 0.0;
muz = 0.0;
ulapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*uu;
ulapy = ulapx;
ulapz = ulapx;
vlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*vv;
vlapy = vlapx;
vlapz = vlapx;
wlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*ww;
wlapy = wlapx;
wlapz = wlapx;
tlapx = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*tt;
tlapy = tlapx;
tlapz = tlapx;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
clapl = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax)];
ux = ux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
ulapx = ulapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vx = vx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
vlapx = vlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wx = wx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
wlapx = wlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
tx = tx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
tlapx = tlapx + clapl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
mux = mux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,7)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,7)]);
}
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
clapl = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax)];
uy = uy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
ulapy = ulapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
vy = vy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
vlapy = vlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
wy = wy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
wlapy = wlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
ty = ty+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
tlapy = tlapy + clapl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
muy = muy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,7)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,7)]);
}
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
clapl = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax)];
uz = uz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
ulapz = ulapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vz = vz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
vlapz = vlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wz = wz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
wlapz = wlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tz = tz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
tlapz = tlapz + clapl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
muz = muz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,7)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,7)]);
}
ux = ux*dcsidx_gpu[__I1_DCSIDX(i)];
vx = vx*dcsidx_gpu[__I1_DCSIDX(i)];
wx = wx*dcsidx_gpu[__I1_DCSIDX(i)];
tx = tx*dcsidx_gpu[__I1_DCSIDX(i)];
mux = mux*dcsidx_gpu[__I1_DCSIDX(i)];
uy = uy*detady_gpu[__I1_DETADY(j)];
vy = vy*detady_gpu[__I1_DETADY(j)];
wy = wy*detady_gpu[__I1_DETADY(j)];
ty = ty*detady_gpu[__I1_DETADY(j)];
muy = muy*detady_gpu[__I1_DETADY(j)];
uz = uz*dzitdz_gpu[__I1_DZITDZ(k)];
vz = vz*dzitdz_gpu[__I1_DZITDZ(k)];
wz = wz*dzitdz_gpu[__I1_DZITDZ(k)];
tz = tz*dzitdz_gpu[__I1_DZITDZ(k)];
muz = muz*dzitdz_gpu[__I1_DZITDZ(k)];
ulapx = ulapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+ux*dcsidx2_gpu[__I1_DCSIDX2(i)];
vlapx = vlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+vx*dcsidx2_gpu[__I1_DCSIDX2(i)];
wlapx = wlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+wx*dcsidx2_gpu[__I1_DCSIDX2(i)];
tlapx = tlapx*dcsidxs_gpu[__I1_DCSIDXS(i)]+tx*dcsidx2_gpu[__I1_DCSIDX2(i)];
ulapy = ulapy*detadys_gpu[__I1_DETADYS(j)]+uy*detady2_gpu[__I1_DETADY2(j)];
vlapy = vlapy*detadys_gpu[__I1_DETADYS(j)]+vy*detady2_gpu[__I1_DETADY2(j)];
wlapy = wlapy*detadys_gpu[__I1_DETADYS(j)]+wy*detady2_gpu[__I1_DETADY2(j)];
tlapy = tlapy*detadys_gpu[__I1_DETADYS(j)]+ty*detady2_gpu[__I1_DETADY2(j)];
ulapz = ulapz*dzitdzs_gpu[__I1_DZITDZS(k)]+uz*dzitdz2_gpu[__I1_DZITDZ2(k)];
vlapz = vlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+vz*dzitdz2_gpu[__I1_DZITDZ2(k)];
wlapz = wlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+wz*dzitdz2_gpu[__I1_DZITDZ2(k)];
tlapz = tlapz*dzitdzs_gpu[__I1_DZITDZS(k)]+tz*dzitdz2_gpu[__I1_DZITDZ2(k)];
ulap = ulapx+ulapy+ulapz;
vlap = vlapx+vlapy+vlapz;
wlap = wlapx+wlapy+wlapz;
tlap = tlapx+tlapy+tlapz;
div = ux+vy+wz;
div3l = div/3.0;
w_aux_gpu[__I4_W_AUX(i,j,k,10)] = div3l;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
w_aux_gpu[__I4_W_AUX(i,j,k,9)] = sqrt(omod2);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
}
sig11 = 2.0*(ux-div3l);
sig12 = uy+vx;
sig13 = uz+wx;
sig22 = 2.0*(vy-div3l);
sig23 = vz+wy;
sig33 = 2.0*(wz-div3l);
sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap;
sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap;
sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap;
sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl;
sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu;
sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;
if (j==1) {
wallprop_gpu[__I3_WALLPROP(i,k,2)] = mu*uy;
wallprop_gpu[__I3_WALLPROP(i,k,3)] = mu*wy;
wallprop_gpu[__I3_WALLPROP(i,k,4)] = mu*ty*cploc/prandtl;
}

}
}


extern "C"{
void visflx_nosensor_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int lmax,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dcsidx2_gpu,real *detady2_gpu,real *dzitdz2_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_nosensor_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,lmax,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dcsidx2_gpu,detady2_gpu,dzitdz2_gpu,cp_coeff_gpu);
}
}



__global__ void  visflx_nosensor_c2_kernel(int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int ortho,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int lmax,real prandtl,real u0,real l0,real t0,int *vis_tag_gpu,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dzitdz2_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *g1_gpu,real *g2_gpu,real *g12_gpu,real *jac_gpu,real *cp_coeff_gpu){
//Kernel for visflx_nosensor_c2_kernel
real ccl;real clapl;real cli;
real clj;real clk;real clapi;
real clapj;real clapk;real sig11;
real sig12;real sig13;real sig22;
real sig23;real sig33;real uu;
real vv;real ww;real tt;
real mu;real ux;real uy;
real uz;real vx;real vy;
real vz;real wx;real wy;
real wz;real tx;real ty;
real tz;real mux;real muy;
real muz;real ucsi;real ueta;
real uzit;real vcsi;real veta;
real vzit;real wcsi;real weta;
real wzit;real tcsi;real teta;
real tzit;real mucsi;real mueta;
real muzit;real ulapcsi;real ulapeta;
real ulapzit;real vlapcsi;real vlapeta;
real vlapzit;real wlapcsi;real wlapeta;
real wlapzit;real ulapcsieta;real vlapcsieta;
real wlapcsieta;real tlapcsieta;real dg1;
real dg2;real dg12csi;real dg12eta;
real tlapcsi;real tlapeta;real tlapzit;
real ulap;real ulapx;real ulapy;
real ulapz;real vlap;real vlapx;
real vlapy;real vlapz;real wlap;
real wlapx;real wlapy;real wlapz;
real tlap;real tlapx;real tlapy;
real tlapz;real sigq;real sigx;
real sigy;real sigz;real sigqt;
real sigah;real div;real div3l;
real omegax;real omegay;real omegaz;
real omod2;real cploc;
int i;int j;int k;
int l;int ll;int iercuda;
int lmaxi;int lmaxj;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
lmaxi = lmax;
lmaxj = lmax;
if (j == 1) lmaxi = vis_tag_gpu[__I1_VIS_TAG(i)];
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) lmaxj = min(j,lmax);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
}
ucsi = 0.0;
vcsi = 0.0;
wcsi = 0.0;
tcsi = 0.0;
mucsi = 0.0;
ueta = 0.0;
veta = 0.0;
weta = 0.0;
teta = 0.0;
mueta = 0.0;
uzit = 0.0;
vzit = 0.0;
wzit = 0.0;
tzit = 0.0;
muzit = 0.0;
dg1 = 0.0;
dg2 = 0.0;
dg12csi = 0.0;
dg12eta = 0.0;
for(int l=1; l<lmax+1; l++){
cli = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxi)];
clj = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxj)];
clk = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax )];
ucsi = ucsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vcsi = vcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wcsi = wcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
tcsi = tcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
mucsi = mucsi+cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,7)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,7)]);
ueta = ueta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
veta = veta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
weta = weta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
teta = teta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
mueta = mueta+clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,7)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,7)]);
uzit = uzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vzit = vzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wzit = wzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tzit = tzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
muzit = muzit+clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,7)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,7)]);
dg1 = dg1 +cli*(g1_gpu[__I2_G1(i+l,j)]-g1_gpu[__I2_G1(i-l,j)]);
dg12csi=dg12csi+cli*(g12_gpu[__I2_G12(i+l,j)]-g12_gpu[__I2_G12(i-l,j)]);
dg2 = dg2 +clj*(g2_gpu[__I2_G2(i,j+l)]-g2_gpu[__I2_G2(i,j-l)]);
dg12eta=dg12eta+clj*(g12_gpu[__I2_G12(i,j+l)]-g12_gpu[__I2_G12(i,j-l)]);
}
ux = ucsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + ueta *detadxc2_gpu[__I2_DETADXC2(i,j)];
vx = vcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + veta *detadxc2_gpu[__I2_DETADXC2(i,j)];
wx = wcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + weta *detadxc2_gpu[__I2_DETADXC2(i,j)];
tx = tcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + teta *detadxc2_gpu[__I2_DETADXC2(i,j)];
mux = mucsi*dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + mueta*detadxc2_gpu[__I2_DETADXC2(i,j)];
uy = ucsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + ueta *detadyc2_gpu[__I2_DETADYC2(i,j)];
vy = vcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + veta *detadyc2_gpu[__I2_DETADYC2(i,j)];
wy = wcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + weta *detadyc2_gpu[__I2_DETADYC2(i,j)];
ty = tcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + teta *detadyc2_gpu[__I2_DETADYC2(i,j)];
muy = mucsi*dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + mueta*detadyc2_gpu[__I2_DETADYC2(i,j)];
uz = uzit *dzitdz_gpu[__I1_DZITDZ(k)];
vz = vzit *dzitdz_gpu[__I1_DZITDZ(k)];
wz = wzit *dzitdz_gpu[__I1_DZITDZ(k)];
tz = tzit *dzitdz_gpu[__I1_DZITDZ(k)];
muz = muzit*dzitdz_gpu[__I1_DZITDZ(k)];
if (j==1) {
wallprop_gpu[__I3_WALLPROP(i,k,3)] = mu*weta;
wallprop_gpu[__I3_WALLPROP(i,k,4)] = mu*teta*cploc/prandtl;
}
ulapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*uu;
ulapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*uu;
ulapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*uu;
vlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*vv;
vlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*vv;
vlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*vv;
wlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*ww;
wlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*ww;
wlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*ww;
tlapcsi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxi)]*tt;
tlapeta = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmaxj)]*tt;
tlapzit = coeff_deriv2_gpu[__I2_COEFF_DERIV2(0,lmax)]*tt;
for(int l=1; l<lmax+1; l++){
clapi = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmaxi)];
clapj = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmaxj)];
clapk = coeff_deriv2_gpu[__I2_COEFF_DERIV2(l,lmax )];
ulapcsi = ulapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
ulapeta = ulapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
ulapzit = ulapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vlapcsi = vlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
vlapeta = vlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
vlapzit = vlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wlapcsi = wlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
wlapeta = wlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
wlapzit = wlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
tlapcsi = tlapcsi + clapi*(w_aux_gpu[__I4_W_AUX(i+l,j,k,6)]+w_aux_gpu[__I4_W_AUX(i-l,j,k,6)]);
tlapeta = tlapeta + clapj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,6)]+w_aux_gpu[__I4_W_AUX(i,j-l,k,6)]);
tlapzit = tlapzit + clapk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,6)]+w_aux_gpu[__I4_W_AUX(i,j,k-l,6)]);
}
if(ortho == 1) {
ulap = g1_gpu[__I2_G1(i,j)]*ulapcsi+ucsi*dg1+g2_gpu[__I2_G2(i,j)]*ulapeta+ueta*dg2;
vlap = g1_gpu[__I2_G1(i,j)]*vlapcsi+vcsi*dg1+g2_gpu[__I2_G2(i,j)]*vlapeta+veta*dg2;
wlap = g1_gpu[__I2_G1(i,j)]*wlapcsi+wcsi*dg1+g2_gpu[__I2_G2(i,j)]*wlapeta+weta*dg2;
tlap = g1_gpu[__I2_G1(i,j)]*tlapcsi+tcsi*dg1+g2_gpu[__I2_G2(i,j)]*tlapeta+teta*dg2;
}else {
ulapcsieta = 0.250*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,2)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,2)]);
vlapcsieta = 0.250*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,3)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,3)]);
wlapcsieta = 0.250*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,4)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,4)]);
tlapcsieta = 0.250*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,6)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,6)]);
if (iblock==ite_rank_x&&i==ite_l-1&&j==1) {
ulapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,2)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,2)]);
vlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,3)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,3)]);
wlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,4)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,4)]);
tlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i-1,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,6)]+w_aux_gpu[__I4_W_AUX(i-1,j-1,k,6)]);
}
if (iblock==itu_rank_x&&i==itu_l+1&&j==1) {
ulapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,2)]+w_aux_gpu[__I4_W_AUX(i,j-1,k,2)]);
vlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,3)]+w_aux_gpu[__I4_W_AUX(i,j-1,k,3)]);
wlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,4)]+w_aux_gpu[__I4_W_AUX(i,j-1,k,4)]);
tlapcsieta = 0.50*(w_aux_gpu[__I4_W_AUX(i+1,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i+1,j-1,k,6)]+w_aux_gpu[__I4_W_AUX(i,j-1,k,6)]);
}
ulap = g1_gpu[__I2_G1(i,j)]*ulapcsi+ucsi*dg1+2.0*g12_gpu[__I2_G12(i,j)]*ulapcsieta+dg12csi*ueta+dg12eta*ucsi+g2_gpu[__I2_G2(i,j)]*ulapeta+ueta*dg2;
vlap = g1_gpu[__I2_G1(i,j)]*vlapcsi+vcsi*dg1+2.0*g12_gpu[__I2_G12(i,j)]*vlapcsieta+dg12csi*veta+dg12eta*vcsi+g2_gpu[__I2_G2(i,j)]*vlapeta+veta*dg2;
wlap = g1_gpu[__I2_G1(i,j)]*wlapcsi+wcsi*dg1+2.0*g12_gpu[__I2_G12(i,j)]*wlapcsieta+dg12csi*weta+dg12eta*wcsi+g2_gpu[__I2_G2(i,j)]*wlapeta+weta*dg2;
tlap = g1_gpu[__I2_G1(i,j)]*tlapcsi+tcsi*dg1+2.0*g12_gpu[__I2_G12(i,j)]*tlapcsieta+dg12csi*teta+dg12eta*tcsi+g2_gpu[__I2_G2(i,j)]*tlapeta+teta*dg2;
}
ulap = ulap * jac_gpu[__I2_JAC(i,j)];
vlap = vlap * jac_gpu[__I2_JAC(i,j)];
wlap = wlap * jac_gpu[__I2_JAC(i,j)];
tlap = tlap * jac_gpu[__I2_JAC(i,j)];
ulapz = ulapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+uz*dzitdz2_gpu[__I1_DZITDZ2(k)];
vlapz = vlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+vz*dzitdz2_gpu[__I1_DZITDZ2(k)];
wlapz = wlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+wz*dzitdz2_gpu[__I1_DZITDZ2(k)];
tlapz = tlapzit*dzitdzs_gpu[__I1_DZITDZS(k)]+tz*dzitdz2_gpu[__I1_DZITDZ2(k)];
ulap = ulap+ulapz;
vlap = vlap+vlapz;
wlap = wlap+wlapz;
tlap = tlap+tlapz;
div = ux+vy+wz;
div3l = div/3.0;
w_aux_gpu[__I4_W_AUX(i,j,k,10)] = div3l;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
w_aux_gpu[__I4_W_AUX(i,j,k,9)] = sqrt(omod2);
sig11 = 2.0*(ux-div3l);
sig12 = uy+vx;
sig13 = uz+wx;
sig22 = 2.0*(vy-div3l);
sig23 = vz+wy;
sig33 = 2.0*(wz-div3l);
sigx = mux*sig11 + muy*sig12 + muz*sig13 + mu*ulap;
sigy = mux*sig12 + muy*sig22 + muz*sig23 + mu*vlap;
sigz = mux*sig13 + muy*sig23 + muz*sig33 + mu*wlap;
sigqt = (mux*tx+muy*ty+muz*tz+mu*tlap)*cploc/prandtl;
sigah = (sig11*ux+sig12*uy+sig13*uz+sig12*vx+sig22*vy+sig23*vz+sig13*wx+sig23*wy+sig33*wz)*mu;
sigq = sigx*uu+sigy*vv+sigz*ww+sigah+sigqt;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_nosensor_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int calorically_perfect,int indx_cp_l,int indx_cp_r,int ortho,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int lmax,real prandtl,real u0,real l0,real t0,int *vis_tag_gpu,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *coeff_deriv2_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dzitdz2_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *g1_gpu,real *g2_gpu,real *g12_gpu,real *jac_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_nosensor_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,calorically_perfect,indx_cp_l,indx_cp_r,ortho,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,lmax,prandtl,u0,l0,t0,vis_tag_gpu,wall_tag_gpu,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,coeff_deriv1_gpu,coeff_deriv2_gpu,dzitdz_gpu,dzitdzs_gpu,dzitdz2_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,g1_gpu,g2_gpu,g12_gpu,jac_gpu,cp_coeff_gpu);
}
}



__global__ void  sponge_kernel(int nx,int ny,int nz,int ng,int nv,int j_sponge,real *w_gpu,real *fln_gpu,real *wfar_gpu,real *f_sponge_gpu){
//Kernel for sponge_kernel
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,j_sponge);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
fln_gpu[__I4_FLN(i,j,k,m)] = fln_gpu[__I4_FLN(i,j,k,m)] - f_sponge_gpu[__I1_F_SPONGE(j)]*(w_gpu[__I4_W(i,j,k,m)] - wfar_gpu[__I2_WFAR(i,m)]);
}

}
}


extern "C"{
void sponge_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int j_sponge,real *w_gpu,real *fln_gpu,real *wfar_gpu,real *f_sponge_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(j_sponge)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((sponge_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,j_sponge,w_gpu,fln_gpu,wfar_gpu,f_sponge_gpu);
}
}

__device__ real get_temperature_from_e_dev_limiter_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_limiter_kernel1_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}





__global__ void  limiter_kernel1_n_limited_rho(int nx,int ny,int nz,int ng,int iblock,int kblock,int calorically_perfect,int indx_cp_l,int indx_cp_r,real t0,real tol_iter_nr,real rho_lim,real tem_lim,real rho_lim_rescale,real tem_lim_rescale,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *n_limited_rho,real *n_limited_tem,real *redn_3d_gpu){
//Kernel for limiter_kernel1_n_limited_rho
real rho;real tem;real eei;
real uu;real vv;real ww;
real rhoe;real qq;real ee;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
rho = w_gpu[__I4_W(i,j,k,1)];
uu = w_gpu[__I4_W(i,j,k,2)]/rho;
vv = w_gpu[__I4_W(i,j,k,3)]/rho;
ww = w_gpu[__I4_W(i,j,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,j,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tem = get_temperature_from_e_dev_limiter_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
if((rho < rho_lim)||(tem < tem_lim)) {
if(rho < rho_lim) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1.0;
}
if(tem < tem_lim) {
}
}

}
}


__global__ void  limiter_kernel1_n_limited_tem(int nx,int ny,int nz,int ng,int iblock,int kblock,int calorically_perfect,int indx_cp_l,int indx_cp_r,real t0,real tol_iter_nr,real rho_lim,real tem_lim,real rho_lim_rescale,real tem_lim_rescale,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *n_limited_rho,real *n_limited_tem,real *redn_3d_gpu){
//Kernel for limiter_kernel1_n_limited_tem
real rho;real tem;real eei;
real uu;real vv;real ww;
real rhoe;real qq;real ee;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
rho = w_gpu[__I4_W(i,j,k,1)];
uu = w_gpu[__I4_W(i,j,k,2)]/rho;
vv = w_gpu[__I4_W(i,j,k,3)]/rho;
ww = w_gpu[__I4_W(i,j,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,j,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tem = get_temperature_from_e_dev_limiter_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
if((rho < rho_lim)||(tem < tem_lim)) {
if(rho < rho_lim) {
}
if(tem < tem_lim) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +1.0;
}
}

}
}



extern "C"{
void limiter_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int iblock,int kblock,int calorically_perfect,int indx_cp_l,int indx_cp_r,real t0,real tol_iter_nr,real rho_lim,real tem_lim,real rho_lim_rescale,real tem_lim_rescale,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *n_limited_rho,real *n_limited_tem,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((limiter_kernel1_n_limited_rho),grid,block,0,stream,nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,w_gpu,w_aux_gpu,cv_coeff_gpu,n_limited_rho,n_limited_tem,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, n_limited_rho);

hipLaunchKernelGGL((limiter_kernel1_n_limited_tem),grid,block,0,stream,nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,w_gpu,w_aux_gpu,cv_coeff_gpu,n_limited_rho,n_limited_tem,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, n_limited_tem);


}
}
__device__ real get_temperature_from_e_dev_limiter_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_limiter_kernel2_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_limiter_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_limiter_kernel2_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  limiter_kernel2(int nx,int ny,int nz,int ng,int iblock,int kblock,int calorically_perfect,int indx_cp_l,int indx_cp_r,real t0,real tol_iter_nr,real rho_lim,real tem_lim,real rho_lim_rescale,real tem_lim_rescale,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
//Kernel for limiter_kernel2
real rho;real tem;real eei;
real uu;real vv;real ww;
real rhoe;real qq;real ee;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(i,j,k,1)];
uu = w_gpu[__I4_W(i,j,k,2)]/rho;
vv = w_gpu[__I4_W(i,j,k,3)]/rho;
ww = w_gpu[__I4_W(i,j,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,j,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tem = get_temperature_from_e_dev_limiter_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
if((rho < rho_lim)||(tem < tem_lim)) {
if(rho < rho_lim) {
rho = rho_lim*rho_lim_rescale;
printf("limiter rhofix:  %d,%d,%d,%d,%d : ",iblock,kblock,i,j,k);
}
if(tem < tem_lim) {
tem = tem_lim*tem_lim_rescale;
printf("limiter temfix:  %d,%d,%d,%d,%d : ",iblock,kblock,i,j,k);
}
w_gpu[__I4_W(i,j,k,1)] = rho;
w_gpu[__I4_W(i,j,k,2)] = rho*uu;
w_gpu[__I4_W(i,j,k,3)] = rho*vv;
w_gpu[__I4_W(i,j,k,4)] = rho*ww;
eei = get_e_from_temperature_dev_limiter_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tem,cv_coeff_gpu);
w_gpu[__I4_W(i,j,k,5)]=rho*eei+0.50*(((w_gpu[__I4_W(i,j,k,2)])*(w_gpu[__I4_W(i,j,k,2)]))+((w_gpu[__I4_W(i,j,k,3)])*(w_gpu[__I4_W(i,j,k,3)]))+((w_gpu[__I4_W(i,j,k,4)])*(w_gpu[__I4_W(i,j,k,4)])))/rho;
}

}
}


extern "C"{
void limiter_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int iblock,int kblock,int calorically_perfect,int indx_cp_l,int indx_cp_r,real t0,real tol_iter_nr,real rho_lim,real tem_lim,real rho_lim_rescale,real tem_lim_rescale,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((limiter_kernel2),grid,block,0,stream,nx,ny,nz,ng,iblock,kblock,calorically_perfect,indx_cp_l,indx_cp_r,t0,tol_iter_nr,rho_lim,tem_lim,rho_lim_rescale,tem_lim_rescale,w_gpu,w_aux_gpu,cv_coeff_gpu);
}
}

__device__ real get_temperature_from_e_dev_filter_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_filter_kernel1_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_temperature_from_e_dev_filter_kernel1_1(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_filter_kernel1_1
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_filter_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_filter_kernel1_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  filter_kernel1(int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int jfilter,real t0,real tol_iter_nr,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *coeff_filter_gpu){
//Kernel for filter_kernel1
real cfilt;real rho;real eei;
real uu;real vv;real ww;
real rhoe;real qq;real ee;
real rho_p;real rho_m;real rho_f;
real tem_p;real tem_m;real tem_f;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1+4);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny-4,1)&&loop_cond(k,nz,1)){
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1 && j < jfilter) {
}else {
rho_f = 0.0;
tem_f = 0.0;
for(int l=0; l<4+1; l++){
rho_p = w_gpu[__I4_W(i+l,j,k,1)];
uu = w_gpu[__I4_W(i+l,j,k,2)]/rho_p;
vv = w_gpu[__I4_W(i+l,j,k,3)]/rho_p;
ww = w_gpu[__I4_W(i+l,j,k,4)]/rho_p;
rhoe = w_gpu[__I4_W(i+l,j,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho_p-qq;
tem_p = get_temperature_from_e_dev_filter_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i+l,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
rho_m = w_gpu[__I4_W(i-l,j,k,1)];
uu = w_gpu[__I4_W(i-l,j,k,2)]/rho_m;
vv = w_gpu[__I4_W(i-l,j,k,3)]/rho_m;
ww = w_gpu[__I4_W(i-l,j,k,4)]/rho_m;
rhoe = w_gpu[__I4_W(i-l,j,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho_m-qq;
tem_m = get_temperature_from_e_dev_filter_kernel1_1(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i-l,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
cfilt = coeff_filter_gpu[__I1_COEFF_FILTER(l)];
rho_f = rho_f + cfilt*(rho_p+rho_m);
tem_f = tem_f + cfilt*(tem_p+tem_m);
}
uu = w_gpu[__I4_W(i,j,k,2)]/w_gpu[__I4_W(i,j,k,1)];
vv = w_gpu[__I4_W(i,j,k,3)]/w_gpu[__I4_W(i,j,k,1)];
ww = w_gpu[__I4_W(i,j,k,4)]/w_gpu[__I4_W(i,j,k,1)];
w_gpu[__I4_W(i,j,k,1)] = rho_f;
w_gpu[__I4_W(i,j,k,2)] = rho_f*uu;
w_gpu[__I4_W(i,j,k,3)] = rho_f*vv;
w_gpu[__I4_W(i,j,k,4)] = rho_f*ww;
eei = get_e_from_temperature_dev_filter_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tem_f,cv_coeff_gpu);
w_gpu[__I4_W(i,j,k,5)]=rho_f*eei+0.50*(((w_gpu[__I4_W(i,j,k,2)])*(w_gpu[__I4_W(i,j,k,2)]))+((w_gpu[__I4_W(i,j,k,3)])*(w_gpu[__I4_W(i,j,k,3)]))+((w_gpu[__I4_W(i,j,k,4)])*(w_gpu[__I4_W(i,j,k,4)])))/rho_f;
w_aux_gpu[__I4_W_AUX(i,j,k,6)] = tem_f;
}

}
}


extern "C"{
void filter_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int jfilter,real t0,real tol_iter_nr,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *coeff_filter_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny-4)-(1+4)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((filter_kernel1),grid,block,0,stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter,t0,tol_iter_nr,wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu);
}
}
__device__ real get_temperature_from_e_dev_filter_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_filter_kernel2_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_temperature_from_e_dev_filter_kernel2_1(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_filter_kernel2_1
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_filter_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_filter_kernel2_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  filter_kernel2(int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int jfilter,real t0,real tol_iter_nr,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *coeff_filter_gpu){
//Kernel for filter_kernel2
real cfilt;real rho;real eei;
real uu;real vv;real ww;
real rhoe;real qq;real ee;
real rho_p;real rho_m;real rho_f;
real tem_p;real tem_m;real tem_f;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1+4);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny-4,1)&&loop_cond(k,nz,1)){
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1 && j < jfilter) {
}else {
rho_f = 0.0;
tem_f = 0.0;
for(int l=0; l<4+1; l++){
rho_p = w_gpu[__I4_W(i,j+l,k,1)];
uu = w_gpu[__I4_W(i,j+l,k,2)]/rho_p;
vv = w_gpu[__I4_W(i,j+l,k,3)]/rho_p;
ww = w_gpu[__I4_W(i,j+l,k,4)]/rho_p;
rhoe = w_gpu[__I4_W(i,j+l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho_p-qq;
tem_p = get_temperature_from_e_dev_filter_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i+l,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
rho_m = w_gpu[__I4_W(i,j-l,k,1)];
uu = w_gpu[__I4_W(i,j-l,k,2)]/rho_m;
vv = w_gpu[__I4_W(i,j-l,k,3)]/rho_m;
ww = w_gpu[__I4_W(i,j-l,k,4)]/rho_m;
rhoe = w_gpu[__I4_W(i,j-l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho_m-qq;
tem_m = get_temperature_from_e_dev_filter_kernel2_1(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i-l,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
cfilt = coeff_filter_gpu[__I1_COEFF_FILTER(l)];
rho_f = rho_f + cfilt*(rho_p+rho_m);
tem_f = tem_f + cfilt*(tem_p+tem_m);
}
uu = w_gpu[__I4_W(i,j,k,2)]/w_gpu[__I4_W(i,j,k,1)];
vv = w_gpu[__I4_W(i,j,k,3)]/w_gpu[__I4_W(i,j,k,1)];
ww = w_gpu[__I4_W(i,j,k,4)]/w_gpu[__I4_W(i,j,k,1)];
w_gpu[__I4_W(i,j,k,1)] = rho_f;
w_gpu[__I4_W(i,j,k,2)] = rho_f*uu;
w_gpu[__I4_W(i,j,k,3)] = rho_f*vv;
w_gpu[__I4_W(i,j,k,4)] = rho_f*ww;
eei = get_e_from_temperature_dev_filter_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tem_f,cv_coeff_gpu);
w_gpu[__I4_W(i,j,k,5)]=rho_f*eei+0.50*(((w_gpu[__I4_W(i,j,k,2)])*(w_gpu[__I4_W(i,j,k,2)]))+((w_gpu[__I4_W(i,j,k,3)])*(w_gpu[__I4_W(i,j,k,3)]))+((w_gpu[__I4_W(i,j,k,4)])*(w_gpu[__I4_W(i,j,k,4)])))/rho_f;
w_aux_gpu[__I4_W_AUX(i,j,k,6)] = tem_f;
}

}
}


extern "C"{
void filter_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int jfilter,real t0,real tol_iter_nr,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *coeff_filter_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny-4)-(1+4)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((filter_kernel2),grid,block,0,stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,jfilter,t0,tol_iter_nr,wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,coeff_filter_gpu);
}
}



__global__ void  visflx_reduced_ord2_kernel(int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int update_sensor,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *x_gpu,real *y_gpu,real *z_gpu,real *cp_coeff_gpu){
//Kernel for visflx_reduced_ord2_kernel
real dxl;real dyl;real dzl;
real sig11;real sig12;real sig13;
real sig21;real sig22;real sig23;
real sig31;real sig32;real sig33;
real uu;real vv;real ww;
real tt;real mu;real ux;
real uy;real uz;real vx;
real vy;real vz;real wx;
real wy;real wz;real tx;
real ty;real tz;real mux;
real muy;real muz;real sigq;
real sigx;real sigy;real sigz;
real sigqt;real sigah;real div;
real div3l;real omegax;real omegay;
real omegaz;real omod2;real cploc;
int i;int j;int k;
int l;int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
}
dxl = 1.0/(x_gpu[__I1_X(i+1)]-x_gpu[__I1_X(i-1)]);
dyl = 1.0/(y_gpu[__I1_Y(j+1)]-y_gpu[__I1_Y(j-1)]);
dzl = 1.0/(z_gpu[__I1_Z(k+1)]-z_gpu[__I1_Z(k-1)]);
ux = dxl*(w_aux_gpu[__I4_W_AUX(i+1,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-1,j,k,2)]);
vx = dxl*(w_aux_gpu[__I4_W_AUX(i+1,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-1,j,k,3)]);
wx = dxl*(w_aux_gpu[__I4_W_AUX(i+1,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-1,j,k,4)]);
mux = dxl*(w_aux_gpu[__I4_W_AUX(i+1,j,k,7)]-w_aux_gpu[__I4_W_AUX(i-1,j,k,7)]);
uy = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,2)]);
vy = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,3)]);
wy = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,4)]);
ty = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,6)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,6)]);
muy = dyl*(w_aux_gpu[__I4_W_AUX(i,j+1,k,7)]-w_aux_gpu[__I4_W_AUX(i,j-1,k,7)]);
uz = dzl*(w_aux_gpu[__I4_W_AUX(i,j,k+1,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-1,2)]);
vz = dzl*(w_aux_gpu[__I4_W_AUX(i,j,k+1,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-1,3)]);
wz = dzl*(w_aux_gpu[__I4_W_AUX(i,j,k+1,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-1,4)]);
muz = dzl*(w_aux_gpu[__I4_W_AUX(i,j,k+1,7)]-w_aux_gpu[__I4_W_AUX(i,j,k-1,7)]);
if (j==1) {
wallprop_gpu[__I3_WALLPROP(i,k,2)] = mu*uy;
wallprop_gpu[__I3_WALLPROP(i,k,3)] = mu*wy;
wallprop_gpu[__I3_WALLPROP(i,k,4)] = mu*ty*cploc/prandtl;
}
div = ux+vy+wz;
div3l = div/3.0;
w_aux_gpu[__I4_W_AUX(i,j,k,10)] = div3l;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
w_aux_gpu[__I4_W_AUX(i,j,k,9)] = sqrt(omod2);
if (update_sensor == 1) {
w_aux_gpu[__I4_W_AUX(i,j,k,8)]=(((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0)))*((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0))));
}
sig11 = ux-2.0*div3l;
sig12 = vx;
sig13 = wx;
sig21 = uy;
sig22 = vy-2.0*div3l;
sig23 = wy;
sig31 = uz;
sig32 = vz;
sig33 = wz-2.0*div3l;
sigx = mux*sig11 + muy*sig12 + muz*sig13;
sigy = mux*sig21 + muy*sig22 + muz*sig23;
sigz = mux*sig31 + muy*sig32 + muz*sig33;
sigah = (sig11*ux+sig12*uy+sig13*uz+sig21*vx+sig22*vy+sig23*vz+sig31*wx+sig32*wy+sig33*wz)*mu;
sigq = sigx*uu+sigy*vv+sigz*ww+sigah;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - sigx;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - sigy;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] - sigz;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] - sigq;

}
}


extern "C"{
void visflx_reduced_ord2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,int update_sensor,real prandtl,real u0,real l0,real t0,real *w_gpu,real *w_aux_gpu,real *wallprop_gpu,real *fl_gpu,real *x_gpu,real *y_gpu,real *z_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_reduced_ord2_kernel),grid,block,0,stream,nx,ny,nz,ng,calorically_perfect,indx_cp_l,indx_cp_r,update_sensor,prandtl,u0,l0,t0,w_gpu,w_aux_gpu,wallprop_gpu,fl_gpu,x_gpu,y_gpu,z_gpu,cp_coeff_gpu);
}
}



__global__ void  sensor_kernel(int nx,int ny,int nz,int ng,real u0,real l0,real *w_aux_gpu){
//Kernel for sensor_kernel
real div;real omod2;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
omod2=((w_aux_gpu[__I4_W_AUX(i,j,k,9)])*(w_aux_gpu[__I4_W_AUX(i,j,k,9)]));
div = 3.0*w_aux_gpu[__I4_W_AUX(i,j,k,10)];
w_aux_gpu[__I4_W_AUX(i,j,k,8)]=(((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0)))*((max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0))));

}
}


extern "C"{
void sensor_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real u0,real l0,real *w_aux_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((sensor_kernel),grid,block,0,stream,nx,ny,nz,ng,u0,l0,w_aux_gpu);
}
}



__global__ void  sensor_c2_kernel(int nx,int ny,int nz,int ng,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int jweno,real u0,real l0_ducros,real teshk,real theta_threshold,real *w_aux_gpu,real *theta_ij_gpu){
//Kernel for sensor_c2_kernel
real div;real omod2;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
omod2=((w_aux_gpu[__I4_W_AUX(i,j,k,9)])*(w_aux_gpu[__I4_W_AUX(i,j,k,9)]));
div = 3.0*w_aux_gpu[__I4_W_AUX(i,j,k,10)];
w_aux_gpu[__I4_W_AUX(i,j,k,8)]=max(-div/sqrt(omod2+((div)*(div))+(((u0/l0_ducros))*((u0/l0_ducros)))),0.0);
if( j == 1) {
if(iblock == ite_rank_x && i == ite_l) w_aux_gpu[__I4_W_AUX(i,j,k,8)] = teshk;
if(iblock == itu_rank_x && i == itu_l) w_aux_gpu[__I4_W_AUX(i,j,k,8)] = teshk;
}
if(theta_ij_gpu[__I2_THETA_IJ(i,j)] > theta_threshold || j > jweno) {
w_aux_gpu[__I4_W_AUX(i,j,k,8)] = 1000.0;
}

}
}


extern "C"{
void sensor_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int iblock,int ite_rank_x,int itu_rank_x,int ite_l,int itu_l,int jweno,real u0,real l0_ducros,real teshk,real theta_threshold,real *w_aux_gpu,real *theta_ij_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((sensor_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,iblock,ite_rank_x,itu_rank_x,ite_l,itu_l,jweno,u0,l0_ducros,teshk,theta_threshold,w_aux_gpu,theta_ij_gpu);
}
}



__global__ void  visflx_x_kernel1(int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *fhat_gpu,real *x_gpu,real *cp_coeff_gpu){
//Kernel for visflx_x_kernel1
real uu;real vv;real ww;
real tt;real mu;real qq;
real uup;real vvp;real wwp;
real ttp;real mup;real qqp;
real sigq;real sigx;real sigy;
real sigz;real sigq_tt;real sigq_qq;
real dxhl;real fl2o;real fl3o;
real fl4o;real fl5o;real ttf;
real muf;real cploc;
int i;int j;int k;
int iv;int ll;int iercuda;

i = __GIDX(x,0);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i ,j,k,2)];
uup = w_aux_gpu[__I4_W_AUX(i+1,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i ,j,k,3)];
vvp = w_aux_gpu[__I4_W_AUX(i+1,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i ,j,k,4)];
wwp = w_aux_gpu[__I4_W_AUX(i+1,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i ,j,k,6)];
ttp = w_aux_gpu[__I4_W_AUX(i+1,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i ,j,k,7)];
mup = w_aux_gpu[__I4_W_AUX(i+1,j,k,7)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
qqp = 0.50*(uup*uup+vvp*vvp+wwp*wwp);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
ttf = 0.50*(tt+ttp);
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((ttf/t0),ll);
}
}
sigx = uup-uu;
sigy = vvp-vv;
sigz = wwp-ww;
sigq_tt = ttp-tt;
sigq_qq = qqp-qq;
muf = mu+mup;
muf = 0.50*muf/(x_gpu[__I1_X(i+1)]-x_gpu[__I1_X(i)]);
sigx = sigx*muf;
sigy = sigy*muf;
sigz = sigz*muf;
sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf;
fhat_gpu[__I4_FHAT(i,j,k,2)] = - sigx;
fhat_gpu[__I4_FHAT(i,j,k,3)] = - sigy;
fhat_gpu[__I4_FHAT(i,j,k,4)] = - sigz;
fhat_gpu[__I4_FHAT(i,j,k,5)] = - sigq;

}
}


extern "C"{
void visflx_x_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *fhat_gpu,real *x_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(0)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_x_kernel1),grid,block,0,stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu);
}
}


__global__ void  visflx_x_kernel2(int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *fhat_gpu,real *x_gpu,real *cp_coeff_gpu){
//Kernel for visflx_x_kernel2
real uu;real vv;real ww;
real tt;real mu;real qq;
real uup;real vvp;real wwp;
real ttp;real mup;real qqp;
real sigq;real sigx;real sigy;
real sigz;real sigq_tt;real sigq_qq;
real dxhl;real fl2o;real fl3o;
real fl4o;real fl5o;real ttf;
real muf;real cploc;
int i;int j;int k;
int iv;int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
dxhl = 2.0/(x_gpu[__I1_X(i+1)]-x_gpu[__I1_X(i-1)]);
for(int iv=2; iv<nv+1; iv++){
fl_gpu[__I4_FL(i,j,k,iv)] = fl_gpu[__I4_FL(i,j,k,iv)] + dxhl*(fhat_gpu[__I4_FHAT(i,j,k,iv)]-fhat_gpu[__I4_FHAT(i-1,j,k,iv)]);
}

}
}


extern "C"{
void visflx_x_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *fhat_gpu,real *x_gpu,real *cp_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((visflx_x_kernel2),grid,block,0,stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,fhat_gpu,x_gpu,cp_coeff_gpu);
}
}



__global__ void  visflx_y_kernel(int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *y_gpu,real *cp_coeff_gpu){
//Kernel for visflx_y_kernel
real uu;real vv;real ww;
real tt;real mu;real qq;
real uup;real vvp;real wwp;
real ttp;real mup;real qqp;
real sigq;real sigx;real sigy;
real sigz;real sigq_tt;real sigq_qq;
real dyhl;real fl2o;real fl3o;
real fl4o;real fl5o;real ttf;
real muf;real cploc;
int i;int j;int k;
int iv;int ll;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=0; j<ny+1; j++){
uu = w_aux_gpu[__I4_W_AUX(i,j ,k,2)];
uup = w_aux_gpu[__I4_W_AUX(i,j+1,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j ,k,3)];
vvp = w_aux_gpu[__I4_W_AUX(i,j+1,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j ,k,4)];
wwp = w_aux_gpu[__I4_W_AUX(i,j+1,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j ,k,6)];
ttp = w_aux_gpu[__I4_W_AUX(i,j+1,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j ,k,7)];
mup = w_aux_gpu[__I4_W_AUX(i,j+1,k,7)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
qqp = 0.50*(uup*uup+vvp*vvp+wwp*wwp);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
ttf = 0.50*(tt+ttp);
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((ttf/t0),ll);
}
}
sigx = uup-uu;
sigy = vvp-vv;
sigz = wwp-ww;
sigq_tt = ttp-tt;
sigq_qq = qqp-qq;
muf = mu+mup;
muf = 0.50*muf/(y_gpu[__I1_Y(j+1)]-y_gpu[__I1_Y(j)]);
sigx = sigx*muf;
sigy = sigy*muf;
sigz = sigz*muf;
sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf;
if (j>0) {
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] + fl2o-sigx*dyhl;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] + fl3o-sigy*dyhl;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] + fl4o-sigz*dyhl;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] + fl5o-sigq*dyhl;
}
if (j<ny) {
dyhl = 2.0/(y_gpu[__I1_Y(j+2)]-y_gpu[__I1_Y(j)]);
fl2o = sigx*dyhl;
fl3o = sigy*dyhl;
fl4o = sigz*dyhl;
fl5o = sigq*dyhl;
}
}

}
}


extern "C"{
void visflx_y_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *y_gpu,real *cp_coeff_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((visflx_y_kernel),grid,block,0,stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,y_gpu,cp_coeff_gpu);
}
}



__global__ void  visflx_z_kernel(int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *z_gpu,real *cp_coeff_gpu){
//Kernel for visflx_z_kernel
real uu;real vv;real ww;
real tt;real mu;real qq;
real uup;real vvp;real wwp;
real ttp;real mup;real qqp;
real sigq;real sigx;real sigy;
real sigz;real sigq_tt;real sigq_qq;
real dzhl;real fl2o;real fl3o;
real fl4o;real fl5o;real ttf;
real muf;real cploc;
int i;int j;int k;
int iv;int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int k=0; k<nz+1; k++){
uu = w_aux_gpu[__I4_W_AUX(i,j,k ,2)];
uup = w_aux_gpu[__I4_W_AUX(i,j,k+1,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k ,3)];
vvp = w_aux_gpu[__I4_W_AUX(i,j,k+1,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k ,4)];
wwp = w_aux_gpu[__I4_W_AUX(i,j,k+1,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(i,j,k+1,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k ,7)];
mup = w_aux_gpu[__I4_W_AUX(i,j,k+1,7)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
qqp = 0.50*(uup*uup+vvp*vvp+wwp*wwp);
if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
ttf = 0.50*(tt+ttp);
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((ttf/t0),ll);
}
}
sigx = uup-uu;
sigy = vvp-vv;
sigz = wwp-ww;
sigq_tt = ttp-tt;
sigq_qq = qqp-qq;
muf = mu+mup;
muf = 0.50*muf/(z_gpu[__I1_Z(k+1)]-z_gpu[__I1_Z(k)]);
sigx = sigx*muf;
sigy = sigy*muf;
sigz = sigz*muf;
sigq = (sigq_tt*cploc/prandtl+sigq_qq)*muf;
if (k>0) {
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] + fl2o-sigx*dzhl;
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] + fl3o-sigy*dzhl;
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] + fl4o-sigz*dzhl;
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] + fl5o-sigq*dzhl;
}
if (k<nz) {
dzhl = 2.0/(z_gpu[__I1_Z(k+2)]-z_gpu[__I1_Z(k)]);
fl2o = sigx*dzhl;
fl3o = sigy*dzhl;
fl4o = sigz*dzhl;
fl5o = sigq*dzhl;
}
}

}
}


extern "C"{
void visflx_z_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int calorically_perfect,int indx_cp_l,int indx_cp_r,real prandtl,real t0,real *fl_gpu,real *w_aux_gpu,real *z_gpu,real *cp_coeff_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((visflx_z_kernel),grid,block,0,stream,nx,ny,nz,nv,ng,calorically_perfect,indx_cp_l,indx_cp_r,prandtl,t0,fl_gpu,w_aux_gpu,z_gpu,cp_coeff_gpu);
}
}



__global__ void  recyc_exchange_kernel_1(int irecyc,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf1s_gpu){
//Kernel for recyc_exchange_kernel_1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
wbuf1s_gpu[__I4_WBUF1S(i,j,k,m)] = w_gpu[__I4_W(irecyc+1-i,j,k,m)];
}

}
}


extern "C"{
void recyc_exchange_kernel_1_wrapper(hipStream_t stream,int irecyc,int nx,int ny,int nz,int ng,int nv,real *w_gpu,real *wbuf1s_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((recyc_exchange_kernel_1),grid,block,0,stream,irecyc,nx,ny,nz,ng,nv,w_gpu,wbuf1s_gpu);
}
}



__global__ void  recyc_exchange_kernel_2(int n1_start_recv,int n1_start_send,int n1_end_recv,int nx,int ny,int nz,int ng,int nv,real *wbuf1r_gpu,real *wrecyc_gpu){
//Kernel for recyc_exchange_kernel_2
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,n1_start_recv);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,n1_end_recv,1)){
for(int m=1; m<nv+1; m++){
wrecyc_gpu[__I4_WRECYC(i,j,k,m)] = wbuf1r_gpu[__I4_WBUF1R(i,j,k-n1_start_recv+n1_start_send,m)];
}

}
}


extern "C"{
void recyc_exchange_kernel_2_wrapper(hipStream_t stream,int n1_start_recv,int n1_start_send,int n1_end_recv,int nx,int ny,int nz,int ng,int nv,real *wbuf1r_gpu,real *wrecyc_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((n1_end_recv)-(n1_start_recv)+1,block.z));

hipLaunchKernelGGL((recyc_exchange_kernel_2),grid,block,0,stream,n1_start_recv,n1_start_send,n1_end_recv,nx,ny,nz,ng,nv,wbuf1r_gpu,wrecyc_gpu);
}
}



__global__ void  recyc_exchange_kernel_3(int n2_start_recv,int n2_start_send,int n2_end_recv,int nx,int ny,int nz,int ng,int nv,real *wbuf2r_gpu,real *wrecyc_gpu){
//Kernel for recyc_exchange_kernel_3
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,n2_start_recv);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,n2_end_recv,1)){
for(int m=1; m<nv+1; m++){
wrecyc_gpu[__I4_WRECYC(i,j,k,m)] = wbuf2r_gpu[__I4_WBUF2R(i,j,k-n2_start_recv+n2_start_send,m)];
}

}
}


extern "C"{
void recyc_exchange_kernel_3_wrapper(hipStream_t stream,int n2_start_recv,int n2_start_send,int n2_end_recv,int nx,int ny,int nz,int ng,int nv,real *wbuf2r_gpu,real *wrecyc_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((n2_end_recv)-(n2_start_recv)+1,block.z));

hipLaunchKernelGGL((recyc_exchange_kernel_3),grid,block,0,stream,n2_start_recv,n2_start_send,n2_end_recv,nx,ny,nz,ng,nv,wbuf2r_gpu,wrecyc_gpu);
}
}

__device__ real get_e_from_temperature_dev_bcextr_sub_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcextr_sub_kernel1_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcextr_sub_kernel1(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
//Kernel for bcextr_sub_kernel1
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int i;int j;int k;
int l;int m;int iercuda;

l = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(1,j,k,1)];
rhou = w_gpu[__I4_W(1,j,k,2)];
rhov = w_gpu[__I4_W(1,j,k,3)];
rhow = w_gpu[__I4_W(1,j,k,4)];
tt = p0/rho/rgas0;
ee = get_e_from_temperature_dev_bcextr_sub_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(1-l,j,k,1)] = rho;
w_gpu[__I4_W(1-l,j,k,2)] = rhou;
w_gpu[__I4_W(1-l,j,k,3)] = rhov;
w_gpu[__I4_W(1-l,j,k,4)] = rhow;
w_gpu[__I4_W(1-l,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bcextr_sub_kernel1_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcextr_sub_kernel1),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu);
}
}
__device__ real get_e_from_temperature_dev_bcextr_sub_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcextr_sub_kernel2_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcextr_sub_kernel2(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
//Kernel for bcextr_sub_kernel2
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int i;int j;int k;
int l;int m;int iercuda;

l = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(nx,j,k,1)];
rhou = w_gpu[__I4_W(nx,j,k,2)];
rhov = w_gpu[__I4_W(nx,j,k,3)];
rhow = w_gpu[__I4_W(nx,j,k,4)];
tt = p0/rho/rgas0;
ee = get_e_from_temperature_dev_bcextr_sub_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(nx+l,j,k,1)] = rho;
w_gpu[__I4_W(nx+l,j,k,2)] = rhou;
w_gpu[__I4_W(nx+l,j,k,3)] = rhov;
w_gpu[__I4_W(nx+l,j,k,4)] = rhow;
w_gpu[__I4_W(nx+l,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bcextr_sub_kernel2_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcextr_sub_kernel2),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu);
}
}
__device__ real get_e_from_temperature_dev_bcextr_sub_kernel3_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcextr_sub_kernel3_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcextr_sub_kernel3(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
//Kernel for bcextr_sub_kernel3
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int i;int j;int k;
int l;int m;int iercuda;

l = __GIDX(x,1);
i = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(i,ny,k,1)];
rhou = w_gpu[__I4_W(i,ny,k,2)];
rhov = w_gpu[__I4_W(i,ny,k,3)];
rhow = w_gpu[__I4_W(i,ny,k,4)];
tt = p0/rho/rgas0;
ee = get_e_from_temperature_dev_bcextr_sub_kernel3_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(i,ny+l,k,1)] = rho;
w_gpu[__I4_W(i,ny+l,k,2)] = rhou;
w_gpu[__I4_W(i,ny+l,k,3)] = rhov;
w_gpu[__I4_W(i,ny+l,k,4)] = rhow;
w_gpu[__I4_W(i,ny+l,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bcextr_sub_kernel3_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcextr_sub_kernel3),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu);
}
}
__device__ real get_e_from_temperature_dev_bcextr_sub_kernel4_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcextr_sub_kernel4_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcextr_sub_kernel4(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
//Kernel for bcextr_sub_kernel4
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
l = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(l,ng,1)){
rho = w_gpu[__I4_W(i,j,1,1)];
rhou = w_gpu[__I4_W(i,j,1,2)];
rhov = w_gpu[__I4_W(i,j,1,3)];
rhow = w_gpu[__I4_W(i,j,1,4)];
tt = p0/rho/rgas0;
ee = get_e_from_temperature_dev_bcextr_sub_kernel4_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(i,j,1-l,1)] = rho;
w_gpu[__I4_W(i,j,1-l,2)] = rhou;
w_gpu[__I4_W(i,j,1-l,3)] = rhov;
w_gpu[__I4_W(i,j,1-l,4)] = rhow;
w_gpu[__I4_W(i,j,1-l,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bcextr_sub_kernel4_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcextr_sub_kernel4),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu);
}
}
__device__ real get_e_from_temperature_dev_bcextr_sub_kernel5_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcextr_sub_kernel5_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcextr_sub_kernel5(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
//Kernel for bcextr_sub_kernel5
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
l = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(l,ng,1)){
rho = w_gpu[__I4_W(i,j,nz,1)];
rhou = w_gpu[__I4_W(i,j,nz,2)];
rhov = w_gpu[__I4_W(i,j,nz,3)];
rhow = w_gpu[__I4_W(i,j,nz,4)];
tt = p0/rho/rgas0;
ee = get_e_from_temperature_dev_bcextr_sub_kernel5_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(i,j,nz+l,1)] = rho;
w_gpu[__I4_W(i,j,nz+l,2)] = rhou;
w_gpu[__I4_W(i,j,nz+l,3)] = rhov;
w_gpu[__I4_W(i,j,nz+l,4)] = rhow;
w_gpu[__I4_W(i,j,nz+l,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bcextr_sub_kernel5_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real p0,real t0,real rgas0,real *cv_coeff_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((ng)-(1)+1,block.z));

hipLaunchKernelGGL((bcextr_sub_kernel5),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,p0,t0,rgas0,cv_coeff_gpu,w_gpu);
}
}



__global__ void  bcrecyc_kernel_1(int nx,int ny,int nz,int ng,int nv,real *wrecycav_gpu,real *wrecyc_gpu){
//Kernel for bcrecyc_kernel_1
int i;int j;int k;
int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
m = __GIDX(z,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)&&loop_cond(m,nv,1)){
wrecycav_gpu[__I3_WRECYCAV(i,j,m)] = 0.0;
for(int k=1; k<nz+1; k++){
wrecycav_gpu[__I3_WRECYCAV(i,j,m)] = wrecycav_gpu[__I3_WRECYCAV(i,j,m)]+wrecyc_gpu[__I4_WRECYC(i,j,k,m)];
}

}
}


extern "C"{
void bcrecyc_kernel_1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real *wrecycav_gpu,real *wrecyc_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nv)-(1)+1,block.z));

hipLaunchKernelGGL((bcrecyc_kernel_1),grid,block,0,stream,nx,ny,nz,ng,nv,wrecycav_gpu,wrecyc_gpu);
}
}



__global__ void  bcrecyc_kernel_2(int nx,int ny,int nz,int nzmax,int ng,real *wrecycav_gpu,real *wrecyc_gpu){
//Kernel for bcrecyc_kernel_2
real ufav;real vfav;real wfav;
real rhom;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,ng,1)&&loop_cond(j,ny,1)){
ufav = wrecycav_gpu[__I3_WRECYCAV(i,j,2)]/wrecycav_gpu[__I3_WRECYCAV(i,j,1)];
vfav = wrecycav_gpu[__I3_WRECYCAV(i,j,3)]/wrecycav_gpu[__I3_WRECYCAV(i,j,1)];
wfav = wrecycav_gpu[__I3_WRECYCAV(i,j,4)]/wrecycav_gpu[__I3_WRECYCAV(i,j,1)];
rhom = wrecycav_gpu[__I3_WRECYCAV(i,j,1)]/nzmax;
for(int k=1; k<nz+1; k++){
wrecyc_gpu[__I4_WRECYC(i,j,k,2)] = wrecyc_gpu[__I4_WRECYC(i,j,k,2)]/wrecyc_gpu[__I4_WRECYC(i,j,k,1)]-ufav;
wrecyc_gpu[__I4_WRECYC(i,j,k,3)] = wrecyc_gpu[__I4_WRECYC(i,j,k,3)]/wrecyc_gpu[__I4_WRECYC(i,j,k,1)]-vfav;
wrecyc_gpu[__I4_WRECYC(i,j,k,4)] = wrecyc_gpu[__I4_WRECYC(i,j,k,4)]/wrecyc_gpu[__I4_WRECYC(i,j,k,1)]-wfav;
wrecyc_gpu[__I4_WRECYC(i,j,k,1)] = wrecyc_gpu[__I4_WRECYC(i,j,k,1)]-rhom;
}

}
}


extern "C"{
void bcrecyc_kernel_2_wrapper(hipStream_t stream,int nx,int ny,int nz,int nzmax,int ng,real *wrecycav_gpu,real *wrecyc_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcrecyc_kernel_2),grid,block,0,stream,nx,ny,nz,nzmax,ng,wrecycav_gpu,wrecyc_gpu);
}
}

__device__ real get_e_from_temperature_dev_bcrecyc_kernel_3_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcrecyc_kernel_3_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcrecyc_kernel_3(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int rand_type,real p0,real rgas0,real betarecyc,real glund1,real t0,real u0,real u0_02,int *map_j_inn_gpu,int *map_j_out_gpu,int *map_j_out_blend_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu,real *wmean_gpu,real *wrecyc_gpu,real *inflow_random_plane_gpu,real *w_gpu,real *weta_inflow_gpu,real *yplus_inflow_gpu,real *eta_inflow_gpu,real *yplus_recyc_gpu,real *eta_recyc_gpu,real *eta_recyc_blend_gpu){
//Kernel for bcrecyc_kernel_3
real eta;real weta;real weta1;
real bdamp;real disty_inn;real disty_out;
real rhofluc;real ufluc;real vfluc;
real wfluc;real rhof_inn;real rhof_out;
real uf_inn;real uf_out;real vf_inn;
real vf_out;real wf_inn;real wf_out;
real etamin;real rhomean;real uumean;
real vvmean;real wwmean;real tmean;
real rho;real uu;real vv;
real ww;real rhou;real rhov;
real rhow;real tt;real ee;
real tfluc;
int i;int j;int k;
int iercuda;int j_inn;int j_out;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
eta = eta_inflow_gpu[__I1_ETA_INFLOW(j)];
etamin = min(eta,1.0);
weta = weta_inflow_gpu[__I1_WETA_INFLOW(j)];
weta1 = 1.0-weta;
bdamp = 0.50*(1.0-tanh(4.0*(eta_inflow_gpu[__I1_ETA_INFLOW(j)]-2.0)));
j_inn = map_j_inn_gpu[__I1_MAP_J_INN(j)];
j_out = map_j_out_gpu[__I1_MAP_J_OUT(j)];
disty_inn = (yplus_inflow_gpu[__I1_YPLUS_INFLOW(j)]-yplus_recyc_gpu[__I1_YPLUS_RECYC(j_inn)])/(yplus_recyc_gpu[__I1_YPLUS_RECYC(j_inn+1)]-yplus_recyc_gpu[__I1_YPLUS_RECYC(j_inn)]);
disty_out = (eta_inflow_gpu[__I1_ETA_INFLOW(j)]-eta_recyc_gpu[__I1_ETA_RECYC(j_out)])/(eta_recyc_gpu[__I1_ETA_RECYC(j_out+1)]-eta_recyc_gpu[__I1_ETA_RECYC(j_out)]);
for(int i=1; i<ng+1; i++){
rhomean = wmean_gpu[__I3_WMEAN(1-i,j,1)];
uumean = wmean_gpu[__I3_WMEAN(1-i,j,2)]/rhomean;
vvmean = wmean_gpu[__I3_WMEAN(1-i,j,3)]/rhomean;
wwmean = wmean_gpu[__I3_WMEAN(1-i,j,4)]/rhomean;
tmean = p0/rhomean/rgas0;
if (j==1||j_inn>=ny||j_out>=ny) {
rhofluc = 0.0;
ufluc = 0.0;
vfluc = 0.0;
wfluc = 0.0;
}else {
rhof_inn = wrecyc_gpu[__I4_WRECYC(i,j_inn,k,1)]*(1.0-disty_inn)+wrecyc_gpu[__I4_WRECYC(i,j_inn+1,k,1)]*disty_inn;
rhof_out = wrecyc_gpu[__I4_WRECYC(i,j_out,k,1)]*(1.0-disty_out)+wrecyc_gpu[__I4_WRECYC(i,j_out+1,k,1)]*disty_out;
uf_inn = wrecyc_gpu[__I4_WRECYC(i,j_inn,k,2)]*(1.0-disty_inn)+wrecyc_gpu[__I4_WRECYC(i,j_inn+1,k,2)]*disty_inn;
uf_out = wrecyc_gpu[__I4_WRECYC(i,j_out,k,2)]*(1.0-disty_out)+wrecyc_gpu[__I4_WRECYC(i,j_out+1,k,2)]*disty_out;
vf_inn = wrecyc_gpu[__I4_WRECYC(i,j_inn,k,3)]*(1.0-disty_inn)+wrecyc_gpu[__I4_WRECYC(i,j_inn+1,k,3)]*disty_inn;
vf_out = wrecyc_gpu[__I4_WRECYC(i,j_out,k,3)]*(1.0-disty_out)+wrecyc_gpu[__I4_WRECYC(i,j_out+1,k,3)]*disty_out;
wf_inn = wrecyc_gpu[__I4_WRECYC(i,j_inn,k,4)]*(1.0-disty_inn)+wrecyc_gpu[__I4_WRECYC(i,j_inn+1,k,4)]*disty_inn;
wf_out = wrecyc_gpu[__I4_WRECYC(i,j_out,k,4)]*(1.0-disty_out)+wrecyc_gpu[__I4_WRECYC(i,j_out+1,k,4)]*disty_out;
rhofluc = rhof_inn*weta1+rhof_out*weta;
ufluc = uf_inn*weta1+ uf_out*weta;
vfluc = vf_inn*weta1+ vf_out*weta;
wfluc = wf_inn*weta1+ wf_out*weta;
rhofluc = rhofluc*bdamp;
ufluc = ufluc *bdamp*betarecyc;
vfluc = vfluc *bdamp*betarecyc;
wfluc = wfluc *bdamp*betarecyc;
ufluc = ufluc+u0_02*(inflow_random_plane_gpu[__I3_INFLOW_RANDOM_PLANE(j,k,1)]-0.50)*etamin;
vfluc = vfluc+u0_02*(inflow_random_plane_gpu[__I3_INFLOW_RANDOM_PLANE(j,k,2)]-0.50)*etamin;
wfluc = wfluc+u0_02*(inflow_random_plane_gpu[__I3_INFLOW_RANDOM_PLANE(j,k,3)]-0.50)*etamin;
}
rhofluc = max(-0.750*rhomean,rhofluc);
rhofluc = min( 4.00*rhomean,rhofluc);
rho = rhomean + rhofluc;
tt = p0/rho/rgas0;
uu = uumean + ufluc;
vv = vvmean + vfluc;
ww = wwmean + wfluc;
rhou = rho*uu;
rhov = rho*vv;
rhow = rho*ww;
ee = get_e_from_temperature_dev_bcrecyc_kernel_3_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(1-i,j,k,1)] = rho;
w_gpu[__I4_W(1-i,j,k,2)] = rhou;
w_gpu[__I4_W(1-i,j,k,3)] = rhov;
w_gpu[__I4_W(1-i,j,k,4)] = rhow;
w_gpu[__I4_W(1-i,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;
}

}
}


extern "C"{
void bcrecyc_kernel_3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int rand_type,real p0,real rgas0,real betarecyc,real glund1,real t0,real u0,real u0_02,int *map_j_inn_gpu,int *map_j_out_gpu,int *map_j_out_blend_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu,real *wmean_gpu,real *wrecyc_gpu,real *inflow_random_plane_gpu,real *w_gpu,real *weta_inflow_gpu,real *yplus_inflow_gpu,real *eta_inflow_gpu,real *yplus_recyc_gpu,real *eta_recyc_gpu,real *eta_recyc_blend_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcrecyc_kernel_3),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rand_type,p0,rgas0,betarecyc,glund1,t0,u0,u0_02,map_j_inn_gpu,map_j_out_gpu,map_j_out_blend_gpu,cv_coeff_gpu,cp_coeff_gpu,wmean_gpu,wrecyc_gpu,inflow_random_plane_gpu,w_gpu,weta_inflow_gpu,yplus_inflow_gpu,eta_inflow_gpu,yplus_recyc_gpu,eta_recyc_gpu,eta_recyc_blend_gpu);
}
}

__device__ real get_e_from_temperature_dev_bclam_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bclam_kernel_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bclam_kernel(int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real p0,real rgas0,real t0,real *cv_coeff_gpu,real *w_gpu,real *wmean_gpu){
//Kernel for bclam_kernel
real rho;real rhou;real rhov;
real rhow;real tt;real ee;
int j;int k;int l;
int iercuda;

l = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
rho = wmean_gpu[__I3_WMEAN(1-l,j,1)];
rhou = wmean_gpu[__I3_WMEAN(1-l,j,2)];
rhov = wmean_gpu[__I3_WMEAN(1-l,j,3)];
rhow = wmean_gpu[__I3_WMEAN(1-l,j,4)];
tt = p0/rho/rgas0;
w_gpu[__I4_W(1-l,j,k,1)] = rho;
w_gpu[__I4_W(1-l,j,k,2)] = rhou;
w_gpu[__I4_W(1-l,j,k,3)] = rhov;
w_gpu[__I4_W(1-l,j,k,4)] = rhow;
ee = get_e_from_temperature_dev_bclam_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
w_gpu[__I4_W(1-l,j,k,5)]=rho*ee+0.50*(((rhou)*(rhou))+((rhov)*(rhov))+((rhow)*(rhow)))/rho;

}
}


extern "C"{
void bclam_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real p0,real rgas0,real t0,real *cv_coeff_gpu,real *w_gpu,real *wmean_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bclam_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,ilat,p0,rgas0,t0,cv_coeff_gpu,w_gpu,wmean_gpu);
}
}



__global__ void  bcfree_kernel1(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel1
int i;int j;int k;
int l;int m;int iercuda;

l = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(1-l,j,k,m)] = winf_gpu[__I1_WINF(m)];
}

}
}


extern "C"{
void bcfree_kernel1_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcfree_kernel1),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}


__global__ void  bcfree_kernel2(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel2
int i;int j;int k;
int l;int m;int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+l,j,k,m)] = winf_gpu[__I1_WINF(m)];
}
}

}
}


extern "C"{
void bcfree_kernel2_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_kernel2),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}


__global__ void  bcfree_kernel3(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel3
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,1-l,k,m)] = winf_gpu[__I1_WINF(m)];
}
}

}
}


extern "C"{
void bcfree_kernel3_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_kernel3),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}


__global__ void  bcfree_kernel4(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel4
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+l,k,m)] = winf_gpu[__I1_WINF(m)];
}
}

}
}


extern "C"{
void bcfree_kernel4_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_kernel4),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}


__global__ void  bcfree_kernel5(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel5
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,1-l,m)] = winf_gpu[__I1_WINF(m)];
}
}

}
}


extern "C"{
void bcfree_kernel5_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_kernel5),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}


__global__ void  bcfree_kernel6(int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
//Kernel for bcfree_kernel6
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,nz+l,m)] = winf_gpu[__I1_WINF(m)];
}
}

}
}


extern "C"{
void bcfree_kernel6_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int nv,real *winf_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_kernel6),grid,block,0,stream,ilat,nx,ny,nz,ng,nv,winf_gpu,w_gpu);
}
}

__device__ real get_gamloc_dev_bcfree_sub_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bcfree_sub_kernel1_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}

__device__ real get_temperature_from_e_dev_bcfree_sub_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcfree_sub_kernel1_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcfree_sub_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcfree_sub_kernel1_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcfree_sub_kernel1(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real ptot0,real ttot0,real rgas0,real aoa,real t0,real tol_iter_nr,real cosangle,real sinangle,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu){
//Kernel for bcfree_sub_kernel1
real rho;real rhou;real rhov;
real rhow;real rhoe;real ri;
real uu;real vv;real ww;
real qq;real ee;real tt;
real pp;real del;real gamloc;
real rmaf;real rml;real vel_mod;
real rmfac;
int i;int j;int k;
int l;int m;int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(1,j,k,1)];
rhou = w_gpu[__I4_W(1,j,k,2)];
rhov = w_gpu[__I4_W(1,j,k,3)];
rhow = w_gpu[__I4_W(1,j,k,4)];
rhoe = w_gpu[__I4_W(1,j,k,5)];
ri = 1.0/rho;
uu = rhou*ri;
vv = rhov*ri;
ww = rhow*ri;
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcfree_sub_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(1,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
gamloc = get_gamloc_dev_bcfree_sub_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
del = 0.50*(gamloc-1.0);
rmfac=pow((pp/ptot0),(-(gamloc-1.0)/gamloc));
rml = sqrt((rmfac-1.0)/del);
tt = ttot0/rmfac;
vel_mod = rml*sqrt(gamloc*rgas0*tt);
uu = vel_mod*cosangle;
vv = vel_mod*sinangle;
ww = 0.0;
rho = pp/tt/rgas0;
rhou = rho*uu;
rhov = rho*vv;
rhow = rho*ww;
ee = get_e_from_temperature_dev_bcfree_sub_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(1-l,j,k,1)] = rho;
w_gpu[__I4_W(1-l,j,k,2)] = rhou;
w_gpu[__I4_W(1-l,j,k,3)] = rhov;
w_gpu[__I4_W(1-l,j,k,4)] = rhow;
w_gpu[__I4_W(1-l,j,k,5)] = rho*ee+0.50*(rhou*rhou+rhov*rhov+rhow*rhow)/rho;
}

}
}


extern "C"{
void bcfree_sub_kernel1_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real ptot0,real ttot0,real rgas0,real aoa,real t0,real tol_iter_nr,real cosangle,real sinangle,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_sub_kernel1),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu);
}
}
__device__ real get_gamloc_dev_bcfree_sub_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bcfree_sub_kernel2_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}

__device__ real get_temperature_from_e_dev_bcfree_sub_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcfree_sub_kernel2_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcfree_sub_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcfree_sub_kernel2_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcfree_sub_kernel2(int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real ptot0,real ttot0,real rgas0,real aoa,real t0,real tol_iter_nr,real cosangle,real sinangle,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu){
//Kernel for bcfree_sub_kernel2
real rho;real rhou;real rhov;
real rhow;real rhoe;real ri;
real uu;real vv;real ww;
real qq;real ee;real tt;
real pp;real del;real gamloc;
real rmaf;real rml;real vel_mod;
real rmfac;
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(i,ny,k,1)];
rhou = w_gpu[__I4_W(i,ny,k,2)];
rhov = w_gpu[__I4_W(i,ny,k,3)];
rhow = w_gpu[__I4_W(i,ny,k,4)];
rhoe = w_gpu[__I4_W(i,ny,k,5)];
ri = 1.0/rho;
uu = rhou*ri;
vv = rhov*ri;
ww = rhow*ri;
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcfree_sub_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,ny,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
gamloc = get_gamloc_dev_bcfree_sub_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
del = 0.50*(gamloc-1.0);
rmfac=pow((pp/ptot0),(-(gamloc-1.0)/gamloc));
rml = sqrt((rmfac-1.0)/del);
tt = ttot0/rmfac;
vel_mod = rml*sqrt(gamloc*rgas0*tt);
uu = vel_mod*cosangle;
vv = vel_mod*sinangle;
ww = 0.0;
rho = pp/tt/rgas0;
rhou = rho*uu;
rhov = rho*vv;
rhow = rho*ww;
ee = get_e_from_temperature_dev_bcfree_sub_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,ny+l,k,1)] = rho;
w_gpu[__I4_W(i,ny+l,k,2)] = rhou;
w_gpu[__I4_W(i,ny+l,k,3)] = rhov;
w_gpu[__I4_W(i,ny+l,k,4)] = rhow;
w_gpu[__I4_W(i,ny+l,k,5)] = rho*ee+0.50*(rhou*rhou+rhov*rhov+rhow*rhow)/rho;
}

}
}


extern "C"{
void bcfree_sub_kernel2_wrapper(hipStream_t stream,int ilat,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real ptot0,real ttot0,real rgas0,real aoa,real t0,real tol_iter_nr,real cosangle,real sinangle,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *cp_coeff_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcfree_sub_kernel2),grid,block,0,stream,ilat,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ptot0,ttot0,rgas0,aoa,t0,tol_iter_nr,cosangle,sinangle,w_gpu,w_aux_gpu,cv_coeff_gpu,cp_coeff_gpu);
}
}



__global__ void  bcshock_kernel(int nx,int ny,int nz,int ng,int nv,int ilat,real xshock_imp,real shock_angle,real tanhfacs,real tanhlen,real *winf_gpu,real *winf_past_shock_gpu,real *w_gpu,real *x_gpu,real *y_gpu){
//Kernel for bcshock_kernel
real xsh;real xx;real tanhf;
real dwinf;
int i;int k;int l;
int m;int iercuda;

l = __GIDX(x,0);
i = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(l,ng,1)&&loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
xsh = xshock_imp-y_gpu[__I1_Y(ny+l)]/tan(shock_angle);
xx = x_gpu[__I1_X(i)]-xsh;
if (abs(xx)>tanhlen) {
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+l,k,m)] = w_gpu[__I4_W(i,ny,k,m)];
}
}else {
tanhf = 0.50*(1.0+tanh(xx/tanhfacs));
for(int m=1; m<nv+1; m++){
dwinf = winf_past_shock_gpu[__I1_WINF_PAST_SHOCK(m)]-winf_gpu[__I1_WINF(m)];
w_gpu[__I4_W(i,ny+l,k,m)] = winf_gpu[__I1_WINF(m)]+dwinf*tanhf;
}
}

}
}


extern "C"{
void bcshock_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real xshock_imp,real shock_angle,real tanhfacs,real tanhlen,real *winf_gpu,real *winf_past_shock_gpu,real *w_gpu,real *x_gpu,real *y_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((ng)-(0)+1,block.x),divideAndRoundUp((nx)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcshock_kernel),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,xshock_imp,shock_angle,tanhfacs,tanhlen,winf_gpu,winf_past_shock_gpu,w_gpu,x_gpu,y_gpu);
}
}



__global__ void  bcextr_var_kernel1(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel1
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(1-l,j,k,1)] = w_gpu[__I4_W(1,j,k,1)];
}

}
}


extern "C"{
void bcextr_var_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel1),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}


__global__ void  bcextr_var_kernel2(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel2
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(nx+l,j,k,1)] = w_gpu[__I4_W(nx,j,k,1)];
}

}
}


extern "C"{
void bcextr_var_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel2),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}


__global__ void  bcextr_var_kernel3(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel3
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,1-l,k,1)] = w_gpu[__I4_W(i,1,k,1)];
}

}
}


extern "C"{
void bcextr_var_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel3),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}


__global__ void  bcextr_var_kernel4(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel4
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,ny+l,k,1)] = w_gpu[__I4_W(i,ny,k,1)];
}

}
}


extern "C"{
void bcextr_var_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel4),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}


__global__ void  bcextr_var_kernel5(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel5
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,j,1-l,1)] = w_gpu[__I4_W(i,j,1,1)];
}

}
}


extern "C"{
void bcextr_var_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel5),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}


__global__ void  bcextr_var_kernel6(int nx,int ny,int nz,int ng,real *w_gpu){
//Kernel for bcextr_var_kernel6
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,j,nz+l,1)] = w_gpu[__I4_W(i,j,nz,1)];
}

}
}


extern "C"{
void bcextr_var_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_var_kernel6),grid,block,0,stream,nx,ny,nz,ng,w_gpu);
}
}



__global__ void  bcextr_airfoil_var_kernel1(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel1
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(1-l,j,k,1)] = 2.0*w_gpu[__I4_W(2-l,j,k,1)]-w_gpu[__I4_W(3-l,j,k,1)];
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel1),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}


__global__ void  bcextr_airfoil_var_kernel2(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel2
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(nx+l,j,k,1)] = 2.0*w_gpu[__I4_W(nx+l-1,j,k,1)]-w_gpu[__I4_W(nx+l-2,j,k,1)];
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel2),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}


__global__ void  bcextr_airfoil_var_kernel3(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel3
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
w_gpu[__I4_W(i,1-l,k,1)] = 2.0*w_gpu[__I4_W(i,2-l,k,1)]-w_gpu[__I4_W(i,3-l,k,1)];
}
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel3),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}


__global__ void  bcextr_airfoil_var_kernel4(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel4
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,ny+l,k,1)] = 2.0*w_gpu[__I4_W(i,ny+l-1,k,1)]-w_gpu[__I4_W(i,ny+l-2,k,1)];
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel4),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}


__global__ void  bcextr_airfoil_var_kernel5(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel5
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,j,1-l,1)] = 2.0*w_gpu[__I4_W(i,j,2-l,1)]-w_gpu[__I4_W(i,j,3-l,1)];
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel5),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}


__global__ void  bcextr_airfoil_var_kernel6(int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
//Kernel for bcextr_airfoil_var_kernel6
int ilat;int i;int j;
int k;int l;int m;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,j,nz+l,1)] = 2.0*w_gpu[__I4_W(i,j,nz+l-1,1)]-w_gpu[__I4_W(i,j,nz+l-2,1)];
}

}
}


extern "C"{
void bcextr_airfoil_var_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndim,int ileftx,int irightx,int ileftz,int irightz,int *wall_tag_gpu,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_airfoil_var_kernel6),grid,block,0,stream,nx,ny,nz,ng,ndim,ileftx,irightx,ileftz,irightz,wall_tag_gpu,w_gpu);
}
}



__global__ void  bcextr_kernel1(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel1
int i;int j;int k;
int l;int m;int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(1-l,j,k,m)] = w_gpu[__I4_W(1,j,k,m)];
}
}

}
}


extern "C"{
void bcextr_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel1),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}


__global__ void  bcextr_kernel2(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel2
int i;int j;int k;
int l;int m;int iercuda;

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(nx+l,j,k,m)] = w_gpu[__I4_W(nx,j,k,m)];
}
}

}
}


extern "C"{
void bcextr_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel2),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}


__global__ void  bcextr_kernel3(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel3
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,1-l,k,m)] = w_gpu[__I4_W(i,1,k,m)];
}
}

}
}


extern "C"{
void bcextr_kernel3_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel3),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}


__global__ void  bcextr_kernel4(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel4
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,ny+l,k,m)] = w_gpu[__I4_W(i,ny,k,m)];
}
}

}
}


extern "C"{
void bcextr_kernel4_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel4),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}


__global__ void  bcextr_kernel5(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel5
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,1-l,m)] = w_gpu[__I4_W(i,j,1,m)];
}
}

}
}


extern "C"{
void bcextr_kernel5_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel5),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}


__global__ void  bcextr_kernel6(int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
//Kernel for bcextr_kernel6
int i;int j;int k;
int l;int m;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int l=1; l<ng+1; l++){
for(int m=1; m<nv+1; m++){
w_gpu[__I4_W(i,j,nz+l,m)] = w_gpu[__I4_W(i,j,nz,m)];
}
}

}
}


extern "C"{
void bcextr_kernel6_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bcextr_kernel6),grid,block,0,stream,nx,ny,nz,ng,nv,ilat,w_gpu);
}
}



__global__ void  bcsym_kernel(int nx,int ny,int nz,int ng,int ilat,real *w_gpu){
//Kernel for bcsym_kernel
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
w_gpu[__I4_W(i,1-l,k,1)] = w_gpu[__I4_W(i,1+l,k,1)];
w_gpu[__I4_W(i,1-l,k,2)] = w_gpu[__I4_W(i,1+l,k,2)];
w_gpu[__I4_W(i,1-l,k,3)] = w_gpu[__I4_W(i,1+l,k,3)];
w_gpu[__I4_W(i,1-l,k,4)] = w_gpu[__I4_W(i,1+l,k,4)];
w_gpu[__I4_W(i,1-l,k,5)] = w_gpu[__I4_W(i,1+l,k,5)];
}

}
}


extern "C"{
void bcsym_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ilat,real *w_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcsym_kernel),grid,block,0,stream,nx,ny,nz,ng,ilat,w_gpu);
}
}



__global__ void  bcsym_c2_kernel1(int nx,int ny,int nz,int ng,int ilat,real *w_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
//Kernel for bcsym_c2_kernel1
real abcsym_i31;real abcsym_i32;real abcsym_i41;
real abcsym_i42;real tauvers_i31;real tauvers_i32;
real tauvers_i41;real tauvers_i42;real taumod;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
taumod=sqrt(((dxdcsic2_gpu[__I2_DXDCSIC2(i,1)])*(dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]))+((dydcsic2_gpu[__I2_DYDCSIC2(i,1)])*(dydcsic2_gpu[__I2_DYDCSIC2(i,1)])));
tauvers_i31 = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/taumod;
tauvers_i32 = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/taumod;
abcsym_i31=((tauvers_i31)*(tauvers_i31))-((tauvers_i32)*(tauvers_i32));
abcsym_i32 = 2.0*tauvers_i31*tauvers_i32;
w_gpu[__I4_W(i,1-l,k,1)] = w_gpu[__I4_W(i,1+l,k,1)];
w_gpu[__I4_W(i,1-l,k,4)] = w_gpu[__I4_W(i,1+l,k,4)];
w_gpu[__I4_W(i,1-l,k,5)] = w_gpu[__I4_W(i,1+l,k,5)];
w_gpu[__I4_W(i,1-l,k,2)] = abcsym_i31*w_gpu[__I4_W(i,1+l,k,2)] + abcsym_i32*w_gpu[__I4_W(i,1+l,k,3)];
w_gpu[__I4_W(i,1-l,k,3)] = abcsym_i32*w_gpu[__I4_W(i,1+l,k,2)] - abcsym_i31*w_gpu[__I4_W(i,1+l,k,3)];
}

}
}


extern "C"{
void bcsym_c2_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ilat,real *w_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcsym_c2_kernel1),grid,block,0,stream,nx,ny,nz,ng,ilat,w_gpu,dxdcsic2_gpu,dydcsic2_gpu);
}
}


__global__ void  bcsym_c2_kernel2(int nx,int ny,int nz,int ng,int ilat,real *w_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
//Kernel for bcsym_c2_kernel2
real abcsym_i31;real abcsym_i32;real abcsym_i41;
real abcsym_i42;real tauvers_i31;real tauvers_i32;
real tauvers_i41;real tauvers_i42;real taumod;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ng+1; l++){
taumod=sqrt(((dxdcsic2_gpu[__I2_DXDCSIC2(i,ny)])*(dxdcsic2_gpu[__I2_DXDCSIC2(i,ny)]))+((dydcsic2_gpu[__I2_DYDCSIC2(i,ny)])*(dydcsic2_gpu[__I2_DYDCSIC2(i,ny)])));
tauvers_i41 = -dxdcsic2_gpu[__I2_DXDCSIC2(i,ny)]/taumod;
tauvers_i42 = -dydcsic2_gpu[__I2_DYDCSIC2(i,ny)]/taumod;
abcsym_i41=((tauvers_i41)*(tauvers_i41))-((tauvers_i42)*(tauvers_i42));
abcsym_i42 = 2.0*tauvers_i41*tauvers_i42;
w_gpu[__I4_W(i,ny+l,k,1)] = w_gpu[__I4_W(i,ny-l,k,1)];
w_gpu[__I4_W(i,ny+l,k,4)] = w_gpu[__I4_W(i,ny-l,k,4)];
w_gpu[__I4_W(i,ny+l,k,5)] = w_gpu[__I4_W(i,ny-l,k,5)];
w_gpu[__I4_W(i,ny+l,k,2)] = abcsym_i41*w_gpu[__I4_W(i,ny-l,k,2)] + abcsym_i42*w_gpu[__I4_W(i,ny-l,k,3)];
w_gpu[__I4_W(i,ny+l,k,3)] = abcsym_i42*w_gpu[__I4_W(i,ny-l,k,2)] - abcsym_i41*w_gpu[__I4_W(i,ny-l,k,3)];
}

}
}


extern "C"{
void bcsym_c2_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ilat,real *w_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcsym_c2_kernel2),grid,block,0,stream,nx,ny,nz,ng,ilat,w_gpu,dxdcsic2_gpu,dydcsic2_gpu);
}
}

__device__ real get_temperature_from_e_dev_bcwall_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcwall_kernel_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real tt,real t0,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_kernel_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_kernel_1(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_kernel_1
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcwall_kernel(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
//Kernel for bcwall_kernel
real rho;real uu;real vv;
real ww;real qq;real pp;
real tt;real rhoe;real ee;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
w_gpu[__I4_W(i,1,k,2)] = 0.0;
w_gpu[__I4_W(i,1,k,3)] = 0.0;
w_gpu[__I4_W(i,1,k,4)] = 0.0;
ee = get_e_from_temperature_dev_bcwall_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,twall,t0,cv_coeff_gpu);
w_gpu[__I4_W(i,1,k,5)] = w_gpu[__I4_W(i,1,k,1)]*ee;
for(int l=1; l<ng+1; l++){
rho = w_gpu[__I4_W(i,1+l,k,1)];
uu = w_gpu[__I4_W(i,1+l,k,2)]/rho;
vv = w_gpu[__I4_W(i,1+l,k,3)]/rho;
ww = w_gpu[__I4_W(i,1+l,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,1+l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcwall_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,1+l,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
tt = 2.0*twall-tt;
ee = get_e_from_temperature_dev_bcwall_kernel_1(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
rho = pp/tt/rgas0;
w_gpu[__I4_W(i,1-l,k,1)] = rho;
w_gpu[__I4_W(i,1-l,k,2)] = -rho*uu;
w_gpu[__I4_W(i,1-l,k,3)] = -rho*vv;
w_gpu[__I4_W(i,1-l,k,4)] = -rho*ww;
w_gpu[__I4_W(i,1-l,k,5)] = rho*(ee+qq);
}

}
}


extern "C"{
void bcwall_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcwall_kernel),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu);
}
}

__device__ real get_temperature_from_e_dev_bcwall_airfoil_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcwall_airfoil_kernel_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_airfoil_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real tt,real t0,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_airfoil_kernel_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_airfoil_kernel_1(int indx_cp_l,int indx_cp_r,int calorically_perfect,real tt,real t0,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_airfoil_kernel_1
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_airfoil_kernel_2(int indx_cp_l,int indx_cp_r,int calorically_perfect,real tt,real t0,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_airfoil_kernel_2
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_airfoil_kernel_3(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_airfoil_kernel_3
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcwall_airfoil_kernel(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real u0,real twall,real t0,real rgas0,real tol_iter_nr,real a_tw,real v_bs,real thic,real kx_tw,real om_tw,real time,real xtw1,real xtw2,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *xc2_gpu,real *yc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
//Kernel for bcwall_airfoil_kernel
real rho;real uu;real vv;
real ww;real qq;real pp;
real tt;real rhoe;real ee;
real dx1;real dx2;real ftanh;
real tww;real ug;real vg;
real wg;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
dx1 = (xc2_gpu[__I2_XC2(i,1)]-xtw1);
dx2 = (xc2_gpu[__I2_XC2(i,1)]-xtw2);
if(a_tw > 0.0 && wall_tag_gpu[__I1_WALL_TAG(i)] == 0) {
ftanh = 0.50*(tanh(dx1/thic)-tanh(dx2/thic));
tww = ftanh*(a_tw*sin(kx_tw*xc2_gpu[__I2_XC2(i,1)]-om_tw*time));
w_gpu[__I4_W(i,1,k,2)] = 0.0;
w_gpu[__I4_W(i,1,k,3)] = 0.0;
w_gpu[__I4_W(i,1,k,4)] = w_gpu[__I4_W(i,1,k,1)]*tww;
ee = get_e_from_temperature_dev_bcwall_airfoil_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,twall,t0,cv_coeff_gpu);
w_gpu[__I4_W(i,1,k,5)] = w_gpu[__I4_W(i,1,k,1)]*ee+w_gpu[__I4_W(i,1,k,1)]*0.50*(tww*tww);
}else if(abs(v_bs) > 0.0 && wall_tag_gpu[__I1_WALL_TAG(i)] == 0) {
ftanh = 0.50*(tanh(dx1/thic)-tanh(dx2/thic));
tww = ftanh*v_bs*u0;
w_gpu[__I4_W(i,1,k,2)] = w_gpu[__I4_W(i,1,k,1)]*tww*dxdetanc2_gpu[__I2_DXDETANC2(i,1)];
w_gpu[__I4_W(i,1,k,3)] = w_gpu[__I4_W(i,1,k,1)]*tww*dydetanc2_gpu[__I2_DYDETANC2(i,1)];
w_gpu[__I4_W(i,1,k,4)] = 0.0;
ee = get_e_from_temperature_dev_bcwall_airfoil_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,twall,t0,cv_coeff_gpu);
w_gpu[__I4_W(i,1,k,5)] = w_gpu[__I4_W(i,1,k,1)]*ee+w_gpu[__I4_W(i,1,k,1)]*0.50*(tww*tww);
}else {
w_gpu[__I4_W(i,1,k,2)] = 0.0;
w_gpu[__I4_W(i,1,k,3)] = 0.0;
w_gpu[__I4_W(i,1,k,4)] = 0.0;
ee = get_e_from_temperature_dev_bcwall_airfoil_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,twall,t0,cv_coeff_gpu);
w_gpu[__I4_W(i,1,k,5)] = w_gpu[__I4_W(i,1,k,1)]*ee;
}
for(int l=1; l<ng+1; l++){
rho = w_gpu[__I4_W(i,1+l,k,1)];
uu = w_gpu[__I4_W(i,1+l,k,2)]/rho;
vv = w_gpu[__I4_W(i,1+l,k,3)]/rho;
ww = w_gpu[__I4_W(i,1+l,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,1+l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcwall_airfoil_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,1+l,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
tt = 2.0*twall-tt;
ee = get_e_from_temperature_dev_bcwall_airfoil_kernel_3(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
rho = pp/tt/rgas0;
ug = 2.0*w_gpu[__I4_W(i,1,k,2)]/w_gpu[__I4_W(i,1,k,1)] - uu;
vg = 2.0*w_gpu[__I4_W(i,1,k,3)]/w_gpu[__I4_W(i,1,k,1)] - vv;
wg = 2.0*w_gpu[__I4_W(i,1,k,4)]/w_gpu[__I4_W(i,1,k,1)] - ww;
qq = 0.50*(ug*ug+vg*vg+wg*wg);
w_gpu[__I4_W(i,1-l,k,1)] = rho;
w_gpu[__I4_W(i,1-l,k,2)] = rho*ug;
w_gpu[__I4_W(i,1-l,k,3)] = rho*vg;
w_gpu[__I4_W(i,1-l,k,4)] = rho*wg;
w_gpu[__I4_W(i,1-l,k,5)] = rho*(ee+qq);
}
}

}
}


extern "C"{
void bcwall_airfoil_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real u0,real twall,real t0,real rgas0,real tol_iter_nr,real a_tw,real v_bs,real thic,real kx_tw,real om_tw,real time,real xtw1,real xtw2,int *wall_tag_gpu,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu,real *xc2_gpu,real *yc2_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bcwall_airfoil_kernel),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat,u0,twall,t0,rgas0,tol_iter_nr,a_tw,v_bs,thic,kx_tw,om_tw,time,xtw1,xtw2,wall_tag_gpu,w_gpu,w_aux_gpu,cv_coeff_gpu,xc2_gpu,yc2_gpu,dxdetanc2_gpu,dydetanc2_gpu);
}
}

__device__ real get_temperature_from_e_dev_bcwall_staggered_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcwall_staggered_kernel1_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_staggered_kernel1_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_staggered_kernel1_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcwall_staggered_kernel1(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
//Kernel for bcwall_staggered_kernel1
real rho;real uu;real vv;
real ww;real qq;real pp;
real tt;real rhoe;real ee;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
l = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(l,ng,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(i,l,k,1)];
uu = w_gpu[__I4_W(i,l,k,2)]/rho;
vv = w_gpu[__I4_W(i,l,k,3)]/rho;
ww = w_gpu[__I4_W(i,l,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcwall_staggered_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,l,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
tt = 2.0*twall-tt;
ee = get_e_from_temperature_dev_bcwall_staggered_kernel1_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
rho = pp/tt/rgas0;
w_gpu[__I4_W(i,1-l,k,1)] = rho;
w_gpu[__I4_W(i,1-l,k,2)] = -rho*uu;
w_gpu[__I4_W(i,1-l,k,3)] = -rho*vv;
w_gpu[__I4_W(i,1-l,k,4)] = -rho*ww;
w_gpu[__I4_W(i,1-l,k,5)] = rho*(ee+qq);

}
}


extern "C"{
void bcwall_staggered_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcwall_staggered_kernel1),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu);
}
}
__device__ real get_temperature_from_e_dev_bcwall_staggered_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_bcwall_staggered_kernel2_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}

__device__ real get_e_from_temperature_dev_bcwall_staggered_kernel2_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real tt,real *cv_coeff_gpu){
//Device kernel for get_e_from_temperature_dev_bcwall_staggered_kernel2_0
real get_e_from_temperature_dev;real ee;
int l;


ee = cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)];
if (calorically_perfect==1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(0)]*(tt/t0-1.0);
}else {
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
ee = ee+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(tt/t0);
}else {
ee=ee+cv_coeff_gpu[__I1_CV_COEFF(l)]/(l+1.0)*(pow((tt/t0),(l+1))-1.0);
}
}
}
ee = ee*t0;
get_e_from_temperature_dev = ee;


return get_e_from_temperature_dev;
}




__global__ void  bcwall_staggered_kernel2(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
//Kernel for bcwall_staggered_kernel2
real rho;real uu;real vv;
real ww;real qq;real pp;
real tt;real rhoe;real ee;
int i;int k;int l;
int iercuda;

i = __GIDX(x,1);
l = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(l,ng,1)&&loop_cond(k,nz,1)){
rho = w_gpu[__I4_W(i,ny+1-l,k,1)];
uu = w_gpu[__I4_W(i,ny+1-l,k,2)]/rho;
vv = w_gpu[__I4_W(i,ny+1-l,k,3)]/rho;
ww = w_gpu[__I4_W(i,ny+1-l,k,4)]/rho;
rhoe = w_gpu[__I4_W(i,ny+1-l,k,5)];
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_bcwall_staggered_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,ny+1-l,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
tt = 2.0*twall-tt;
ee = get_e_from_temperature_dev_bcwall_staggered_kernel2_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,tt,cv_coeff_gpu);
rho = pp/tt/rgas0;
w_gpu[__I4_W(i,ny+l,k,1)] = rho;
w_gpu[__I4_W(i,ny+l,k,2)] = -rho*uu;
w_gpu[__I4_W(i,ny+l,k,3)] = -rho*vv;
w_gpu[__I4_W(i,ny+l,k,4)] = -rho*ww;
w_gpu[__I4_W(i,ny+l,k,5)] = rho*(ee+qq);

}
}


extern "C"{
void bcwall_staggered_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,int ilat,real twall,real t0,real rgas0,real tol_iter_nr,real *w_gpu,real *w_aux_gpu,real *cv_coeff_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ng)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((bcwall_staggered_kernel2),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,ilat,twall,t0,rgas0,tol_iter_nr,w_gpu,w_aux_gpu,cv_coeff_gpu);
}
}




__global__ void  compute_residual_kernel_residual_rhou(int nx,int ny,int nz,int ng,int nv,real dt,int *fluid_mask_gpu,real *fln_gpu,real *residual_rhou,real *redn_3d_gpu){
//Kernel for compute_residual_kernel_residual_rhou
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=+(((fln_gpu[__I4_FLN(i,j,k,2)]/dt))*((fln_gpu[__I4_FLN(i,j,k,2)]/dt)));
}

}
}



extern "C"{
void compute_residual_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real dt,int *fluid_mask_gpu,real *fln_gpu,real *residual_rhou,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_residual_kernel_residual_rhou),grid,block,0,stream,nx,ny,nz,ng,nv,dt,fluid_mask_gpu,fln_gpu,residual_rhou,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, residual_rhou);


}
}




__global__ void  compute_airfoil_forces_runtime_kernel_n(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_n
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +dn*ds;
}

}
}


__global__ void  compute_airfoil_forces_runtime_kernel_a(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_a
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +da*ds;
}

}
}


__global__ void  compute_airfoil_forces_runtime_kernel_pn(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_pn
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +dpn*ds;
}

}
}


__global__ void  compute_airfoil_forces_runtime_kernel_pa(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_pa
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +dpa*ds;
}

}
}


__global__ void  compute_airfoil_forces_runtime_kernel_tn(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_tn
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +dtn*ds;
}

}
}


__global__ void  compute_airfoil_forces_runtime_kernel_ta(int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
//Kernel for compute_airfoil_forces_runtime_kernel_ta
real dudy;real dudyw;real tauw;
real pw;real cf;real cp;
real pwf;real tauwf;real ds;
real costh;real sinth;real dn;
real da;real dpn;real dpa;
real dtn;real dta;real ut1;
real ut2;real ut3;real ut4;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,1,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1 + 36.0*ut2 - 18.0*ut3 + 4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
tauw = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
pw = w_aux_gpu[__I4_W_AUX(i,1,k,1)]*w_aux_gpu[__I4_W_AUX(i,1,k,6)]*rgas0;
cf = tauw/(0.50*u0*u0);
cp = (pw-p0)/(0.50*u0*u0);
pwf = pw;
tauwf = tauw;
dn = -costh*pwf+sinth*tauwf;
da = +sinth*pwf+costh*tauwf;
dpn = -costh*pwf;
dpa = +sinth*pwf;
dtn = sinth*tauwf;
dta = costh*tauwf;
if(wall_tag_gpu[__I1_WALL_TAG(i)] < -1) ds = ds/2.0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = +dta*ds;
}

}
}



extern "C"{
void compute_airfoil_forces_runtime_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int nv,real p0,real u0,real rgas0,int *wall_tag_gpu,real *w_aux_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu,real *n,real *a,real *pn,real *pa,real *tn,real *ta,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((1)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_n),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, n);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_a),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, a);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_pn),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, pn);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_pa),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, pa);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_tn),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, tn);

hipLaunchKernelGGL((compute_airfoil_forces_runtime_kernel_ta),grid,block,0,stream,nx,ny,nz,ng,nv,p0,u0,rgas0,wall_tag_gpu,w_aux_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu,n,a,pn,pa,tn,ta,redn_3d_gpu);
reduce<real, reduce_op_add>(redn_3d_gpu, nz*ny*nx, ta);


}
}




__global__ void  compute_rho_t_p_minmax_kernel1_rhomin(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomin,real *tmin,real *pmin,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel1_rhomin
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = rho;
}

}
}


__global__ void  compute_rho_t_p_minmax_kernel1_tmin(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomin,real *tmin,real *pmin,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel1_tmin
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = tt;
}

}
}


__global__ void  compute_rho_t_p_minmax_kernel1_pmin(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomin,real *tmin,real *pmin,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel1_pmin
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = pp;
}

}
}



extern "C"{
void compute_rho_t_p_minmax_kernel1_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomin,real *tmin,real *pmin,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel1_rhomin),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomin,tmin,pmin,redn_3d_gpu);
reduce<real, reduce_op_min>(redn_3d_gpu, nz*ny*nx, rhomin);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel1_tmin),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomin,tmin,pmin,redn_3d_gpu);
reduce<real, reduce_op_min>(redn_3d_gpu, nz*ny*nx, tmin);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel1_pmin),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomin,tmin,pmin,redn_3d_gpu);
reduce<real, reduce_op_min>(redn_3d_gpu, nz*ny*nx, pmin);


}
}



__global__ void  compute_rho_t_p_minmax_kernel2_rhomax(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomax,real *tmax,real *pmax,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel2_rhomax
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = rho;
}

}
}


__global__ void  compute_rho_t_p_minmax_kernel2_tmax(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomax,real *tmax,real *pmax,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel2_tmax
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = tt;
}

}
}


__global__ void  compute_rho_t_p_minmax_kernel2_pmax(int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomax,real *tmax,real *pmax,real *redn_3d_gpu){
//Kernel for compute_rho_t_p_minmax_kernel2_pmax
real rho;real tt;real pp;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
pp = rho*tt*rgas0;
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = pp;
}

}
}



extern "C"{
void compute_rho_t_p_minmax_kernel2_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,real rgas0,int *fluid_mask_gpu,real *w_aux_gpu,real *rhomax,real *tmax,real *pmax,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel2_rhomax),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomax,tmax,pmax,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, rhomax);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel2_tmax),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomax,tmax,pmax,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, tmax);

hipLaunchKernelGGL((compute_rho_t_p_minmax_kernel2_pmax),grid,block,0,stream,nx,ny,nz,ng,rgas0,fluid_mask_gpu,w_aux_gpu,rhomax,tmax,pmax,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, pmax);


}
}

__device__ real get_gamloc_dev_compute_dt_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_compute_dt_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}





__global__ void  compute_dt_kernel_dtxi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtxi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxi;
}

}
}


__global__ void  compute_dt_kernel_dtyi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtyi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyi;
}

}
}


__global__ void  compute_dt_kernel_dtzi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtzi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzi;
}

}
}


__global__ void  compute_dt_kernel_dtxv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtxv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxv;
}

}
}


__global__ void  compute_dt_kernel_dtyv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtyv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyv;
}

}
}


__global__ void  compute_dt_kernel_dtzv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtzv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzv;
}

}
}


__global__ void  compute_dt_kernel_dtxk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtxk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxk;
}

}
}


__global__ void  compute_dt_kernel_dtyk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtyk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyk;
}

}
}


__global__ void  compute_dt_kernel_dtzk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_kernel_dtzk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
dtxi = (abs(uu)+c)*dcsidx_gpu[__I1_DCSIDX(i)];
dtyi = (abs(vv)+c)*detady_gpu[__I1_DETADY(j)];
dtzi = (abs(ww)+c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyv = nu*detadys_gpu[__I1_DETADYS(j)];
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*dcsidxs_gpu[__I1_DCSIDXS(i)];
dtyk = k_over_rhocp*detadys_gpu[__I1_DETADYS(j)];
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzk;
}

}
}



extern "C"{
void compute_dt_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *dcsidxs_gpu,real *detadys_gpu,real *dzitdzs_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_dt_kernel_dtxi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxi_max);

hipLaunchKernelGGL((compute_dt_kernel_dtyi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyi_max);

hipLaunchKernelGGL((compute_dt_kernel_dtzi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzi_max);

hipLaunchKernelGGL((compute_dt_kernel_dtxv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxv_max);

hipLaunchKernelGGL((compute_dt_kernel_dtyv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyv_max);

hipLaunchKernelGGL((compute_dt_kernel_dtzv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzv_max);

hipLaunchKernelGGL((compute_dt_kernel_dtxk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxk_max);

hipLaunchKernelGGL((compute_dt_kernel_dtyk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyk_max);

hipLaunchKernelGGL((compute_dt_kernel_dtzk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,dcsidxs_gpu,detadys_gpu,dzitdzs_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzk_max);


}
}

__device__ real get_gamloc_dev_compute_dt_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_compute_dt_c2_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}





__global__ void  compute_dt_c2_kernel_dtxi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtxi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxi;
}

}
}


__global__ void  compute_dt_c2_kernel_dtyi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtyi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyi;
}

}
}


__global__ void  compute_dt_c2_kernel_dtzi_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtzi_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzi;
}

}
}


__global__ void  compute_dt_c2_kernel_dtxv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtxv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxv;
}

}
}


__global__ void  compute_dt_c2_kernel_dtyv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtyv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyv;
}

}
}


__global__ void  compute_dt_c2_kernel_dtzv_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtzv_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzv;
}

}
}


__global__ void  compute_dt_c2_kernel_dtxk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtxk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtxk;
}

}
}


__global__ void  compute_dt_c2_kernel_dtyk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtyk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtyk;
}

}
}


__global__ void  compute_dt_c2_kernel_dtzk_max(int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
//Kernel for compute_dt_c2_kernel_dtzk_max
real dtxi;real dtyi;real dtzi;
real dtxv;real dtyv;real dtzv;
real dtxk;real dtyk;real dtzk;
real rho;real ri;real uu;
real vv;real ww;real tt;
real mu;real nu;real k_over_rhocp;
real c;real gamloc;real csis1;
real etas1;real util;real vtil;
int i;int j;int k;
int ll;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
redn_3d_gpu[__I3_REDN_3D(i,j,k)]=0.0;
if (fluid_mask_gpu[__I3_FLUID_MASK(i,j,k)]==0) {
rho = w_gpu[__I4_W(i,j,k,1)];
ri = 1.0/rho;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
mu = w_aux_gpu[__I4_W_AUX(i,j,k,7)];
nu = ri*mu;
k_over_rhocp = nu/prandtl;
gamloc = get_gamloc_dev_compute_dt_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
k_over_rhocp = nu/prandtl*(gamloc)/(gamloc-1);
c = sqrt (gamloc*rgas0*tt);
csis1=((mcsi_gpu[__I2_MCSI(i,j)])*(mcsi_gpu[__I2_MCSI(i,j)]));
etas1=((meta_gpu[__I2_META(i,j)])*(meta_gpu[__I2_META(i,j)]));
util = uu*dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)]+vv*dcsidync2_gpu[__I2_DCSIDYNC2(i,j)];
vtil = uu*detadxnc2_gpu[__I2_DETADXNC2(i,j)]+vv*detadync2_gpu[__I2_DETADYNC2(i,j)];
dtxi = (abs(util)+c)*mcsi_gpu[__I2_MCSI(i,j)];
dtyi = (abs(vtil)+c)*meta_gpu[__I2_META(i,j)];
dtzi = (abs(ww) +c)*dzitdz_gpu[__I1_DZITDZ(k)];
dtxv = nu*csis1;
dtyv = nu*etas1;
dtzv = nu*dzitdzs_gpu[__I1_DZITDZS(k)];
dtxk = k_over_rhocp*csis1;
dtyk = k_over_rhocp*etas1;
dtzk = k_over_rhocp*dzitdzs_gpu[__I1_DZITDZS(k)];
redn_3d_gpu[__I3_REDN_3D(i,j,k)] = dtzk;
}

}
}



extern "C"{
void compute_dt_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real prandtl,int *fluid_mask_gpu,real *w_gpu,real *w_aux_gpu,real *cp_coeff_gpu,real *dzitdz_gpu,real *dzitdzs_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *mcsi_gpu,real *meta_gpu,real *dtxi_max,real *dtyi_max,real *dtzi_max,real *dtxv_max,real *dtyv_max,real *dtzv_max,real *dtxk_max,real *dtyk_max,real *dtzk_max,real *redn_3d_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));


dim3 block0(THREE_X,THREE_Y,THREE_Z);
dim3 grid0(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));
hipLaunchKernelGGL((reduce_init_kernel),grid0,block0,0,stream,nx,ny,nz,redn_3d_gpu);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtxi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxi_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtyi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyi_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtzi_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzi_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtxv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxv_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtyv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyv_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtzv_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzv_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtxk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtxk_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtyk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtyk_max);

hipLaunchKernelGGL((compute_dt_c2_kernel_dtzk_max),grid,block,0,stream,nx,ny,nz,ng,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,prandtl,fluid_mask_gpu,w_gpu,w_aux_gpu,cp_coeff_gpu,dzitdz_gpu,dzitdzs_gpu,dcsidxnc2_gpu,dcsidync2_gpu,detadxnc2_gpu,detadync2_gpu,mcsi_gpu,meta_gpu,dtxi_max,dtyi_max,dtzi_max,dtxv_max,dtyv_max,dtzv_max,dtxk_max,dtyk_max,dtzk_max,redn_3d_gpu);
reduce<real, reduce_op_max>(redn_3d_gpu, nz*ny*nx, dtzk_max);


}
}

__device__ real get_temperature_from_e_dev_eval_aux_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t_start,real t0,real tol_iter_nr,real ee,real *cv_coeff_gpu){
//Device kernel for get_temperature_from_e_dev_eval_aux_kernel_0
real get_temperature_from_e_dev;real tt;real t_old;
real ebar;real den;real num;
real t_pow;real t_powp;
int l;int iter;int max_iter;


max_iter = 50;
ebar = ee - cv_coeff_gpu[__I1_CV_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+ebar/cv_coeff_gpu[__I1_CV_COEFF(0)];
}else {
t_old = t_start;
for(int iter=1; iter<max_iter+1; iter++){
den = 0.0;
num = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
if (l==-1) {
t_pow=pow((t_old/t0),l);
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*log(t_old/t0);
}else {
t_pow=pow((t_old/t0),l);
t_powp = (t_old/t0)*t_pow;
den = den+cv_coeff_gpu[__I1_CV_COEFF(l)]*t_pow;
num = num+cv_coeff_gpu[__I1_CV_COEFF(l)]*(t_powp-1.0)/(l+1.0);
}
}
num = num*t0;
tt = t_old+(ebar-num)/den;
if (abs(tt-t_old) < tol_iter_nr)  break;
t_old = tt;
}
}
get_temperature_from_e_dev = tt;


return get_temperature_from_e_dev;
}




__global__ void  eval_aux_kernel(int nx,int ny,int nz,int ng,int visc_model,int istart,int iend,int jstart,int jend,int kstart,int kend,int visc_power,int visc_sutherland,int visc_no,int indx_cp_l,int indx_cp_r,int calorically_perfect,real mu0,real t0,real sutherland_s,real t_ref_dim,real powerlaw_vtexp,real rgas0,real tol_iter_nr,real *cv_coeff_gpu,real *w_gpu,real *w_aux_gpu){
//Kernel for eval_aux_kernel
real rho;real rhou;real rhov;
real rhow;real rhoe;real ri;
real uu;real vv;real ww;
real qq;real pp;real tt;
real mu;real ee;
int i;int j;int k;
int iercuda;

i = __GIDX(x,istart);
j = __GIDX(y,jstart);
k = __GIDX(z,kstart);


if(loop_cond(i,iend,1)&&loop_cond(j,jend,1)&&loop_cond(k,kend,1)){
rho = w_gpu[__I4_W(i,j,k,1)];
rhou = w_gpu[__I4_W(i,j,k,2)];
rhov = w_gpu[__I4_W(i,j,k,3)];
rhow = w_gpu[__I4_W(i,j,k,4)];
rhoe = w_gpu[__I4_W(i,j,k,5)];
ri = 1.0/rho;
uu = rhou*ri;
vv = rhov*ri;
ww = rhow*ri;
qq = 0.50*(uu*uu+vv*vv+ww*ww);
ee = rhoe/rho-qq;
tt = get_temperature_from_e_dev_eval_aux_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,w_aux_gpu[__I4_W_AUX(i,j,k,6)],t0,tol_iter_nr,ee,cv_coeff_gpu);
pp = rho*tt*rgas0;
w_aux_gpu[__I4_W_AUX(i,j,k,1)] = rho;
w_aux_gpu[__I4_W_AUX(i,j,k,2)] = uu;
w_aux_gpu[__I4_W_AUX(i,j,k,3)] = vv;
w_aux_gpu[__I4_W_AUX(i,j,k,4)] = ww;
w_aux_gpu[__I4_W_AUX(i,j,k,5)] = (rhoe+pp)/rho;
w_aux_gpu[__I4_W_AUX(i,j,k,6)] = tt;
if (visc_model == visc_power) {
mu=mu0*pow((tt/t0),powerlaw_vtexp);
}else if (visc_model == visc_sutherland) {
mu=mu0*pow((tt/t0),1.50)*(1.0+sutherland_s/t_ref_dim)/(tt/t0+sutherland_s/t_ref_dim);
}else if (visc_model == visc_no) {
mu = 0.0;
}
w_aux_gpu[__I4_W_AUX(i,j,k,7)] = mu;

}
}


extern "C"{
void eval_aux_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_model,int istart,int iend,int jstart,int jend,int kstart,int kend,int visc_power,int visc_sutherland,int visc_no,int indx_cp_l,int indx_cp_r,int calorically_perfect,real mu0,real t0,real sutherland_s,real t_ref_dim,real powerlaw_vtexp,real rgas0,real tol_iter_nr,real *cv_coeff_gpu,real *w_gpu,real *w_aux_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((iend)-(istart)+1,block.x),divideAndRoundUp((jend)-(jstart)+1,block.y),divideAndRoundUp((kend)-(kstart)+1,block.z));

hipLaunchKernelGGL((eval_aux_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_model,istart,iend,jstart,jend,kstart,kend,visc_power,visc_sutherland,visc_no,indx_cp_l,indx_cp_r,calorically_perfect,mu0,t0,sutherland_s,t_ref_dim,powerlaw_vtexp,rgas0,tol_iter_nr,cv_coeff_gpu,w_gpu,w_aux_gpu);
}
}



__global__ void  tripping_pressure_kernel(int nx,int ny,int nz,int nv,int ng,int itr1,int itr2,int i_rank_start,real pi,real x0tr,real y0tr,real x0ts,real y0ts,real lamx,real lamy,real lamz,real lamz1,real phiz,real phiz1,real asl,real bt,int *wall_tag_gpu,real *w_gpu,real *fl_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu,real *xc2_gpu,real *yc2_gpu,real *z_gpu){
//Kernel for tripping_pressure_kernel
real xx;real yy;real zz;
real hzi;real hzi1;real gzt;
real fz;real fzx;real fzy;
int i;int j;int k;
int iercuda;int ig;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
ig = i + i_rank_start;
if (ig > itr1 && ig < itr2) {
xx = xc2_gpu[__I2_XC2(i,j)];
yy = yc2_gpu[__I2_YC2(i,j)];
zz = z_gpu[__I1_Z(k)];
hzi = sin(2.0*pi*lamz *zz+phiz );
hzi1 = sin(2.0*pi*lamz1*zz+phiz1);
gzt = asl*((1-bt)*hzi+bt*hzi1);
fz=gzt*exp(-((((xx-x0tr)/lamx))*(((xx-x0tr)/lamx)))-((((yy-y0tr)/lamy))*(((yy-y0tr)/lamy))));
fzx = 0.0;
fzy = fz;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - fzx*w_gpu[__I4_W(i,j,k,1)];
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - fzy*w_gpu[__I4_W(i,j,k,1)];
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] -(fzx*w_gpu[__I4_W(i,j,k,2)]+fzy*w_gpu[__I4_W(i,j,k,3)]);
}

}
}


extern "C"{
void tripping_pressure_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int itr1,int itr2,int i_rank_start,real pi,real x0tr,real y0tr,real x0ts,real y0ts,real lamx,real lamy,real lamz,real lamz1,real phiz,real phiz1,real asl,real bt,int *wall_tag_gpu,real *w_gpu,real *fl_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu,real *xc2_gpu,real *yc2_gpu,real *z_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((tripping_pressure_kernel),grid,block,0,stream,nx,ny,nz,nv,ng,itr1,itr2,i_rank_start,pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lamz,lamz1,phiz,phiz1,asl,bt,wall_tag_gpu,w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu);
}
}



__global__ void  tripping_suction_kernel(int nx,int ny,int nz,int nv,int ng,int its1,int its2,int i_rank_start,real pi,real x0tr,real y0tr,real x0ts,real y0ts,real lamx,real lamy,real lams,real lams1,real phis,real phis1,real asl,real bt,int *wall_tag_gpu,real *w_gpu,real *fl_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu,real *xc2_gpu,real *yc2_gpu,real *z_gpu){
//Kernel for tripping_suction_kernel
real xx;real yy;real zz;
real hzi;real hzi1;real gzt;
real fz;real fzx;real fzy;
int i;int j;int k;
int iercuda;int ig;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
ig = i + i_rank_start;
if (ig > its1 && ig < its2) {
xx = xc2_gpu[__I2_XC2(i,j)];
yy = yc2_gpu[__I2_YC2(i,j)];
zz = z_gpu[__I1_Z(k)];
hzi = sin(2.0*pi*lams *zz+phis);
hzi1 = sin(2.0*pi*lams1*zz+phis1);
gzt = asl*((1-bt)*hzi+bt*hzi1);
fz=gzt*exp(-((((xx-x0ts)/lamx))*(((xx-x0ts)/lamx)))-((((yy-y0ts)/lamy))*(((yy-y0ts)/lamy))));
fzx = 0.0;
fzy = fz;
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] - fzx*w_gpu[__I4_W(i,j,k,1)];
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] - fzy*w_gpu[__I4_W(i,j,k,1)];
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] -(fzx*w_gpu[__I4_W(i,j,k,2)]+fzy*w_gpu[__I4_W(i,j,k,3)]);
}

}
}


extern "C"{
void tripping_suction_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int nv,int ng,int its1,int its2,int i_rank_start,real pi,real x0tr,real y0tr,real x0ts,real y0ts,real lamx,real lamy,real lams,real lams1,real phis,real phis1,real asl,real bt,int *wall_tag_gpu,real *w_gpu,real *fl_gpu,real *dxdetanc2_gpu,real *dydetanc2_gpu,real *xc2_gpu,real *yc2_gpu,real *z_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((tripping_suction_kernel),grid,block,0,stream,nx,ny,nz,nv,ng,its1,its2,i_rank_start,pi,x0tr,y0tr,x0ts,y0ts,lamx,lamy,lams,lams1,phis,phis1,asl,bt,wall_tag_gpu,w_gpu,fl_gpu,dxdetanc2_gpu,dydetanc2_gpu,xc2_gpu,yc2_gpu,z_gpu);
}
}



__global__ void  insitu_div_kernel(int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
//Kernel for insitu_div_kernel
real ccl;real uu;real vv;
real ww;real ux;real vy;
real wz;real div;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
ux = 0.0;
vy = 0.0;
wz = 0.0;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
ux = ux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vy = vy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
wz = wz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
}
ux = ux*dcsidx_gpu[__I1_DCSIDX(i)];
vy = vy*detady_gpu[__I1_DETADY(j)];
wz = wz*dzitdz_gpu[__I1_DZITDZ(k)];
div = ux+vy+wz;
psi_gpu[__I4_PSI(i,j,k,mpsi)] = div;

}
}


extern "C"{
void insitu_div_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((insitu_div_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu);
}
}



__global__ void  insitu_omega_kernel(int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
//Kernel for insitu_omega_kernel
real ccl;real uu;real vv;
real ww;real uy;real uz;
real vx;real vz;real wx;
real wy;real omegax;real omegay;
real omegaz;real omod2;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
vx = 0.0;
wx = 0.0;
uy = 0.0;
wy = 0.0;
uz = 0.0;
vz = 0.0;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
vx = vx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wx = wx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
uy = uy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
wy = wy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
uz = uz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vz = vz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
}
vx = vx*dcsidx_gpu[__I1_DCSIDX(i)];
wx = wx*dcsidx_gpu[__I1_DCSIDX(i)];
uy = uy*detady_gpu[__I1_DETADY(j)];
wy = wy*detady_gpu[__I1_DETADY(j)];
uz = uz*dzitdz_gpu[__I1_DZITDZ(k)];
vz = vz*dzitdz_gpu[__I1_DZITDZ(k)];
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
psi_gpu[__I4_PSI(i,j,k,mpsi)] = sqrt(omod2);

}
}


extern "C"{
void insitu_omega_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((insitu_omega_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu);
}
}



__global__ void  insitu_ducros_kernel(int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real u0,real l0,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
//Kernel for insitu_ducros_kernel
real ccl;real uu;real vv;
real ww;real ux;real uy;
real uz;real vx;real vy;
real vz;real wx;real wy;
real wz;real div;real omegax;
real omegay;real omegaz;real omod2;
int i;int j;int k;
int l;int iercuda;

i = __GIDX(x,1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
ux = 0.0;
vx = 0.0;
wx = 0.0;
uy = 0.0;
vy = 0.0;
wy = 0.0;
uz = 0.0;
vz = 0.0;
wz = 0.0;
for(int l=1; l<lmax+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)];
ux = ux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vx = vx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wx = wx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
uy = uy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
vy = vy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
wy = wy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
uz = uz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vz = vz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wz = wz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
}
ux = ux*dcsidx_gpu[__I1_DCSIDX(i)];
vx = vx*dcsidx_gpu[__I1_DCSIDX(i)];
wx = wx*dcsidx_gpu[__I1_DCSIDX(i)];
uy = uy*detady_gpu[__I1_DETADY(j)];
vy = vy*detady_gpu[__I1_DETADY(j)];
wy = wy*detady_gpu[__I1_DETADY(j)];
uz = uz*dzitdz_gpu[__I1_DZITDZ(k)];
vz = vz*dzitdz_gpu[__I1_DZITDZ(k)];
wz = wz*dzitdz_gpu[__I1_DZITDZ(k)];
div = ux+vy+wz;
omegax = wy-vz;
omegay = uz-wx;
omegaz = vx-uy;
omod2 = omegax*omegax+omegay*omegay+omegaz*omegaz;
psi_gpu[__I4_PSI(i,j,k,mpsi)]=max(-div/sqrt(omod2+((div)*(div))+(((u0/l0))*((u0/l0)))),0.0);

}
}


extern "C"{
void insitu_ducros_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int visc_order,int npsi,int mpsi,int lmax,real u0,real l0,real *w_aux_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu){
dim3 block(THREE_X,THREE_Y,THREE_Z);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((insitu_ducros_kernel),grid,block,0,stream,nx,ny,nz,ng,visc_order,npsi,mpsi,lmax,u0,l0,w_aux_gpu,coeff_deriv1_gpu,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu);
}
}



__global__ void  probe_interpolation_kernel(int num_probe,int nx,int ny,int nz,int ng,int nv_aux,int *ijk_probe_gpu,real *w_aux_gpu,real *w_aux_probe_gpu,real *probe_coeff_gpu){
//Kernel for probe_interpolation_kernel
real w1;real w2;real w3;
real w4;real w5;real w6;
int i;int j;int k;
int ii;int jj;int kk;
int l;int iercuda;

l = __GIDX(x,1);


if(loop_cond(l,num_probe,1)){
ii = ijk_probe_gpu[__I2_IJK_PROBE(1,l)];
jj = ijk_probe_gpu[__I2_IJK_PROBE(2,l)];
kk = ijk_probe_gpu[__I2_IJK_PROBE(3,l)];
w1 = 0.0;
w2 = 0.0;
w3 = 0.0;
w4 = 0.0;
w5 = 0.0;
w6 = 0.0;
for(int k=1; k<2+1; k++){
for(int j=1; j<2+1; j++){
for(int i=1; i<2+1; i++){
w1 = w1 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,1)];
w2 = w2 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,2)];
w3 = w3 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,3)];
w4 = w4 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,4)];
w5 = w5 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,5)];
w6 = w6 + probe_coeff_gpu[__I4_PROBE_COEFF(i,j,k,l)]*w_aux_gpu[__I4_W_AUX(i+ii-1,j+jj-1,k+kk-1,6)];
}
}
}
w_aux_probe_gpu[__I2_W_AUX_PROBE(1,l)] = w1;
w_aux_probe_gpu[__I2_W_AUX_PROBE(2,l)] = w2;
w_aux_probe_gpu[__I2_W_AUX_PROBE(3,l)] = w3;
w_aux_probe_gpu[__I2_W_AUX_PROBE(4,l)] = w4;
w_aux_probe_gpu[__I2_W_AUX_PROBE(5,l)] = w5;
w_aux_probe_gpu[__I2_W_AUX_PROBE(6,l)] = w6;

}
}


extern "C"{
void probe_interpolation_kernel_wrapper(hipStream_t stream,int num_probe,int nx,int ny,int nz,int ng,int nv_aux,int *ijk_probe_gpu,real *w_aux_gpu,real *w_aux_probe_gpu,real *probe_coeff_gpu){
dim3 block(ONE_X);
dim3 grid(divideAndRoundUp((num_probe)-(1)+1,block.x));

hipLaunchKernelGGL((probe_interpolation_kernel),grid,block,0,stream,num_probe,nx,ny,nz,ng,nv_aux,ijk_probe_gpu,w_aux_gpu,w_aux_probe_gpu,probe_coeff_gpu);
}
}



__global__ void  compute_tspec_kernel(int nx,int ny,int nz,int ng,int ndft,int j_slice,int i_win,int it_win,real pi,real *w_aux_gpu,real *w_tspec_gpu){
//Kernel for compute_tspec_kernel
real wei_real;real wei_imag;real pp;
real win_scale;
int i;int j;int k;
int iercuda;int ll;int l;
int nn;int n;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int l=1; l<ndft+1; l++){
n=it_win;
ll = l-1;
nn = n-1;
wei_real = cos((-2.00*pi*ll*nn)/ndft);
wei_imag = sin((-2.00*pi*ll*nn)/ndft);
win_scale = 1.0;
pp = win_scale*w_aux_gpu[__I4_W_AUX(i,j_slice,k,1)]*w_aux_gpu[__I4_W_AUX(i,j_slice,k,6)];
if(nn==0) {
w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,1)] = pp*wei_real;
w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,2)] = pp*wei_imag;
}else {
w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,1)] = w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,1)] + pp*wei_real;
w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,2)] = w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,2)] + pp*wei_imag;
}
}

}
}


extern "C"{
void compute_tspec_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int ndft,int j_slice,int i_win,int it_win,real pi,real *w_aux_gpu,real *w_tspec_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((compute_tspec_kernel),grid,block,0,stream,nx,ny,nz,ng,ndft,j_slice,i_win,it_win,pi,w_aux_gpu,w_tspec_gpu);
}
}



__global__ void  compute_psd_tspec_kernel(int nx,int ny,int nz,int ndft,int i_win,real dt_tspec,real pi,real *w_tspec_gpu,real *w_psd_tspec_gpu){
//Kernel for compute_psd_tspec_kernel
int i;int j;int k;
int iercuda;int ll;int l;
int nn;int n;

i = __GIDX(x,1);
l = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(l,ndft,1)){
w_psd_tspec_gpu[__I2_W_PSD_TSPEC(i,l)] = 0.0;
for(int k=1; k<nz+1; k++){
w_psd_tspec_gpu[__I2_W_PSD_TSPEC(i,l)]=w_psd_tspec_gpu[__I2_W_PSD_TSPEC(i,l)]+(((w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,1)])*(w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,1)]))+((w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,2)])*(w_tspec_gpu[__I5_W_TSPEC(i,k,l,i_win,2)])));
}

}
}


extern "C"{
void compute_psd_tspec_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ndft,int i_win,real dt_tspec,real pi,real *w_tspec_gpu,real *w_psd_tspec_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ndft)-(1)+1,block.y));

hipLaunchKernelGGL((compute_psd_tspec_kernel),grid,block,0,stream,nx,ny,nz,ndft,i_win,dt_tspec,pi,w_tspec_gpu,w_psd_tspec_gpu);
}
}



__global__ void  compute_wallprop_c2_kernel(int nx,int ny,int nz,int ng,int *wall_tag_gpu,real *w_aux_gpu,real *wallprop_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
//Kernel for compute_wallprop_c2_kernel
real ut1;real ut2;real ut3;
real ut4;real ds;real costh;
real sinth;real dudy;real dudyw;
int i;int j;int k;
int iercuda;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
if(wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
ds = csimod_gpu[__I2_CSIMOD(i,1)];
costh = dxdcsic2_gpu[__I2_DXDCSIC2(i,1)]/ds;
sinth = dydcsic2_gpu[__I2_DYDCSIC2(i,1)]/ds;
ut1 = w_aux_gpu[__I4_W_AUX(i,1,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,1,k,3)]*sinth;
ut2 = w_aux_gpu[__I4_W_AUX(i,2,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,2,k,3)]*sinth;
ut3 = w_aux_gpu[__I4_W_AUX(i,3,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,3,k,3)]*sinth;
ut4 = w_aux_gpu[__I4_W_AUX(i,4,k,2)]*costh+w_aux_gpu[__I4_W_AUX(i,4,k,3)]*sinth;
dudy = -22.0*ut1+36.0*ut2-18.0*ut3+4.0*ut4;
dudyw = dudy*meta_gpu[__I2_META(i,1)]/12.0;
wallprop_gpu[__I3_WALLPROP(i,k,2)] = w_aux_gpu[__I4_W_AUX(i,1,k,7)]*dudyw;
}else {
wallprop_gpu[__I3_WALLPROP(i,k,2)] = w_aux_gpu[__I4_W_AUX(i,1,k,2)];
}

}
}


extern "C"{
void compute_wallprop_c2_kernel_wrapper(hipStream_t stream,int nx,int ny,int nz,int ng,int *wall_tag_gpu,real *w_aux_gpu,real *wallprop_gpu,real *meta_gpu,real *csimod_gpu,real *dxdcsic2_gpu,real *dydcsic2_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((compute_wallprop_c2_kernel),grid,block,0,stream,nx,ny,nz,ng,wall_tag_gpu,w_aux_gpu,wallprop_gpu,meta_gpu,csimod_gpu,dxdcsic2_gpu,dydcsic2_gpu);
}
}

__device__ void wenorec_1d_euler_x_fluxes_hybrid_c2_kernel_0(int nvar,int weno_version,int iweno,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_euler_x_fluxes_hybrid_c2_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = betascale[__LI_BETASCALE(1)];
betascale[__LI_BETASCALE(3)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(3)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(1)];
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ void eigenvectors_x_c2_euler_x_fluxes_hybrid_c2_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real ut,real dcsidxav,real dcsidyav,real *el,real *er){
//Device kernel for eigenvectors_x_c2_euler_x_fluxes_hybrid_c2_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ut * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + dcsidxav * ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + dcsidyav * ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2*uu;
el[__LI_EL(3,2)] = b2*vv;
el[__LI_EL(4,2)] = b2*ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = dcsidyav * uu - dcsidxav * vv;
el[__LI_EL(2,3)] = - dcsidyav;
el[__LI_EL(3,3)] = dcsidxav;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ut * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - dcsidxav * ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - dcsidyav * ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - dcsidxav * c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = - dcsidyav;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + dcsidxav * c;
er[__LI_ER(1,3)] = vv - dcsidyav * c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = dcsidxav;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + dcsidyav * c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - ut * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = - dcsidyav * uu + dcsidxav * vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + ut * c;


}

__device__ void compute_roe_average_euler_x_fluxes_hybrid_c2_kernel_0(int nx,int ny,int nz,int ng,int ip,int indx_cp_l,int indx_cp_r,int calorically_perfect,int i,int j,int jp,int k,int kp,real rgas0,real tol_iter_nr,real t0,real &b1,real &b2,real &b3,real &c,real &ci,real &h,real &uu,real &vv,real &ww,real *w_aux_gpu,real *cp_coeff_gpu){
//Device kernel for compute_roe_average_euler_x_fluxes_hybrid_c2_kernel_0
real up;real vp;real wp;
real qqp;real hp;real r;
real rp1;real cc;real qq;
real gam;real gm1;real tt;
real gm1loc;real hbar;real ttp;
real told;real num;real den;
real gamloc;real cploc;real tpow;
real tpowp;real p_rho;real p_e;
real etot;real rho;
int ll;int iter;int max_iter;


max_iter = 50;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
up = w_aux_gpu[__I4_W_AUX(ip,jp,kp,2)];
vp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,3)];
wp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,4)];
qqp = 0.50 * (up*up +vp*vp +wp*wp);
hp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,5)];
r = w_aux_gpu[__I4_W_AUX(ip,jp,kp,1)]/w_aux_gpu[__I4_W_AUX(i,j,k,1)];
r = sqrt(r);
rho = r*w_aux_gpu[__I4_W_AUX(i,j,k,1)];
rp1 = 1.0/(r +1.0);
uu = (r*up +uu)*rp1;
vv = (r*vp +vv)*rp1;
ww = (r*wp +ww)*rp1;
h = (r*hp +h)*rp1;
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
hbar = h - qq - cp_coeff_gpu[__I1_CP_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+hbar/cp_coeff_gpu[__I1_CP_COEFF(0)];
gamloc = cp_coeff_gpu[__I1_CP_COEFF(0)]/(cp_coeff_gpu[__I1_CP_COEFF(0)]-rgas0);
}else {
tt = w_aux_gpu[__I4_W_AUX(i ,j ,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,6)];
tt = (r*ttp +tt)*rp1;
told = tt;
for(int iter=1; iter<max_iter+1; iter++){
num = 0.0;
den = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
if (ll==-1) {
tpow=pow((told/t0),ll);
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*log(told/t0);
}else {
tpow=pow((told/t0),ll);
tpowp = (told/t0)*tpow;
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*(tpowp-1.0)/(ll+1.0);
}
}
num = num*t0;
tt = told+(hbar-num)/den;
if (abs(tt-told) < tol_iter_nr)  break;
told = tt;
}
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
gamloc = cploc/(cploc-rgas0);
}
cc = gamloc*tt*rgas0;
gm1loc = gamloc-1.0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*gm1loc;
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);


}

__device__ real get_gamloc_dev_euler_x_fluxes_hybrid_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_x_fluxes_hybrid_c2_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) euler_x_fluxes_hybrid_c2_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_imin,int eul_imax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,int *lmax_tag_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidxc2_gpu,real *dcsidxnc2_gpu,real *dcsidyc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *mcsijac1_gpu){
//Kernel for euler_x_fluxes_hybrid_c2_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real uti;real vvi;
real wwi;real ppi;real enti;
real rhoi;real tti;real dcsidxi;
real dcsidyi;real uuip;real utip;
real vvip;real wwip;real ppip;
real entip;real rhoip;real ttip;
real dcsidxip;real dcsidyip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real ft7;
real uvs1;real uvs2;real uvs3;
real uvs4;real uvs5;real uvs6;
real uvs7;real uv_part;real b1;
real b2;real b3;real c;
real ci;real h;real uu;
real vv;real ww;real rho;
real pp;real wc;real gc;
real rhou;real ut;real rhov;
real rhow;real rhoh;real dcsidxav;
real dcsidyav;real gcsixav;real gcsiyav;
real dcsimod;real qq;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int ii;
int lmax;int wenorec_ord;int lmaxi;
int weno_scheme_i;int weno_size_i;int ishk;
int ll;int mm;int n;
int n2;int iii;int jjj;
int kkk;
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))
real evmax[5];
#undef __LI_EVMAX
#define __LI_EVMAX(i) (i-(1))
real fi[5];
#undef __LI_FI
#define __LI_FI(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))
real eler[5*5];
#undef __LI_ELER
#define __LI_ELER(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,+eul_imin-2+1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,eul_imax,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
weno_scheme_i = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
if(j == 1) {
weno_scheme_i = max(lmax_tag_gpu[__I1_LMAX_TAG(i)]+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
lmax = weno_scheme_i;
}
weno_size = 2*weno_scheme_i;
ishk = 0;
for(int ii=i-weno_scheme_i+1; ii<i+weno_scheme_i+1; ii++){
if (w_aux_gpu[__I4_W_AUX(ii,j,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
ft7 = 0.0;
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
uvs7 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i-m,j,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i-m,j,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i-m,j,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i-m,j,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i-m,j,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i-m,j,k,6)];
ppi = tti*rhoi*rgas0;
uti = (uui*dcsidxc2_gpu[__I2_DCSIDXC2(i-m,j)] + vvi*dcsidyc2_gpu[__I2_DCSIDYC2(i-m,j)])/jac_gpu[__I2_JAC(i-m,j)];
dcsidxi = dcsidxc2_gpu[__I2_DCSIDXC2(i-m,j)]/jac_gpu[__I2_JAC(i-m,j)];
dcsidyi = dcsidyc2_gpu[__I2_DCSIDYC2(i-m,j)]/jac_gpu[__I2_JAC(i-m,j)];
rhoip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,6)];
ppip = ttip*rhoip*rgas0;
utip = (uuip*dcsidxc2_gpu[__I2_DCSIDXC2(i-m+l,j)] + vvip*dcsidyc2_gpu[__I2_DCSIDYC2(i-m+l,j)])/jac_gpu[__I2_JAC(i-m+l,j)];
dcsidxip = dcsidxc2_gpu[__I2_DCSIDXC2(i-m+l,j)]/jac_gpu[__I2_JAC(i-m+l,j)];
dcsidyip = dcsidyc2_gpu[__I2_DCSIDYC2(i-m+l,j)]/jac_gpu[__I2_JAC(i-m+l,j)];
rhom = rhoi+rhoip;
uv_part = (uti+utip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip)*(dcsidxi+dcsidxip);
uvs7 = uvs7 + (2.0)*(ppi+ppip)*(dcsidyi+dcsidyip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
ft7 = ft7 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs7;
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((i==0 && force_zero_flux_min == 1)||(i==nx && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh2 = fh2 + 0.250*ft6;
fh3 = fh3 + 0.250*ft7;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
compute_roe_average_euler_x_fluxes_hybrid_c2_kernel_0(nx,ny,nz,ng,i+1,indx_cp_l,indx_cp_r,calorically_perfect,i,j,j,k,k,rgas0,tol_iter_nr,t0,b1,b2,b3,c,ci,h,uu,vv,ww,w_aux_gpu,cp_coeff_gpu);
dcsidxav = 0.50 * (dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)] + dcsidxnc2_gpu[__I2_DCSIDXNC2(i+1,j)]);
dcsidyav = 0.50 * (dcsidync2_gpu[__I2_DCSIDYNC2(i,j)] + dcsidync2_gpu[__I2_DCSIDYNC2(i+1,j)]);
dcsimod=1.0/sqrt(((dcsidxav)*(dcsidxav))+((dcsidyav)*(dcsidyav)));
dcsidxav = dcsidxav * dcsimod;
dcsidyav = dcsidyav * dcsimod;
ut = dcsidxav * uu + dcsidyav * vv;
eigenvectors_x_c2_euler_x_fluxes_hybrid_c2_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,ut,dcsidxav,dcsidyav,el,er);
for(int m=1; m<5+1; m++){
evmax[__LI_EVMAX(m)] = -1.0;
}
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme_i;
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(ll,j,k,3)];
ut = dcsidxnc2_gpu[__I2_DCSIDXNC2(ll,j)]*uu+dcsidync2_gpu[__I2_DCSIDYNC2(ll,j)]*vv;
tt = w_aux_gpu[__I4_W_AUX(ll,j,k,6)];
gamloc = get_gamloc_dev_euler_x_fluxes_hybrid_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evmax[__LI_EVMAX(1)] = max(abs(ut-c),evmax[__LI_EVMAX(1)]);
evmax[__LI_EVMAX(2)] = max(abs(ut ),evmax[__LI_EVMAX(2)]);
evmax[__LI_EVMAX(3)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(4)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(5)] = max(abs(ut+c),evmax[__LI_EVMAX(5)]);
}
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme_i;
rho = w_aux_gpu[__I4_W_AUX(ll,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(ll,j,k,3)];
ut = dcsidxnc2_gpu[__I2_DCSIDXNC2(ll,j)]*uu+dcsidync2_gpu[__I2_DCSIDYNC2(ll,j)]*vv;
ww = w_aux_gpu[__I4_W_AUX(ll,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(ll,j,k,5)];
pp = rho*w_aux_gpu[__I4_W_AUX(ll,j,k,6)]*rgas0;
fi[__LI_FI(1)] = rho * ut;
fi[__LI_FI(2)] = (rho * uu * ut + pp * dcsidxnc2_gpu[__I2_DCSIDXNC2(ll,j)]);
fi[__LI_FI(3)] = (rho * vv * ut + pp * dcsidync2_gpu[__I2_DCSIDYNC2(ll,j)]);
fi[__LI_FI(4)] = rho * ww * ut;
fi[__LI_FI(5)] = rho * h * ut;
for(int m=1; m<5+1; m++){
wc = 0.0;
gc = 0.0;
wc = wc + el[__LI_EL(1,m)] * rho * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
gc = gc + el[__LI_EL(1,m)] * fi[__LI_FI(1)] * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
wc = wc + el[__LI_EL(2,m)] * rho*uu * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
gc = gc + el[__LI_EL(2,m)] * fi[__LI_FI(2)] * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
wc = wc + el[__LI_EL(3,m)] * rho*vv * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
gc = gc + el[__LI_EL(3,m)] * fi[__LI_FI(3)] * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
wc = wc + el[__LI_EL(4,m)] * rho*ww * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
gc = gc + el[__LI_EL(4,m)] * fi[__LI_FI(4)] * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
wc = wc + el[__LI_EL(5,m)] * (rho*h-pp) * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
gc = gc + el[__LI_EL(5,m)] * fi[__LI_FI(5)] * mcsijac1_gpu[__I2_MCSIJAC1(ll,j)];
c = 0.50 * (gc + evmax[__LI_EVMAX(m)] * wc);
gp[__LI_GP(m,l)] = c;
gm[__LI_GM(m,l)] = gc - c;
}
}
wenorec_1d_euler_x_fluxes_hybrid_c2_kernel_0(nv,weno_version,weno_scheme_i,weno_scheme_i,rho0,u0,gp,gm,fi);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fhat_gpu[__I4_FHAT(i,j,k,m)] + er[__LI_ER(mm,m)] * fi[__LI_FI(mm)];
}
}
}

}
}


extern "C"{
void euler_x_fluxes_hybrid_c2_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_imin,int eul_imax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,int *lmax_tag_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidxc2_gpu,real *dcsidxnc2_gpu,real *dcsidyc2_gpu,real *dcsidync2_gpu,real *jac_gpu,real *mcsijac1_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((eul_imax)-(+eul_imin-2+1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_x_fluxes_hybrid_c2_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_imin,eul_imax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,lmax_tag_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidxc2_gpu,dcsidxnc2_gpu,dcsidyc2_gpu,dcsidync2_gpu,jac_gpu,mcsijac1_gpu);
}
}

__device__ void wenorec_1d_euler_x_fluxes_hybrid_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_euler_x_fluxes_hybrid_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = betascale[__LI_BETASCALE(1)];
betascale[__LI_BETASCALE(3)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(3)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(1)];
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ void eigenvectors_x_euler_x_fluxes_hybrid_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_x_euler_x_fluxes_hybrid_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + uu * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv );
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2*uu;
el[__LI_EL(3,2)] = b2*vv;
el[__LI_EL(4,2)] = b2*ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -vv;
el[__LI_EL(2,3)] = 0.0;
el[__LI_EL(3,3)] = 1.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - uu * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv );
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 0.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + c;
er[__LI_ER(1,3)] = vv;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 1.0;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - uu * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + uu * c;


}

__device__ void compute_roe_average_euler_x_fluxes_hybrid_kernel_0(int nx,int ny,int nz,int ng,int ip,int indx_cp_l,int indx_cp_r,int calorically_perfect,int i,int j,int jp,int k,int kp,real rgas0,real tol_iter_nr,real t0,real &b1,real &b2,real &b3,real &c,real &ci,real &h,real &uu,real &vv,real &ww,real *w_aux_gpu,real *cp_coeff_gpu){
//Device kernel for compute_roe_average_euler_x_fluxes_hybrid_kernel_0
real up;real vp;real wp;
real qqp;real hp;real r;
real rp1;real cc;real qq;
real gam;real gm1;real tt;
real gm1loc;real hbar;real ttp;
real told;real num;real den;
real gamloc;real cploc;real tpow;
real tpowp;real p_rho;real p_e;
real etot;real rho;
int ll;int iter;int max_iter;


max_iter = 50;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
up = w_aux_gpu[__I4_W_AUX(ip,jp,kp,2)];
vp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,3)];
wp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,4)];
qqp = 0.50 * (up*up +vp*vp +wp*wp);
hp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,5)];
r = w_aux_gpu[__I4_W_AUX(ip,jp,kp,1)]/w_aux_gpu[__I4_W_AUX(i,j,k,1)];
r = sqrt(r);
rho = r*w_aux_gpu[__I4_W_AUX(i,j,k,1)];
rp1 = 1.0/(r +1.0);
uu = (r*up +uu)*rp1;
vv = (r*vp +vv)*rp1;
ww = (r*wp +ww)*rp1;
h = (r*hp +h)*rp1;
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
hbar = h - qq - cp_coeff_gpu[__I1_CP_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+hbar/cp_coeff_gpu[__I1_CP_COEFF(0)];
gamloc = cp_coeff_gpu[__I1_CP_COEFF(0)]/(cp_coeff_gpu[__I1_CP_COEFF(0)]-rgas0);
}else {
tt = w_aux_gpu[__I4_W_AUX(i ,j ,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,6)];
tt = (r*ttp +tt)*rp1;
told = tt;
for(int iter=1; iter<max_iter+1; iter++){
num = 0.0;
den = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
if (ll==-1) {
tpow=pow((told/t0),ll);
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*log(told/t0);
}else {
tpow=pow((told/t0),ll);
tpowp = (told/t0)*tpow;
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*(tpowp-1.0)/(ll+1.0);
}
}
num = num*t0;
tt = told+(hbar-num)/den;
if (abs(tt-told) < tol_iter_nr)  break;
told = tt;
}
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
gamloc = cploc/(cploc-rgas0);
}
cc = gamloc*tt*rgas0;
gm1loc = gamloc-1.0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*gm1loc;
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);


}

__device__ real get_gamloc_dev_euler_x_fluxes_hybrid_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_x_fluxes_hybrid_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) euler_x_fluxes_hybrid_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int istart_face,int iend_face,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu){
//Kernel for euler_x_fluxes_hybrid_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real rho;real pp;real wc;
real gc;real rhou;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int ii;
int lmax;int wenorec_ord;int ishk;
int ll;int mm;int n;
int n2;
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))
real evmax[5];
#undef __LI_EVMAX
#define __LI_EVMAX(i) (i-(1))
real fi[5];
#undef __LI_FI
#define __LI_FI(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,+istart_face-1+1);
j = __GIDX(y,1);
k = __GIDX(z,1);


if(loop_cond(i,iend_face,1)&&loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
ishk = 0;
for(int ii=i-weno_scheme+1; ii<i+weno_scheme+1; ii++){
if (w_aux_gpu[__I4_W_AUX(ii,j,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i-m,j,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i-m,j,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i-m,j,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i-m,j,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i-m,j,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i-m,j,k,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi+rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (uui+uuip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(uui*ppip+uuip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i-m,j,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i-m,j,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i-m,j,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i-m,j,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i-m,j,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i-m,j,k,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi+rhoip;
uv_part = (uui+uuip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((i==0 && force_zero_flux_min == 1)||(i==nx && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh2 = fh2 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
compute_roe_average_euler_x_fluxes_hybrid_kernel_0(nx,ny,nz,ng,i+1,indx_cp_l,indx_cp_r,calorically_perfect,i,j,j,k,k,rgas0,tol_iter_nr,t0,b1,b2,b3,c,ci,h,uu,vv,ww,w_aux_gpu,cp_coeff_gpu);
eigenvectors_x_euler_x_fluxes_hybrid_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
evmax[__LI_EVMAX(m)] = -1.0;
}
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme;
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
tt = w_aux_gpu[__I4_W_AUX(ll,j,k,6)];
gamloc = get_gamloc_dev_euler_x_fluxes_hybrid_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evmax[__LI_EVMAX(1)] = max(abs(uu-c),evmax[__LI_EVMAX(1)]);
evmax[__LI_EVMAX(2)] = max(abs(uu ),evmax[__LI_EVMAX(2)]);
evmax[__LI_EVMAX(3)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(4)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(5)] = max(abs(uu+c),evmax[__LI_EVMAX(5)]);
}
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(ll,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(ll,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(ll,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(ll,j,k,5)];
rhou = rho*uu;
pp = rho*w_aux_gpu[__I4_W_AUX(ll,j,k,6)]*rgas0;
fi[__LI_FI(1)] = rhou;
fi[__LI_FI(2)] = uu * rhou + pp;
fi[__LI_FI(3)] = vv * rhou;
fi[__LI_FI(4)] = ww * rhou;
fi[__LI_FI(5)] = h * rhou;
for(int m=1; m<5+1; m++){
wc = 0.0;
gc = 0.0;
wc = wc + el[__LI_EL(1,m)] * rho;
gc = gc + el[__LI_EL(1,m)] * fi[__LI_FI(1)];
wc = wc + el[__LI_EL(2,m)] * rho*uu;
gc = gc + el[__LI_EL(2,m)] * fi[__LI_FI(2)];
wc = wc + el[__LI_EL(3,m)] * rho*vv;
gc = gc + el[__LI_EL(3,m)] * fi[__LI_FI(3)];
wc = wc + el[__LI_EL(4,m)] * rho*ww;
gc = gc + el[__LI_EL(4,m)] * fi[__LI_FI(4)];
wc = wc + el[__LI_EL(5,m)] * (rho*h-pp);
gc = gc + el[__LI_EL(5,m)] * fi[__LI_FI(5)];
c = 0.50 * (gc + evmax[__LI_EVMAX(m)] * wc);
gp[__LI_GP(m,l)] = c;
gm[__LI_GM(m,l)] = gc - c;
}
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
wenorec_1d_euler_x_fluxes_hybrid_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fi);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fhat_gpu[__I4_FHAT(i,j,k,m)] + er[__LI_ER(mm,m)] * fi[__LI_FI(mm)];
}
}
}

}
}


extern "C"{
void euler_x_fluxes_hybrid_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int istart_face,int iend_face,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((iend_face)-(+istart_face-1+1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y),divideAndRoundUp((nz)-(1)+1,block.z));

hipLaunchKernelGGL((euler_x_fluxes_hybrid_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,istart_face,iend_face,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu);
}
}

__device__ void wenorec_1d_rusanov_euler_x_fluxes_hybrid_rusanov_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_rusanov_euler_x_fluxes_hybrid_rusanov_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(3)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(2)]/u0_2;
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ real get_gamloc_dev_euler_x_fluxes_hybrid_rusanov_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_x_fluxes_hybrid_rusanov_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void  euler_x_fluxes_hybrid_rusanov_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int istart_face,int iend_face,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu){
//Kernel for euler_x_fluxes_hybrid_rusanov_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real evm;real evmax;real rhoevm;
real rho;real pp;real wc;
real gc;real rhou;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int ii;
int lmax;int wenorec_ord;int ishk;
int ll;int mm;int n;
int n2;
real fi[5];
#undef __LI_FI
#define __LI_FI(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,+istart_face-1+1);
j = __GIDX(y,1);


if(loop_cond(i,iend_face,1)&&loop_cond(j,ny,1)){
for(int k=1; k<nz+1; k++){
ishk = 0;
for(int ii=i-weno_scheme+1; ii<i+weno_scheme+1; ii++){
if (w_aux_gpu[__I4_W_AUX(ii,j,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i-m,j,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i-m,j,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i-m,j,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i-m,j,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i-m,j,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i-m,j,k,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi+rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (uui+uuip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(uui*ppip+uuip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i-m,j,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i-m,j,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i-m,j,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i-m,j,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i-m,j,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i-m,j,k,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i-m+l,j,k,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi+rhoip;
uv_part = (uui+uuip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((i==0 && force_zero_flux_min == 1)||(i==nx && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh2 = fh2 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
evmax = -1.0;
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme;
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
tt = w_aux_gpu[__I4_W_AUX(ll,j,k,6)];
gamloc = get_gamloc_dev_euler_x_fluxes_hybrid_rusanov_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evm = max(abs(uu-c),abs(uu+c));
evmax = max(evm,evmax);
}
for(int l=1; l<weno_size+1; l++){
ll = i + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(ll,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(ll,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(ll,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(ll,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(ll,j,k,5)];
rhou = rho*uu;
pp = rho*w_aux_gpu[__I4_W_AUX(ll,j,k,6)]*rgas0;
rhoevm = rho*evmax;
evm = rhou;
c = 0.50 * (evm + rhoevm);
gp[__LI_GP(1,l)] = c;
gm[__LI_GM(1,l)] = evm - c;
evm = uu * rhou + pp;
c = 0.50 * (evm + rhoevm * uu);
gp[__LI_GP(2,l)] = c;
gm[__LI_GM(2,l)] = evm - c;
evm = vv * rhou;
c = 0.50 * (evm + rhoevm * vv);
gp[__LI_GP(3,l)] = c;
gm[__LI_GM(3,l)] = evm - c;
evm = ww * rhou;
c = 0.50 * (evm + rhoevm * ww);
gp[__LI_GP(4,l)] = c;
gm[__LI_GM(4,l)] = evm - c;
evm = h * rhou;
c = 0.50 * (evm + evmax * (rho*h-pp));
gp[__LI_GP(5,l)] = c;
gm[__LI_GM(5,l)] = evm - c;
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,1)],1);
wenorec_1d_rusanov_euler_x_fluxes_hybrid_rusanov_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fi);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fi[__LI_FI(m)];
}
}
}

}
}


extern "C"{
void euler_x_fluxes_hybrid_rusanov_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int istart_face,int iend_face,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int weno_scheme,int weno_size,int weno_version,int force_zero_flux_min,int force_zero_flux_max,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *coeff_deriv1_gpu,real *dcsidx_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((iend_face)-(+istart_face-1+1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((euler_x_fluxes_hybrid_rusanov_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,istart_face,iend_face,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,weno_scheme,weno_size,weno_version,force_zero_flux_min,force_zero_flux_max,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,coeff_deriv1_gpu,dcsidx_gpu);
}
}

__device__ void wenorec_1d_euler_z_hybrid_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_euler_z_hybrid_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = betascale[__LI_BETASCALE(1)];
betascale[__LI_BETASCALE(3)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(3)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(1)];
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ void eigenvectors_z_euler_z_hybrid_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_z_euler_z_hybrid_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ww * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu );
el[__LI_EL(3,1)] = -0.50 * (b2 * vv );
el[__LI_EL(4,1)] = -0.50 * (b2 * ww + ci);
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -uu;
el[__LI_EL(2,3)] = 1.0;
el[__LI_EL(3,3)] = 0.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -vv;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 1.0;
el[__LI_EL(4,4)] = 0.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ww * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu );
el[__LI_EL(3,5)] = -0.50 * (b2 * vv );
el[__LI_EL(4,5)] = -0.50 * (b2 * ww - ci);
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 1.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu;
er[__LI_ER(1,3)] = vv;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 0.0;
er[__LI_ER(4,3)] = 1.0;
er[__LI_ER(5,3)] = vv;
er[__LI_ER(1,4)] = ww - c;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 0.0;
er[__LI_ER(5,4)] = ww + c;
er[__LI_ER(1,5)] = h - ww * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = uu;
er[__LI_ER(4,5)] = vv;
er[__LI_ER(5,5)] = h + ww * c;


}

__device__ void compute_roe_average_euler_z_hybrid_kernel_0(int nx,int ny,int nz,int ng,int kp,int indx_cp_l,int indx_cp_r,int calorically_perfect,int i,int ip,int j,int jp,int k,real rgas0,real tol_iter_nr,real t0,real &b1,real &b2,real &b3,real &c,real &ci,real &h,real &uu,real &vv,real &ww,real *w_aux_gpu,real *cp_coeff_gpu){
//Device kernel for compute_roe_average_euler_z_hybrid_kernel_0
real up;real vp;real wp;
real qqp;real hp;real r;
real rp1;real cc;real qq;
real gam;real gm1;real tt;
real gm1loc;real hbar;real ttp;
real told;real num;real den;
real gamloc;real cploc;real tpow;
real tpowp;real p_rho;real p_e;
real etot;real rho;
int ll;int iter;int max_iter;


max_iter = 50;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
up = w_aux_gpu[__I4_W_AUX(ip,jp,kp,2)];
vp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,3)];
wp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,4)];
qqp = 0.50 * (up*up +vp*vp +wp*wp);
hp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,5)];
r = w_aux_gpu[__I4_W_AUX(ip,jp,kp,1)]/w_aux_gpu[__I4_W_AUX(i,j,k,1)];
r = sqrt(r);
rho = r*w_aux_gpu[__I4_W_AUX(i,j,k,1)];
rp1 = 1.0/(r +1.0);
uu = (r*up +uu)*rp1;
vv = (r*vp +vv)*rp1;
ww = (r*wp +ww)*rp1;
h = (r*hp +h)*rp1;
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
hbar = h - qq - cp_coeff_gpu[__I1_CP_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+hbar/cp_coeff_gpu[__I1_CP_COEFF(0)];
gamloc = cp_coeff_gpu[__I1_CP_COEFF(0)]/(cp_coeff_gpu[__I1_CP_COEFF(0)]-rgas0);
}else {
tt = w_aux_gpu[__I4_W_AUX(i ,j ,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,6)];
tt = (r*ttp +tt)*rp1;
told = tt;
for(int iter=1; iter<max_iter+1; iter++){
num = 0.0;
den = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
if (ll==-1) {
tpow=pow((told/t0),ll);
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*log(told/t0);
}else {
tpow=pow((told/t0),ll);
tpowp = (told/t0)*tpow;
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*(tpowp-1.0)/(ll+1.0);
}
}
num = num*t0;
tt = told+(hbar-num)/den;
if (abs(tt-told) < tol_iter_nr)  break;
told = tt;
}
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
gamloc = cploc/(cploc-rgas0);
}
cc = gamloc*tt*rgas0;
gm1loc = gamloc-1.0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*gm1loc;
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);


}

__device__ real get_gamloc_dev_euler_z_hybrid_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_z_hybrid_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) euler_z_hybrid_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_kmin,int eul_kmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dzitdz_gpu){
//Kernel for euler_z_hybrid_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real rho;real pp;real wc;
real gc;real rhow;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int kk;
int lmax;int wenorec_ord;int ishk;
int ll;int mm;int n;
int n2;
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))
real evmax[5];
#undef __LI_EVMAX
#define __LI_EVMAX(i) (i-(1))
real fk[5];
#undef __LI_FK
#define __LI_FK(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int k=eul_kmin-1; k<eul_kmax+1; k++){
ishk = 0;
for(int kk=k-weno_scheme+1; kk<k+weno_scheme+1; kk++){
if (w_aux_gpu[__I4_W_AUX(i,j,kk,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,3)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j,k-m,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j,k-m,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j,k-m,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j,k-m,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j,k-m,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j,k-m,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi + rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (wwi+wwip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(wwi*ppip+wwip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j,k-m,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j,k-m,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j,k-m,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j,k-m,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j,k-m,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j,k-m,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi+rhoip;
uv_part = (wwi+wwip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((k==0 && force_zero_flux_min == 1)||(k==nz && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh4 = fh4 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
compute_roe_average_euler_z_hybrid_kernel_0(nx,ny,nz,ng,k+1,indx_cp_l,indx_cp_r,calorically_perfect,i,i,j,j,k,rgas0,tol_iter_nr,t0,b1,b2,b3,c,ci,h,uu,vv,ww,w_aux_gpu,cp_coeff_gpu);
eigenvectors_z_euler_z_hybrid_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
evmax[__LI_EVMAX(m)] = -1.0;
}
for(int l=1; l<weno_size+1; l++){
ll = k + l - weno_scheme;
ww = w_aux_gpu[__I4_W_AUX(i,j,ll,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,ll,6)];
gamloc = get_gamloc_dev_euler_z_hybrid_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evmax[__LI_EVMAX(1)] = max(abs(ww-c),evmax[__LI_EVMAX(1)]);
evmax[__LI_EVMAX(2)] = max(abs(ww ),evmax[__LI_EVMAX(2)]);
evmax[__LI_EVMAX(3)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(4)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(5)] = max(abs(ww+c),evmax[__LI_EVMAX(5)]);
}
for(int l=1; l<weno_size+1; l++){
ll = k + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(i,j,ll,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,ll,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,ll,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,ll,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,ll,5)];
rhow = rho*ww;
pp = rho*w_aux_gpu[__I4_W_AUX(i,j,ll,6)]*rgas0;
fk[__LI_FK(1)] = rhow;
fk[__LI_FK(2)] = uu * rhow;
fk[__LI_FK(3)] = vv * rhow;
fk[__LI_FK(4)] = ww * rhow + pp;
fk[__LI_FK(5)] = h * rhow;
for(int m=1; m<5+1; m++){
wc = 0.0;
gc = 0.0;
wc = wc + el[__LI_EL(1,m)] * rho;
gc = gc + el[__LI_EL(1,m)] * fk[__LI_FK(1)];
wc = wc + el[__LI_EL(2,m)] * rho*uu;
gc = gc + el[__LI_EL(2,m)] * fk[__LI_FK(2)];
wc = wc + el[__LI_EL(3,m)] * rho*vv;
gc = gc + el[__LI_EL(3,m)] * fk[__LI_FK(3)];
wc = wc + el[__LI_EL(4,m)] * rho*ww;
gc = gc + el[__LI_EL(4,m)] * fk[__LI_FK(4)];
wc = wc + el[__LI_EL(5,m)] * (rho*h-pp);
gc = gc + el[__LI_EL(5,m)] * fk[__LI_FK(5)];
c = 0.50 * (gc + evmax[__LI_EVMAX(m)] * wc);
gp[__LI_GP(m,l)] = c;
gm[__LI_GM(m,l)] = gc - c;
}
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,3)],1);
wenorec_1d_euler_z_hybrid_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fk);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fhat_gpu[__I4_FHAT(i,j,k,m)] + er[__LI_ER(mm,m)] * fk[__LI_FK(mm)];
}
}
}
}

}
}


extern "C"{
void euler_z_hybrid_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_kmin,int eul_kmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dzitdz_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((euler_z_hybrid_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu);
}
}

__device__ void wenorec_1d_rusanov_euler_z_hybrid_rusanov_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_rusanov_euler_z_hybrid_rusanov_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(3)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(2)]/u0_2;
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ real get_gamloc_dev_euler_z_hybrid_rusanov_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_z_hybrid_rusanov_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void  euler_z_hybrid_rusanov_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_kmin,int eul_kmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dzitdz_gpu){
//Kernel for euler_z_hybrid_rusanov_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real evm;real evmax;real rhoevm;
real rho;real pp;real wc;
real gc;real rhow;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int kk;
int lmax;int wenorec_ord;int ishk;
int ll;int mm;int n;
int n2;
real fk[5];
#undef __LI_FK
#define __LI_FK(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
for(int k=eul_kmin-1; k<eul_kmax+1; k++){
ishk = 0;
for(int kk=k-weno_scheme+1; kk<k+weno_scheme+1; kk++){
if (w_aux_gpu[__I4_W_AUX(i,j,kk,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,3)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j,k-m,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j,k-m,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j,k-m,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j,k-m,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j,k-m,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j,k-m,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi + rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (wwi+wwip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(wwi*ppip+wwip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j,k-m,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j,k-m,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j,k-m,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j,k-m,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j,k-m,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j,k-m,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j,k-m+l,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi+rhoip;
uv_part = (wwi+wwip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((k==0 && force_zero_flux_min == 1)||(k==nz && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh4 = fh4 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
evmax = -1.0;
for(int l=1; l<weno_size+1; l++){
ll = k + l - weno_scheme;
ww = w_aux_gpu[__I4_W_AUX(i,j,ll,4)];
tt = w_aux_gpu[__I4_W_AUX(i,j,ll,6)];
gamloc = get_gamloc_dev_euler_z_hybrid_rusanov_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evm = max(abs(ww-c),abs(ww+c));
evmax = max(evm,evmax);
}
for(int l=1; l<weno_size+1; l++){
ll = k + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(i,j,ll,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,ll,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,ll,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,ll,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,ll,5)];
rhow = rho*ww;
pp = rho*w_aux_gpu[__I4_W_AUX(i,j,ll,6)]*rgas0;
rhoevm = rho*evmax;
evm = rhow;
c = 0.50 * (evm + rhoevm);
gp[__LI_GP(1,l)] = c;
gm[__LI_GM(1,l)] = evm-c;
evm = uu * rhow;
c = 0.50 * (evm + rhoevm * uu);
gp[__LI_GP(2,l)] = c;
gm[__LI_GM(2,l)] = evm-c;
evm = vv * rhow;
c = 0.50 * (evm + rhoevm * vv);
gp[__LI_GP(3,l)] = c;
gm[__LI_GM(3,l)] = evm-c;
evm = ww * rhow + pp;
c = 0.50 * (evm + rhoevm * ww);
gp[__LI_GP(4,l)] = c;
gm[__LI_GM(4,l)] = evm-c;
evm = h * rhow;
c = 0.50 * (evm + evmax * (rho*h-pp));
gp[__LI_GP(5,l)] = c;
gm[__LI_GM(5,l)] = evm-c;
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,3)],1);
wenorec_1d_rusanov_euler_z_hybrid_rusanov_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fk);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fk[__LI_FK(m)];
}
}
}

}
}


extern "C"{
void euler_z_hybrid_rusanov_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_kmin,int eul_kmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *dzitdz_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((euler_z_hybrid_rusanov_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_kmin,eul_kmax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,dzitdz_gpu);
}
}

__device__ void wenorec_1d_euler_y_hybrid_c2_kernel_0(int nvar,int weno_version,int iweno,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_euler_y_hybrid_c2_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = betascale[__LI_BETASCALE(1)];
betascale[__LI_BETASCALE(3)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(3)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(1)];
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ void eigenvectors_y_c2_euler_y_hybrid_c2_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real ut,real detadxav,real detadyav,real *el,real *er){
//Device kernel for eigenvectors_y_c2_euler_y_hybrid_c2_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ut * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + detadxav * ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + detadyav * ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -detadyav * uu + detadxav * vv;
el[__LI_EL(2,3)] = +detadyav;
el[__LI_EL(3,3)] = -detadxav;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ut * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - detadxav * ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - detadyav * ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - detadxav * c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = detadyav;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + detadxav * c;
er[__LI_ER(1,3)] = vv - detadyav * c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = - detadxav;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + detadyav * c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - ut * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = detadyav * uu - detadxav * vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + ut * c;


}

__device__ void compute_roe_average_euler_y_hybrid_c2_kernel_0(int nx,int ny,int nz,int ng,int jp,int indx_cp_l,int indx_cp_r,int calorically_perfect,int i,int ip,int j,int k,int kp,real rgas0,real tol_iter_nr,real t0,real &b1,real &b2,real &b3,real &c,real &ci,real &h,real &uu,real &vv,real &ww,real *w_aux_gpu,real *cp_coeff_gpu){
//Device kernel for compute_roe_average_euler_y_hybrid_c2_kernel_0
real up;real vp;real wp;
real qqp;real hp;real r;
real rp1;real cc;real qq;
real gam;real gm1;real tt;
real gm1loc;real hbar;real ttp;
real told;real num;real den;
real gamloc;real cploc;real tpow;
real tpowp;real p_rho;real p_e;
real etot;real rho;
int ll;int iter;int max_iter;


max_iter = 50;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
up = w_aux_gpu[__I4_W_AUX(ip,jp,kp,2)];
vp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,3)];
wp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,4)];
qqp = 0.50 * (up*up +vp*vp +wp*wp);
hp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,5)];
r = w_aux_gpu[__I4_W_AUX(ip,jp,kp,1)]/w_aux_gpu[__I4_W_AUX(i,j,k,1)];
r = sqrt(r);
rho = r*w_aux_gpu[__I4_W_AUX(i,j,k,1)];
rp1 = 1.0/(r +1.0);
uu = (r*up +uu)*rp1;
vv = (r*vp +vv)*rp1;
ww = (r*wp +ww)*rp1;
h = (r*hp +h)*rp1;
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
hbar = h - qq - cp_coeff_gpu[__I1_CP_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+hbar/cp_coeff_gpu[__I1_CP_COEFF(0)];
gamloc = cp_coeff_gpu[__I1_CP_COEFF(0)]/(cp_coeff_gpu[__I1_CP_COEFF(0)]-rgas0);
}else {
tt = w_aux_gpu[__I4_W_AUX(i ,j ,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,6)];
tt = (r*ttp +tt)*rp1;
told = tt;
for(int iter=1; iter<max_iter+1; iter++){
num = 0.0;
den = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
if (ll==-1) {
tpow=pow((told/t0),ll);
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*log(told/t0);
}else {
tpow=pow((told/t0),ll);
tpowp = (told/t0)*tpow;
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*(tpowp-1.0)/(ll+1.0);
}
}
num = num*t0;
tt = told+(hbar-num)/den;
if (abs(tt-told) < tol_iter_nr)  break;
told = tt;
}
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
gamloc = cploc/(cploc-rgas0);
}
cc = gamloc*tt*rgas0;
gm1loc = gamloc-1.0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*gm1loc;
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);


}

__device__ real get_gamloc_dev_euler_y_hybrid_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_y_hybrid_c2_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) euler_y_hybrid_c2_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *wall_tag_gpu,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detadxc2_gpu,real *detadxnc2_gpu,real *detadyc2_gpu,real *detadync2_gpu,real *jac_gpu,real *metajac1_gpu){
//Kernel for euler_y_hybrid_c2_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real uti;real vvi;
real wwi;real ppi;real enti;
real rhoi;real tti;real detadxi;
real detadyi;real uuip;real utip;
real vvip;real wwip;real ppip;
real entip;real rhoip;real ttip;
real detadxip;real detadyip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real ft7;
real uvs1;real uvs2;real uvs3;
real uvs4;real uvs5;real uvs6;
real uvs7;real uv_part;real b1;
real b2;real b3;real c;
real ci;real h;real uu;
real vv;real ww;real rho;
real pp;real wc;real gc;
real rhou;real rhov;real rhow;
real rhoe;real ut;real detadxav;
real detadyav;real detamod;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int jj;
int ishk;int ll;int mm;
int lmax;int wenorec_ord;int n;
int n2;int weno_scheme_j;int iii;
int jjj;int kkk;
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))
real evmax[5];
#undef __LI_EVMAX
#define __LI_EVMAX(i) (i-(1))
real fj[5];
#undef __LI_FJ
#define __LI_FJ(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))
real eler[5*5];
#undef __LI_ELER
#define __LI_ELER(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=0; j<eul_jmax+1; j++){
lmax = lmax_base;
weno_scheme_j = weno_scheme;
if (j <= lmax) {
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) {
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
weno_scheme_j = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
}
}else {
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
weno_scheme_j = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
}
weno_size = 2*weno_scheme_j;
ishk = 0;
for(int jj=j-weno_scheme_j+1; jj<j+weno_scheme_j+1; jj++){
if (w_aux_gpu[__I4_W_AUX(i,jj,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
ft7 = 0.0;
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
uvs7 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j-m,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j-m,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j-m,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j-m,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j-m,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j-m,k,6)];
ppi = tti*rhoi*rgas0;
uti = (uui*detadxc2_gpu[__I2_DETADXC2(i,j-m)]+vvi*detadyc2_gpu[__I2_DETADYC2(i,j-m)])/jac_gpu[__I2_JAC(i,j-m)];
detadxi = detadxc2_gpu[__I2_DETADXC2(i,j-m)]/jac_gpu[__I2_JAC(i,j-m)];
detadyi = detadyc2_gpu[__I2_DETADYC2(i,j-m)]/jac_gpu[__I2_JAC(i,j-m)];
rhoip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,6)];
ppip = ttip*rhoip*rgas0;
utip = (uuip*detadxc2_gpu[__I2_DETADXC2(i,j-m+l)]+vvip*detadyc2_gpu[__I2_DETADYC2(i,j-m+l)])/jac_gpu[__I2_JAC(i,j-m+l)];
detadxip = detadxc2_gpu[__I2_DETADXC2(i,j-m+l)]/jac_gpu[__I2_JAC(i,j-m+l)];
detadyip = detadyc2_gpu[__I2_DETADYC2(i,j-m+l)]/jac_gpu[__I2_JAC(i,j-m+l)];
rhom = rhoi + rhoip;
uv_part = (uti+utip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip)*(detadxi+detadxip);
uvs7 = uvs7 + (2.0)*(ppi+ppip)*(detadyi+detadyip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
ft7 = ft7 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs7;
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((j==0 && force_zero_flux_min == 1)||(j==ny && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh2 = fh2 + 0.250*ft6;
fh3 = fh3 + 0.250*ft7;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
compute_roe_average_euler_y_hybrid_c2_kernel_0(nx,ny,nz,ng,j+1,indx_cp_l,indx_cp_r,calorically_perfect,i,i,j,k,k,rgas0,tol_iter_nr,t0,b1,b2,b3,c,ci,h,uu,vv,ww,w_aux_gpu,cp_coeff_gpu);
detadxav = 0.50 * (detadxnc2_gpu[__I2_DETADXNC2(i,j)] + detadxnc2_gpu[__I2_DETADXNC2(i,j+1)]);
detadyav = 0.50 * (detadync2_gpu[__I2_DETADYNC2(i,j)] + detadync2_gpu[__I2_DETADYNC2(i,j+1)]);
detamod=1.0/sqrt(((detadxav)*(detadxav))+((detadyav)*(detadyav)));
detadxav = detadxav * detamod;
detadyav = detadyav * detamod;
ut = detadxav * uu + detadyav * vv;
eigenvectors_y_c2_euler_y_hybrid_c2_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,ut,detadxav,detadyav,el,er);
for(int m=1; m<5+1; m++){
evmax[__LI_EVMAX(m)] = -1.0;
}
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme_j;
uu = w_aux_gpu[__I4_W_AUX(i,ll,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
ut = detadxnc2_gpu[__I2_DETADXNC2(i,ll)]*uu+detadync2_gpu[__I2_DETADYNC2(i,ll)]*vv;
tt = w_aux_gpu[__I4_W_AUX(i,ll,k,6)];
gamloc = get_gamloc_dev_euler_y_hybrid_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evmax[__LI_EVMAX(1)] = max(abs(ut-c),evmax[__LI_EVMAX(1)]);
evmax[__LI_EVMAX(2)] = max(abs(ut ),evmax[__LI_EVMAX(2)]);
evmax[__LI_EVMAX(3)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(4)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(5)] = max(abs(ut+c),evmax[__LI_EVMAX(5)]);
}
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme_j;
rho = w_aux_gpu[__I4_W_AUX(i,ll,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,ll,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
ut = detadxnc2_gpu[__I2_DETADXNC2(i,ll)]*uu+detadync2_gpu[__I2_DETADYNC2(i,ll)]*vv;
ww = w_aux_gpu[__I4_W_AUX(i,ll,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,ll,k,5)];
rhov = rho*vv;
pp = rho*w_aux_gpu[__I4_W_AUX(i,ll,k,6)]*rgas0;
fj[__LI_FJ(1)] = rho * ut;
fj[__LI_FJ(2)] = (rho * uu * ut + pp * detadxnc2_gpu[__I2_DETADXNC2(i,ll)]);
fj[__LI_FJ(3)] = (rho * vv * ut + pp * detadync2_gpu[__I2_DETADYNC2(i,ll)]);
fj[__LI_FJ(4)] = rho * ww * ut;
fj[__LI_FJ(5)] = rho * h * ut;
for(int m=1; m<5+1; m++){
wc = 0.0;
gc = 0.0;
wc = wc + el[__LI_EL(1,m)] * rho * metajac1_gpu[__I2_METAJAC1(i,ll)];
gc = gc + el[__LI_EL(1,m)] * fj[__LI_FJ(1)] * metajac1_gpu[__I2_METAJAC1(i,ll)];
wc = wc + el[__LI_EL(2,m)] * rho*uu * metajac1_gpu[__I2_METAJAC1(i,ll)];
gc = gc + el[__LI_EL(2,m)] * fj[__LI_FJ(2)] * metajac1_gpu[__I2_METAJAC1(i,ll)];
wc = wc + el[__LI_EL(3,m)] * rho*vv * metajac1_gpu[__I2_METAJAC1(i,ll)];
gc = gc + el[__LI_EL(3,m)] * fj[__LI_FJ(3)] * metajac1_gpu[__I2_METAJAC1(i,ll)];
wc = wc + el[__LI_EL(4,m)] * rho*ww * metajac1_gpu[__I2_METAJAC1(i,ll)];
gc = gc + el[__LI_EL(4,m)] * fj[__LI_FJ(4)] * metajac1_gpu[__I2_METAJAC1(i,ll)];
wc = wc + el[__LI_EL(5,m)] * (rho*h-pp) * metajac1_gpu[__I2_METAJAC1(i,ll)];
gc = gc + el[__LI_EL(5,m)] * fj[__LI_FJ(5)] * metajac1_gpu[__I2_METAJAC1(i,ll)];
c = 0.50 * (gc + evmax[__LI_EVMAX(m)] * wc);
gp[__LI_GP(m,l)] = c;
gm[__LI_GM(m,l)] = gc - c;
}
}
wenorec_1d_euler_y_hybrid_c2_kernel_0(nv,weno_version,weno_scheme_j,weno_scheme_j,rho0,u0,gp,gm,fj);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fhat_gpu[__I4_FHAT(i,j,k,m)] + er[__LI_ER(mm,m)] * fj[__LI_FJ(mm)];
}
}
}
}

}
}


extern "C"{
void euler_y_hybrid_c2_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *wall_tag_gpu,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detadxc2_gpu,real *detadxnc2_gpu,real *detadyc2_gpu,real *detadync2_gpu,real *jac_gpu,real *metajac1_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((euler_y_hybrid_c2_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,wall_tag_gpu,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detadxc2_gpu,detadxnc2_gpu,detadyc2_gpu,detadync2_gpu,jac_gpu,metajac1_gpu);
}
}

__device__ void wenorec_1d_euler_y_hybrid_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_euler_y_hybrid_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = betascale[__LI_BETASCALE(1)];
betascale[__LI_BETASCALE(3)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(3)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(1)];
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ void eigenvectors_y_euler_y_hybrid_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_y_euler_y_hybrid_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + vv * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu );
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -uu;
el[__LI_EL(2,3)] = 1.0;
el[__LI_EL(3,3)] = 0.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - vv * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu );
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 1.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu;
er[__LI_ER(1,3)] = vv - c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 0.0;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - vv * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = uu;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + vv * c;


}

__device__ void compute_roe_average_euler_y_hybrid_kernel_0(int nx,int ny,int nz,int ng,int jp,int indx_cp_l,int indx_cp_r,int calorically_perfect,int i,int ip,int j,int k,int kp,real rgas0,real tol_iter_nr,real t0,real &b1,real &b2,real &b3,real &c,real &ci,real &h,real &uu,real &vv,real &ww,real *w_aux_gpu,real *cp_coeff_gpu){
//Device kernel for compute_roe_average_euler_y_hybrid_kernel_0
real up;real vp;real wp;
real qqp;real hp;real r;
real rp1;real cc;real qq;
real gam;real gm1;real tt;
real gm1loc;real hbar;real ttp;
real told;real num;real den;
real gamloc;real cploc;real tpow;
real tpowp;real p_rho;real p_e;
real etot;real rho;
int ll;int iter;int max_iter;


max_iter = 50;
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
up = w_aux_gpu[__I4_W_AUX(ip,jp,kp,2)];
vp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,3)];
wp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,4)];
qqp = 0.50 * (up*up +vp*vp +wp*wp);
hp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,5)];
r = w_aux_gpu[__I4_W_AUX(ip,jp,kp,1)]/w_aux_gpu[__I4_W_AUX(i,j,k,1)];
r = sqrt(r);
rho = r*w_aux_gpu[__I4_W_AUX(i,j,k,1)];
rp1 = 1.0/(r +1.0);
uu = (r*up +uu)*rp1;
vv = (r*vp +vv)*rp1;
ww = (r*wp +ww)*rp1;
h = (r*hp +h)*rp1;
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
hbar = h - qq - cp_coeff_gpu[__I1_CP_COEFF(indx_cp_r+1)]*t0;
if (calorically_perfect==1) {
tt = t0+hbar/cp_coeff_gpu[__I1_CP_COEFF(0)];
gamloc = cp_coeff_gpu[__I1_CP_COEFF(0)]/(cp_coeff_gpu[__I1_CP_COEFF(0)]-rgas0);
}else {
tt = w_aux_gpu[__I4_W_AUX(i ,j ,k ,6)];
ttp = w_aux_gpu[__I4_W_AUX(ip,jp,kp,6)];
tt = (r*ttp +tt)*rp1;
told = tt;
for(int iter=1; iter<max_iter+1; iter++){
num = 0.0;
den = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
if (ll==-1) {
tpow=pow((told/t0),ll);
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*log(told/t0);
}else {
tpow=pow((told/t0),ll);
tpowp = (told/t0)*tpow;
den = den+cp_coeff_gpu[__I1_CP_COEFF(ll)]*tpow;
num = num+cp_coeff_gpu[__I1_CP_COEFF(ll)]*(tpowp-1.0)/(ll+1.0);
}
}
num = num*t0;
tt = told+(hbar-num)/den;
if (abs(tt-told) < tol_iter_nr)  break;
told = tt;
}
cploc = 0.0;
for(int ll=indx_cp_l; ll<indx_cp_r+1; ll++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(ll)]*pow((tt/t0),ll);
}
gamloc = cploc/(cploc-rgas0);
}
cc = gamloc*tt*rgas0;
gm1loc = gamloc-1.0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*gm1loc;
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);


}

__device__ real get_gamloc_dev_euler_y_hybrid_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_y_hybrid_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) euler_y_hybrid_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detady_gpu){
//Kernel for euler_y_hybrid_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real rho;real pp;real wc;
real gc;real rhov;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int jj;
int ishk;int ll;int mm;
int lmax;int wenorec_ord;int n;
int n2;
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))
real evmax[5];
#undef __LI_EVMAX
#define __LI_EVMAX(i) (i-(1))
real fj[5];
#undef __LI_FJ
#define __LI_FJ(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=eul_jmin-1; j<eul_jmax+1; j++){
ishk = 0;
for(int jj=j-weno_scheme+1; jj<j+weno_scheme+1; jj++){
if (w_aux_gpu[__I4_W_AUX(i,jj,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j-m,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j-m,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j-m,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j-m,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j-m,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j-m,k,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi + rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (vvi+vvip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(vvi*ppip+vvip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j-m,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j-m,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j-m,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j-m,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j-m,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j-m,k,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi + rhoip;
uv_part = (vvi+vvip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((j==0 && force_zero_flux_min == 1)||(j==ny && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh3 = fh3 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
compute_roe_average_euler_y_hybrid_kernel_0(nx,ny,nz,ng,j+1,indx_cp_l,indx_cp_r,calorically_perfect,i,i,j,k,k,rgas0,tol_iter_nr,t0,b1,b2,b3,c,ci,h,uu,vv,ww,w_aux_gpu,cp_coeff_gpu);
eigenvectors_y_euler_y_hybrid_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
evmax[__LI_EVMAX(m)] = -1.0;
}
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme;
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
tt = w_aux_gpu[__I4_W_AUX(i,ll,k,6)];
gamloc = get_gamloc_dev_euler_y_hybrid_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evmax[__LI_EVMAX(1)] = max(abs(vv-c),evmax[__LI_EVMAX(1)]);
evmax[__LI_EVMAX(2)] = max(abs(vv ),evmax[__LI_EVMAX(2)]);
evmax[__LI_EVMAX(3)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(4)] = evmax[__LI_EVMAX(2)];
evmax[__LI_EVMAX(5)] = max(abs(vv+c),evmax[__LI_EVMAX(5)]);
}
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(i,ll,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,ll,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,ll,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,ll,k,5)];
rhov = rho*vv;
pp = rho*w_aux_gpu[__I4_W_AUX(i,ll,k,6)]*rgas0;
fj[__LI_FJ(1)] = rhov;
fj[__LI_FJ(2)] = uu * rhov;
fj[__LI_FJ(3)] = vv * rhov + pp;
fj[__LI_FJ(4)] = ww * rhov;
fj[__LI_FJ(5)] = h * rhov;
for(int m=1; m<5+1; m++){
wc = 0.0;
gc = 0.0;
wc = wc + el[__LI_EL(1,m)] * rho;
gc = gc + el[__LI_EL(1,m)] * fj[__LI_FJ(1)];
wc = wc + el[__LI_EL(2,m)] * rho*uu;
gc = gc + el[__LI_EL(2,m)] * fj[__LI_FJ(2)];
wc = wc + el[__LI_EL(3,m)] * rho*vv;
gc = gc + el[__LI_EL(3,m)] * fj[__LI_FJ(3)];
wc = wc + el[__LI_EL(4,m)] * rho*ww;
gc = gc + el[__LI_EL(4,m)] * fj[__LI_FJ(4)];
wc = wc + el[__LI_EL(5,m)] * (rho*h-pp);
gc = gc + el[__LI_EL(5,m)] * fj[__LI_FJ(5)];
c = 0.50 * (gc + evmax[__LI_EVMAX(m)] * wc);
gp[__LI_GP(m,l)] = c;
gm[__LI_GM(m,l)] = gc - c;
}
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
wenorec_1d_euler_y_hybrid_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fj);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fhat_gpu[__I4_FHAT(i,j,k,m)] + er[__LI_ER(mm,m)] * fj[__LI_FJ(mm)];
}
}
}
}

}
}


extern "C"{
void euler_y_hybrid_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detady_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((euler_y_hybrid_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu);
}
}

__device__ void wenorec_1d_rusanov_euler_y_hybrid_rusanov_kernel_0(int nvar,int iweno,int weno_version,int wenorec_ord,real rho0,real u0,real *vp,real *vm,real *vhat){
//Device kernel for wenorec_1d_rusanov_euler_y_hybrid_rusanov_kernel_0
real vminus;real vplus;real c0;
real c1;real c2;real c3;
real c4;real d0;real d1;
real d2;real d3;real summ;
real sump;real tau5p;real tau5m;
real eps40;real u0_2;real rho0_2u0_2;
real rho0_2u0_4;
int i;int l;int m;
#undef __LI_VM
#define __LI_VM(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VP
#define __LI_VP(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_VHAT
#define __LI_VHAT(i) (i-(1))
real dwe[((4)-(-1))+1];
#undef __LI_DWE
#define __LI_DWE(i) (i-(-1))
real betap[((4)-(-1))+1];
#undef __LI_BETAP
#define __LI_BETAP(i) (i-(-1))
real betam[((4)-(-1))+1];
#undef __LI_BETAM
#define __LI_BETAM(i) (i-(-1))
real betascale[5];
#undef __LI_BETASCALE
#define __LI_BETASCALE(i) (i-(1))


u0_2 = u0*u0;
rho0_2u0_2 = rho0*rho0*u0_2;
rho0_2u0_4 = rho0_2u0_2*u0_2;
betascale[__LI_BETASCALE(1)] = 1.0/rho0_2u0_2;
betascale[__LI_BETASCALE(2)] = 1.0/rho0_2u0_4;
betascale[__LI_BETASCALE(3)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(4)] = betascale[__LI_BETASCALE(2)];
betascale[__LI_BETASCALE(5)] = betascale[__LI_BETASCALE(2)]/u0_2;
if (wenorec_ord==1) {
i = iweno;
for(int m=1; m<5+1; m++){
vminus = vp[__LI_VP(m,i)];
vplus = vm[__LI_VM(m,i+1)];
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else if (wenorec_ord==2) {
i = iweno;
dwe[__LI_DWE(1)] = 2.0/3.0;
dwe[__LI_DWE(0)] = 1.0/3.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(0)]=(((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)]))*((vp[__LI_VP(m,i)]-vp[__LI_VP(m,i-1)])));
betap[__LI_BETAP(1)]=(((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)]))*((vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i)])));
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(0)]=(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(1)]=(((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+1)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<1+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(0)] *(-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i )]) + betap[__LI_BETAP(1)] *( vp[__LI_VP(m,i )]+ vp[__LI_VP(m,i+1)]);
vplus = betam[__LI_BETAM(0)] *(-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]) + betam[__LI_BETAM(1)] *( vm[__LI_VM(m,i )]+ vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = 0.50*(vminus+vplus);
}
}else if (wenorec_ord==3) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/10.0;
dwe[__LI_DWE( 1)] = 6.0/10.0;
dwe[__LI_DWE( 2)] = 3.0/10.0;
d0 = 13.0/12.0;
d1 = 1.0/4.0;
c0 = 1.0/3.0;
c1 = 5.0/6.0;
c2 =-1.0/6.0;
c3 =-7.0/6.0;
c4 =11.0/6.0;
if (weno_version==0) {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}else {
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(2)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d0*(((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2.0*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d1*(((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d0*(((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((vp[__LI_VP(m,i)]-2.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])))+d1*(((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)]))*((3.0*vp[__LI_VP(m,i)]-4.0*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i-2)])));
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(2)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d0*(((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2.0*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d1*(((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d0*(((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((vm[__LI_VM(m,i+1)]-2.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])))+d1*(((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)]))*((3.0*vm[__LI_VM(m,i+1)]-4.0*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+3)])));
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
eps40 = 1.e-35;
tau5p = abs(betap[__LI_BETAP(0)]-betap[__LI_BETAP(2)])+eps40;
tau5m = abs(betam[__LI_BETAM(0)]-betam[__LI_BETAM(2)])+eps40;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = (betap[__LI_BETAP(l)]+eps40)/(betap[__LI_BETAP(l)]+tau5p);
betam[__LI_BETAM(l)] = (betam[__LI_BETAM(l)]+eps40)/(betam[__LI_BETAM(l)]+tau5m);
}
sump = 0.0;
summ = 0.0;
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = dwe[__LI_DWE(l)]/betap[__LI_BETAP(l)];
betam[__LI_BETAM(l)] = dwe[__LI_DWE(l)]/betam[__LI_BETAM(l)];
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<2+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(2)]*(c0*vp[__LI_VP(m,i )]+c1*vp[__LI_VP(m,i+1)]+c2*vp[__LI_VP(m,i+2)]) + betap[__LI_BETAP(1)]*(c2*vp[__LI_VP(m,i-1)]+c1*vp[__LI_VP(m,i )]+c0*vp[__LI_VP(m,i+1)]) + betap[__LI_BETAP(0)]*(c0*vp[__LI_VP(m,i-2)]+c3*vp[__LI_VP(m,i-1)]+c4*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(2)]*(c0*vm[__LI_VM(m,i+1)]+c1*vm[__LI_VM(m,i )]+c2*vm[__LI_VM(m,i-1)]) + betam[__LI_BETAM(1)]*(c2*vm[__LI_VM(m,i+2)]+c1*vm[__LI_VM(m,i+1)]+c0*vm[__LI_VM(m,i )]) + betam[__LI_BETAM(0)]*(c0*vm[__LI_VM(m,i+3)]+c3*vm[__LI_VM(m,i+2)]+c4*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = vminus+vplus;
}
}
}else if (wenorec_ord==4) {
i = iweno;
dwe[__LI_DWE( 0)] = 1.0/35.0;
dwe[__LI_DWE( 1)] = 12.0/35.0;
dwe[__LI_DWE( 2)] = 18.0/35.0;
dwe[__LI_DWE( 3)] = 4.0/35.0;
d1 = 1.0/36.0;
d2 = 13.0/12.0;
d3 = 781.0/720.0;
for(int m=1; m<5+1; m++){
betap[__LI_BETAP(3)]=d1*(((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)]))*((-11*vp[__LI_VP(m,i)]+18*vp[__LI_VP(m,i+1)]-9*vp[__LI_VP(m,i+2)]+2*vp[__LI_VP(m,i+3)])))+d2*(((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)]))*((2*vp[__LI_VP(m,i)]-5*vp[__LI_VP(m,i+1)]+4*vp[__LI_VP(m,i+2)]-vp[__LI_VP(m,i+3)])))+d3*(((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)]))*((-vp[__LI_VP(m,i)]+3*vp[__LI_VP(m,i+1)]-3*vp[__LI_VP(m,i+2)]+vp[__LI_VP(m,i+3)])));
betap[__LI_BETAP(2)]=d1*(((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)]))*((-2*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+6*vp[__LI_VP(m,i+1)]-vp[__LI_VP(m,i+2)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)]))*((-vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]-3*vp[__LI_VP(m,i+1)]+vp[__LI_VP(m,i+2)])));
betap[__LI_BETAP(1)]=d1*(((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-2)]-6*vp[__LI_VP(m,i-1)]+3*vp[__LI_VP(m,i)]+2*vp[__LI_VP(m,i+1)])))+d2*(((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((vp[__LI_VP(m,i-1)]-2*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])))+d3*(((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)]))*((-vp[__LI_VP(m,i-2)]+3*vp[__LI_VP(m,i-1)]-3*vp[__LI_VP(m,i)]+vp[__LI_VP(m,i+1)])));
betap[__LI_BETAP(0)]=d1*(((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)]))*((-2*vp[__LI_VP(m,i-3)]+9*vp[__LI_VP(m,i-2)]-18*vp[__LI_VP(m,i-1)]+11*vp[__LI_VP(m,i)])))+d2*(((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+4*vp[__LI_VP(m,i-2)]-5*vp[__LI_VP(m,i-1)]+2*vp[__LI_VP(m,i)])))+d3*(((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)]))*((-vp[__LI_VP(m,i-3)]+3*vp[__LI_VP(m,i-2)]-3*vp[__LI_VP(m,i-1)]+vp[__LI_VP(m,i)])));
betap[__LI_BETAP(3)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(3)];
betap[__LI_BETAP(2)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(2)];
betap[__LI_BETAP(1)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(1)];
betap[__LI_BETAP(0)] = betascale[__LI_BETASCALE(m)]*betap[__LI_BETAP(0)];
betam[__LI_BETAM(3)]=d1*(((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)]))*((-11*vm[__LI_VM(m,i+1)]+18*vm[__LI_VM(m,i)]-9*vm[__LI_VM(m,i-1)]+2*vm[__LI_VM(m,i-2)])))+d2*(((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)]))*((2*vm[__LI_VM(m,i+1)]-5*vm[__LI_VM(m,i)]+4*vm[__LI_VM(m,i-1)]-vm[__LI_VM(m,i-2)])))+d3*(((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)]))*((-vm[__LI_VM(m,i+1)]+3*vm[__LI_VM(m,i)]-3*vm[__LI_VM(m,i-1)]+vm[__LI_VM(m,i-2)])));
betam[__LI_BETAM(2)]=d1*(((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)]))*((-2*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+6*vm[__LI_VM(m,i)]-vm[__LI_VM(m,i-1)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)]))*((-vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]-3*vm[__LI_VM(m,i)]+vm[__LI_VM(m,i-1)])));
betam[__LI_BETAM(1)]=d1*(((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+3)]-6*vm[__LI_VM(m,i+2)]+3*vm[__LI_VM(m,i+1)]+2*vm[__LI_VM(m,i)])))+d2*(((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((vm[__LI_VM(m,i+2)]-2*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])))+d3*(((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)]))*((-vm[__LI_VM(m,i+3)]+3*vm[__LI_VM(m,i+2)]-3*vm[__LI_VM(m,i+1)]+vm[__LI_VM(m,i)])));
betam[__LI_BETAM(0)]=d1*(((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)]))*((-2*vm[__LI_VM(m,i+4)]+9*vm[__LI_VM(m,i+3)]-18*vm[__LI_VM(m,i+2)]+11*vm[__LI_VM(m,i+1)])))+d2*(((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+4*vm[__LI_VM(m,i+3)]-5*vm[__LI_VM(m,i+2)]+2*vm[__LI_VM(m,i+1)])))+d3*(((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)]))*((-vm[__LI_VM(m,i+4)]+3*vm[__LI_VM(m,i+3)]-3*vm[__LI_VM(m,i+2)]+vm[__LI_VM(m,i+1)])));
betam[__LI_BETAM(3)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(3)];
betam[__LI_BETAM(2)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(2)];
betam[__LI_BETAM(1)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(1)];
betam[__LI_BETAM(0)] = betascale[__LI_BETASCALE(m)]*betam[__LI_BETAM(0)];
sump = 0.0;
summ = 0.0;
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betap[__LI_BETAP(l)]))*((0.0000010+betap[__LI_BETAP(l)])));
betam[__LI_BETAM(l)]=dwe[__LI_DWE(l)]/(((0.0000010+betam[__LI_BETAM(l)]))*((0.0000010+betam[__LI_BETAM(l)])));
sump = sump + betap[__LI_BETAP(l)];
summ = summ + betam[__LI_BETAM(l)];
}
for(int l=0; l<3+1; l++){
betap[__LI_BETAP(l)] = betap[__LI_BETAP(l)]/sump;
betam[__LI_BETAM(l)] = betam[__LI_BETAM(l)]/summ;
}
vminus = betap[__LI_BETAP(3)]*( 6*vp[__LI_VP(m,i )]+26*vp[__LI_VP(m,i+1)]-10*vp[__LI_VP(m,i+2)]+ 2*vp[__LI_VP(m,i+3)])+betap[__LI_BETAP(2)]*(-2*vp[__LI_VP(m,i-1)]+14*vp[__LI_VP(m,i )]+14*vp[__LI_VP(m,i+1)]- 2*vp[__LI_VP(m,i+2)])+betap[__LI_BETAP(1)]*( 2*vp[__LI_VP(m,i-2)]-10*vp[__LI_VP(m,i-1)]+26*vp[__LI_VP(m,i )]+ 6*vp[__LI_VP(m,i+1)])+betap[__LI_BETAP(0)]*(-6*vp[__LI_VP(m,i-3)]+26*vp[__LI_VP(m,i-2)]-46*vp[__LI_VP(m,i-1)]+50*vp[__LI_VP(m,i )]);
vplus = betam[__LI_BETAM(3)]*( 6*vm[__LI_VM(m,i+1)]+26*vm[__LI_VM(m,i )]-10*vm[__LI_VM(m,i-1)]+ 2*vm[__LI_VM(m,i-2)])+betam[__LI_BETAM(2)]*(-2*vm[__LI_VM(m,i+2)]+14*vm[__LI_VM(m,i+1)]+14*vm[__LI_VM(m,i )]- 2*vm[__LI_VM(m,i-1)])+betam[__LI_BETAM(1)]*( 2*vm[__LI_VM(m,i+3)]-10*vm[__LI_VM(m,i+2)]+26*vm[__LI_VM(m,i+1)]+ 6*vm[__LI_VM(m,i )])+betam[__LI_BETAM(0)]*(-6*vm[__LI_VM(m,i+4)]+26*vm[__LI_VM(m,i+3)]-46*vm[__LI_VM(m,i+2)]+50*vm[__LI_VM(m,i+1)]);
vhat[__LI_VHAT(m)] = (vminus+vplus)/24.0;
}
}


}

__device__ real get_gamloc_dev_euler_y_hybrid_rusanov_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_euler_y_hybrid_rusanov_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void  euler_y_hybrid_rusanov_kernel(int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detady_gpu){
//Kernel for euler_y_hybrid_rusanov_kernel
real fh1;real fh2;real fh3;
real fh4;real fh5;real rhom;
real uui;real vvi;real wwi;
real ppi;real enti;real rhoi;
real tti;real uuip;real vvip;
real wwip;real ppip;real entip;
real rhoip;real ttip;real ft1;
real ft2;real ft3;real ft4;
real ft5;real ft6;real uvs1;
real uvs2;real uvs3;real uvs4;
real uvs6;real uv_part;real uvs5;
real b1;real b2;real b3;
real c;real ci;real h;
real uu;real vv;real ww;
real evm;real evmax;real rhoevm;
real rho;real pp;real wc;
real gc;real rhov;real tt;
real gamloc;real uvs5_i;real uvs5_k;
real uvs5_p;real eei;real eeip;
real drho;real dee;real eem;
real drhof;real deef;real sumnumrho;
real sumnumee;real sumdenrho;real sumdenee;
real t_sumdenrho;real t_sumdenee;real t2_sumdenrho;
real t2_sumdenee;
int i;int j;int k;
int m;int l;int jj;
int ishk;int ll;int mm;
int lmax;int wenorec_ord;int n;
int n2;
real fj[5];
#undef __LI_FJ
#define __LI_FJ(i) (i-(1))
real gp[5*8];
#undef __LI_GP
#define __LI_GP(i,j) ((i-(1))+5*(j-(1)))
real gm[5*8];
#undef __LI_GM
#define __LI_GM(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=eul_jmin-1; j<eul_jmax+1; j++){
ishk = 0;
for(int jj=j-weno_scheme+1; jj<j+weno_scheme+1; jj++){
if (w_aux_gpu[__I4_W_AUX(i,jj,k,8)] > sensor_threshold) ishk = 1;
}
if (ishk == 0) {
ft1 = 0.0;
ft2 = 0.0;
ft3 = 0.0;
ft4 = 0.0;
ft5 = 0.0;
ft6 = 0.0;
lmax = max(lmax_base+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
if (nkeep>=0) {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5_i = 0.0;
uvs5_k = 0.0;
uvs5_p = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j-m,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j-m,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j-m,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j-m,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j-m,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j-m,k,6)];
ppi = tti*rhoi*rgas0;
eei = enti-ppi/rhoi-0.50*(uui*uui+vvi*vvi+wwi*wwi);
rhoip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,6)];
ppip = ttip*rhoip*rgas0;
eeip = entip-ppip/rhoip-0.50*(uuip*uuip+vvip*vvip+wwip*wwip);
rhom = rhoi + rhoip;
eem = eei + eeip;
if(nkeep == 0) {
drhof = 1.0;
deef = 1.0;
}else {
sumnumrho = 1.0;
drho = 2.0*(rhoip-rhoi)/rhom;
dee = 2.0*(eeip - eei)/eem;
t_sumdenrho = (0.50*drho)*(0.50*drho);
t_sumdenee = (0.50*dee )*(0.50*dee );
t2_sumdenrho = t_sumdenrho;
t2_sumdenee = t_sumdenee;
sumdenrho = 1.0 + t_sumdenrho / (3.0);
sumdenee = 1.0 + t_sumdenee;
sumnumee = 1.0 + t_sumdenee / (3.0);
for(int n=2; n<nkeep+1; n++){
n2 = 2*n;
t_sumdenrho = t2_sumdenrho * t_sumdenrho;
t_sumdenee = t2_sumdenee * t_sumdenee;
sumdenrho = sumdenrho + t_sumdenrho / (1.0+n2);
sumdenee = sumdenee + t_sumdenee;
sumnumee = sumnumee + t_sumdenee / (1.0+n2);
}
drhof = sumnumrho/sumdenrho;
deef = sumnumee /sumdenee;
}
uv_part = (vvi+vvip) * rhom * drhof;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5_i = uvs5_i + uv_part * eem * deef;
uvs5_k = uvs5_k + uv_part * (uui*uuip+vvi*vvip+wwi*wwip);
uvs5_p = uvs5_p + 4.0*(vvi*ppip+vvip*ppi);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*(uvs5_i+uvs5_k+uvs5_p);
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}else {
for(int l=1; l<lmax+1; l++){
uvs1 = 0.0;
uvs2 = 0.0;
uvs3 = 0.0;
uvs4 = 0.0;
uvs5 = 0.0;
uvs6 = 0.0;
for(int m=0; m<l-1+1; m++){
rhoi = w_aux_gpu[__I4_W_AUX(i,j-m,k,1)];
uui = w_aux_gpu[__I4_W_AUX(i,j-m,k,2)];
vvi = w_aux_gpu[__I4_W_AUX(i,j-m,k,3)];
wwi = w_aux_gpu[__I4_W_AUX(i,j-m,k,4)];
enti = w_aux_gpu[__I4_W_AUX(i,j-m,k,5)];
tti = w_aux_gpu[__I4_W_AUX(i,j-m,k,6)];
ppi = tti*rhoi*rgas0;
rhoip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,1)];
uuip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,2)];
vvip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,3)];
wwip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,4)];
entip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,5)];
ttip = w_aux_gpu[__I4_W_AUX(i,j-m+l,k,6)];
ppip = ttip*rhoip*rgas0;
rhom = rhoi + rhoip;
uv_part = (vvi+vvip) * rhom;
uvs1 = uvs1 + uv_part * (2.0);
uvs2 = uvs2 + uv_part * (uui+uuip);
uvs3 = uvs3 + uv_part * (vvi+vvip);
uvs4 = uvs4 + uv_part * (wwi+wwip);
uvs5 = uvs5 + uv_part * (enti+entip);
uvs6 = uvs6 + (2.0)*(ppi+ppip);
}
ft1 = ft1 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs1;
ft2 = ft2 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs2;
ft3 = ft3 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs3;
ft4 = ft4 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs4;
ft5 = ft5 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs5;
ft6 = ft6 + coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmax)]*uvs6;
}
}
fh1 = 0.250*ft1;
fh2 = 0.250*ft2;
fh3 = 0.250*ft3;
fh4 = 0.250*ft4;
fh5 = 0.250*ft5;
if ((j==0 && force_zero_flux_min == 1)||(j==ny && force_zero_flux_max == 1)) {
fh1 = 0.0;
fh2 = 0.0;
fh3 = 0.0;
fh4 = 0.0;
fh5 = 0.0;
}
fh3 = fh3 + 0.50*ft6;
fhat_gpu[__I4_FHAT(i,j,k,1)] = fh1;
fhat_gpu[__I4_FHAT(i,j,k,2)] = fh2;
fhat_gpu[__I4_FHAT(i,j,k,3)] = fh3;
fhat_gpu[__I4_FHAT(i,j,k,4)] = fh4;
fhat_gpu[__I4_FHAT(i,j,k,5)] = fh5;
}else {
evmax = -1.0;
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme;
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
tt = w_aux_gpu[__I4_W_AUX(i,ll,k,6)];
gamloc = get_gamloc_dev_euler_y_hybrid_rusanov_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
c = sqrt (gamloc*rgas0*tt);
evm = max(abs(vv-c),abs(vv+c));
evmax = max(evm,evmax);
}
for(int l=1; l<weno_size+1; l++){
ll = j + l - weno_scheme;
rho = w_aux_gpu[__I4_W_AUX(i,ll,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,ll,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,ll,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,ll,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,ll,k,5)];
rhov = rho*vv;
pp = rho*w_aux_gpu[__I4_W_AUX(i,ll,k,6)]*rgas0;
rhoevm = rho*evmax;
evm = rhov;
c = 0.50 * (evm + rhoevm);
gp[__LI_GP(1,l)] = c;
gm[__LI_GM(1,l)] = evm-c;
evm = uu * rhov;
c = 0.50 * (evm + rhoevm * uu);
gp[__LI_GP(2,l)] = c;
gm[__LI_GM(2,l)] = evm-c;
evm = vv * rhov + pp;
c = 0.50 * (evm + rhoevm * vv);
gp[__LI_GP(3,l)] = c;
gm[__LI_GM(3,l)] = evm-c;
evm = ww * rhov;
c = 0.50 * (evm + rhoevm * ww);
gp[__LI_GP(4,l)] = c;
gm[__LI_GM(4,l)] = evm-c;
evm = h * rhov;
c = 0.50 * (evm + evmax * (rho*h-pp));
gp[__LI_GP(5,l)] = c;
gm[__LI_GM(5,l)] = evm-c;
}
wenorec_ord = max(weno_scheme+ep_ord_change_gpu[__I4_EP_ORD_CHANGE(i,j,k,2)],1);
wenorec_1d_rusanov_euler_y_hybrid_rusanov_kernel_0(nv,weno_scheme,weno_version,wenorec_ord,rho0,u0,gp,gm,fj);
for(int m=1; m<5+1; m++){
fhat_gpu[__I4_FHAT(i,j,k,m)] = fj[__LI_FJ(m)];
}
}
}

}
}


extern "C"{
void euler_y_hybrid_rusanov_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int ng,int nv_aux,int eul_jmin,int eul_jmax,int lmax_base,int nkeep,int indx_cp_l,int indx_cp_r,int calorically_perfect,int force_zero_flux_min,int force_zero_flux_max,int weno_scheme,int weno_size,int weno_version,real sensor_threshold,real rgas0,real tol_iter_nr,real rho0,real u0,real t0,int *ep_ord_change_gpu,real *w_aux_gpu,real *fhat_gpu,real *cp_coeff_gpu,real *fl_gpu,real *coeff_deriv1_gpu,real *detady_gpu){
dim3 block(TWO_X,TWO_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((euler_y_hybrid_rusanov_kernel),grid,block,0,stream,nv,nx,ny,nz,ng,nv_aux,eul_jmin,eul_jmax,lmax_base,nkeep,indx_cp_l,indx_cp_r,calorically_perfect,force_zero_flux_min,force_zero_flux_max,weno_scheme,weno_size,weno_version,sensor_threshold,rgas0,tol_iter_nr,rho0,u0,t0,ep_ord_change_gpu,w_aux_gpu,fhat_gpu,cp_coeff_gpu,fl_gpu,coeff_deriv1_gpu,detady_gpu);
}
}

__device__ void eigenvectors_x_bc_nr_lat_x_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_x_bc_nr_lat_x_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + uu * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv );
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2*uu;
el[__LI_EL(3,2)] = b2*vv;
el[__LI_EL(4,2)] = b2*ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -vv;
el[__LI_EL(2,3)] = 0.0;
el[__LI_EL(3,3)] = 1.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - uu * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv );
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 0.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + c;
er[__LI_ER(1,3)] = vv;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 1.0;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - uu * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + uu * c;


}

__device__ real get_gamloc_dev_bc_nr_lat_x_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bc_nr_lat_x_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) bc_nr_lat_x_kernel(int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dcsidx_gpu,real *winf_gpu){
//Kernel for bc_nr_lat_x_kernel
real w_target;real df;real uu;
real vv;real ww;real h;
real qq;real cc;real c;
real ci;real b2;real b1;
real rho;real tt;real gamloc;
real p_rho;real p_e;real etot;
real b3;
int i;int j;int k;
int l;int m;int mm;
int sgn_dw;
real c_one[3];
#undef __LI_C_ONE
#define __LI_C_ONE(i) (i-(1))
real dw_dn[5];
#undef __LI_DW_DN
#define __LI_DW_DN(i) (i-(1))
real dwc_dn[5];
#undef __LI_DWC_DN
#define __LI_DWC_DN(i) (i-(1))
real ev[5];
#undef __LI_EV
#define __LI_EV(i) (i-(1))
real dw_dn_outer[5];
#undef __LI_DW_DN_OUTER
#define __LI_DW_DN_OUTER(i) (i-(1))
real dwc_dn_outer[5];
#undef __LI_DWC_DN_OUTER
#define __LI_DWC_DN_OUTER(i) (i-(1))
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
c_one[__LI_C_ONE(1)] = -1.50;
c_one[__LI_C_ONE(2)] =  2.0;
c_one[__LI_C_ONE(3)] =  -0.50;
if(start_or_end == 1) {
i = 1;
sgn_dw = 1;
}else if(start_or_end == 2) {
i = nx;
sgn_dw = -1;
}
for(int m=1; m<5+1; m++){
dw_dn[__LI_DW_DN(m)] = 0.0;
for(int l=1; l<3+1; l++){
dw_dn[__LI_DW_DN(m)] = dw_dn[__LI_DW_DN(m)] + sgn_dw * c_one[__LI_C_ONE(l)]*w_gpu[__I4_W(i+sgn_dw*(l-1),j,k,m)];
}
w_target = 0.0;
if (nr_type == 2) {
w_target = w_gpu[__I4_W(i-sgn_dw,j,k,m)];
}else if(nr_type == 3) {
w_target = winf_gpu[__I1_WINF(m)];
}
dw_dn_outer[__LI_DW_DN_OUTER(m)] = sgn_dw * (w_gpu[__I4_W(i,j,k,m)]-w_target);
}
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
gamloc = get_gamloc_dev_bc_nr_lat_x_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
cc = gamloc * tt * rgas0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*(gamloc-1.0);
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);
eigenvectors_x_bc_nr_lat_x_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = 0.0;
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
dwc_dn[__LI_DWC_DN(m)] = dwc_dn[__LI_DWC_DN(m)] + el[__LI_EL(mm,m)] * dw_dn[__LI_DW_DN(mm)];
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)] + el[__LI_EL(mm,m)] * dw_dn_outer[__LI_DW_DN_OUTER(mm)];
}
}
ev[__LI_EV(1)] = uu-c;
ev[__LI_EV(2)] = uu;
ev[__LI_EV(3)] = ev[__LI_EV(2)];
ev[__LI_EV(4)] = ev[__LI_EV(2)];
ev[__LI_EV(5)] = uu+c;
if (nr_type == 1) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = sgn_dw*min(sgn_dw*ev[__LI_EV(l)] ,0.0);
}
}
if (nr_type == 2 || nr_type == 3) {
for(int m=1; m<5+1; m++){
if(sgn_dw*ev[__LI_EV(m)] > 0.0) {
dwc_dn[__LI_DWC_DN(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)];
}
}
}
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = ev[__LI_EV(m)] * dwc_dn[__LI_DWC_DN(m)];
}
if (nr_type == 6) {
dwc_dn[__LI_DWC_DN(2)] = 0.0;
dwc_dn[__LI_DWC_DN(3)] = 0.0;
dwc_dn[__LI_DWC_DN(4)] = 0.0;
if(start_or_end == 1) {
dwc_dn[__LI_DWC_DN(5)] = dwc_dn[__LI_DWC_DN(1)];
}else if(start_or_end == 2) {
dwc_dn[__LI_DWC_DN(1)] = dwc_dn[__LI_DWC_DN(5)];
}
}
for(int m=1; m<5+1; m++){
df = 0.0;
for(int mm=1; mm<5+1; mm++){
df = df + er[__LI_ER(mm,m)] * dwc_dn[__LI_DWC_DN(mm)];
}
fl_gpu[__I4_FL(i,j,k,m)] = fl_gpu[__I4_FL(i,j,k,m)] + df * dcsidx_gpu[__I1_DCSIDX(i)];
}

}
}


extern "C"{
void bc_nr_lat_x_kernel_wrapper(hipStream_t stream,int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dcsidx_gpu,real *winf_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bc_nr_lat_x_kernel),grid,block,0,stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu);
}
}

__device__ void eigenvectors_x_c2_bc_nr_lat_x_c2_kernel_0(real dcsidxav,real dcsidyav,real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real ut,real *el,real *er){
//Device kernel for eigenvectors_x_c2_bc_nr_lat_x_c2_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ut * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + dcsidxav * ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + dcsidyav * ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2*uu;
el[__LI_EL(3,2)] = b2*vv;
el[__LI_EL(4,2)] = b2*ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = dcsidyav * uu - dcsidxav * vv;
el[__LI_EL(2,3)] = - dcsidyav;
el[__LI_EL(3,3)] = dcsidxav;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ut * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - dcsidxav * ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - dcsidyav * ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - dcsidxav * c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = - dcsidyav;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + dcsidxav * c;
er[__LI_ER(1,3)] = vv - dcsidyav * c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = dcsidxav;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + dcsidyav * c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - ut * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = - dcsidyav * uu + dcsidxav * vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + ut * c;


}

__device__ real get_gamloc_dev_bc_nr_lat_x_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bc_nr_lat_x_c2_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) bc_nr_lat_x_c2_kernel(int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dcsidx_gpu,real *winf_gpu,real *jac_gpu,real *mcsi_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *dcsidxc2_gpu,real *dcsidyc2_gpu){
//Kernel for bc_nr_lat_x_c2_kernel
real w_target;real df;real uu;
real vv;real ww;real h;
real qq;real cc;real c;
real ci;real b2;real b1;
real dcsixjdcsi;real dcsiyjdcsi;real pp;
real ut;real rho;real tt;
real gamloc;real p_rho;real p_e;
real etot;real b3;
int i;int j;int k;
int l;int m;int mm;
int sgn_dw;
real c_one[3];
#undef __LI_C_ONE
#define __LI_C_ONE(i) (i-(1))
real dw_dn[5];
#undef __LI_DW_DN
#define __LI_DW_DN(i) (i-(1))
real dwc_dn[5];
#undef __LI_DWC_DN
#define __LI_DWC_DN(i) (i-(1))
real ev[5];
#undef __LI_EV
#define __LI_EV(i) (i-(1))
real dw_dn_outer[5];
#undef __LI_DW_DN_OUTER
#define __LI_DW_DN_OUTER(i) (i-(1))
real dwc_dn_outer[5];
#undef __LI_DWC_DN_OUTER
#define __LI_DWC_DN_OUTER(i) (i-(1))
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))

j = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(j,ny,1)&&loop_cond(k,nz,1)){
c_one[__LI_C_ONE(1)] = -1.50;
c_one[__LI_C_ONE(2)] =  2.0;
c_one[__LI_C_ONE(3)] =  -0.50;
if(start_or_end == 1) {
i = 1;
sgn_dw = 1;
}else if(start_or_end == 2) {
i = nx;
sgn_dw = -1;
}
for(int m=1; m<nv+1; m++){
dw_dn[__LI_DW_DN(m)] = 0.0;
for(int l=1; l<3+1; l++){
dw_dn[__LI_DW_DN(m)] = dw_dn[__LI_DW_DN(m)] + sgn_dw * c_one[__LI_C_ONE(l)]*w_gpu[__I4_W(i+sgn_dw*(l-1),j,k,m)];
}
w_target = 0.0;
if (nr_type == 2) {
w_target = w_gpu[__I4_W(i-sgn_dw,j,k,m)];
}else if(nr_type == 3) {
w_target = winf_gpu[__I1_WINF(m)];
}
dw_dn_outer[__LI_DW_DN_OUTER(m)] = sgn_dw * (w_gpu[__I4_W(i,j,k,m)]-w_target);
}
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
gamloc = get_gamloc_dev_bc_nr_lat_x_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
cc = gamloc * tt * rgas0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*(gamloc-1.0);
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);
ut = dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)] * uu + dcsidync2_gpu[__I2_DCSIDYNC2(i,j)] * vv;
eigenvectors_x_c2_bc_nr_lat_x_c2_kernel_0(dcsidxnc2_gpu[__I2_DCSIDXNC2(i,j)],dcsidync2_gpu[__I2_DCSIDYNC2(i,j)],b1,b2,b3,uu,vv,ww,c,ci,h,ut,el,er);
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = 0.0;
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
dwc_dn[__LI_DWC_DN(m)] = dwc_dn[__LI_DWC_DN(m)] + el[__LI_EL(mm,m)] * dw_dn[__LI_DW_DN(mm)];
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)] + el[__LI_EL(mm,m)] * dw_dn_outer[__LI_DW_DN_OUTER(mm)];
}
}
ev[__LI_EV(1)] = ut-c;
ev[__LI_EV(2)] = ut;
ev[__LI_EV(3)] = ev[__LI_EV(2)];
ev[__LI_EV(4)] = ev[__LI_EV(2)];
ev[__LI_EV(5)] = ut+c;
if (nr_type == 1) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = sgn_dw*min(sgn_dw*ev[__LI_EV(l)] ,0.0);
}
}
if (nr_type == 2 || nr_type == 3) {
for(int m=1; m<5+1; m++){
if(sgn_dw*ev[__LI_EV(m)] > 0.0) {
dwc_dn[__LI_DWC_DN(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)];
}
}
}
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = ev[__LI_EV(m)] * dwc_dn[__LI_DWC_DN(m)];
}
if (nr_type == 6) {
dwc_dn[__LI_DWC_DN(2)] = 0.0;
dwc_dn[__LI_DWC_DN(3)] = 0.0;
dwc_dn[__LI_DWC_DN(4)] = 0.0;
if (start_or_end == 1) {
dwc_dn[__LI_DWC_DN(5)] = dwc_dn[__LI_DWC_DN(1)];
}else if(start_or_end == 2) {
dwc_dn[__LI_DWC_DN(1)] = dwc_dn[__LI_DWC_DN(5)];
}
}
for(int m=1; m<5+1; m++){
df = 0.0;
for(int mm=1; mm<5+1; mm++){
df = df + er[__LI_ER(mm,m)] * dwc_dn[__LI_DWC_DN(mm)];
}
fl_gpu[__I4_FL(i,j,k,m)] = fl_gpu[__I4_FL(i,j,k,m)] + df * mcsi_gpu[__I2_MCSI(i,j)];
}
dcsixjdcsi = sgn_dw * (c_one[__LI_C_ONE(1)]*dcsidxc2_gpu[__I2_DCSIDXC2(i,j)]/jac_gpu[__I2_JAC(i,j)]+c_one[__LI_C_ONE(2)]*dcsidxc2_gpu[__I2_DCSIDXC2(i+1*sgn_dw,j)]/jac_gpu[__I2_JAC(i+1*sgn_dw,j)]+ c_one[__LI_C_ONE(3)]*dcsidxc2_gpu[__I2_DCSIDXC2(i+2*sgn_dw,j)]/jac_gpu[__I2_JAC(i+2*sgn_dw,j)]);
dcsiyjdcsi = sgn_dw * (c_one[__LI_C_ONE(1)]*dcsidyc2_gpu[__I2_DCSIDYC2(i,j)]/jac_gpu[__I2_JAC(i,j)]+c_one[__LI_C_ONE(2)]*dcsidyc2_gpu[__I2_DCSIDYC2(i+1*sgn_dw,j)]/jac_gpu[__I2_JAC(i+1*sgn_dw,j)]+ c_one[__LI_C_ONE(3)]*dcsidyc2_gpu[__I2_DCSIDYC2(i+2*sgn_dw,j)]/jac_gpu[__I2_JAC(i+2*sgn_dw,j)]);
pp = p_rho * rho;
fl_gpu[__I4_FL(i,j,k,1)] = fl_gpu[__I4_FL(i,j,k,1)]+( w_gpu[__I4_W(i,j,k,1)]*uu *dcsixjdcsi+ w_gpu[__I4_W(i,j,k,1)]*vv *dcsiyjdcsi)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)]+((w_gpu[__I4_W(i,j,k,2)]*uu+pp)*dcsixjdcsi+ w_gpu[__I4_W(i,j,k,2)]*vv *dcsiyjdcsi)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)]+( w_gpu[__I4_W(i,j,k,3)]*uu *dcsixjdcsi+(w_gpu[__I4_W(i,j,k,3)]*vv+pp)*dcsiyjdcsi)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)]+( w_gpu[__I4_W(i,j,k,4)]*uu *dcsixjdcsi+ w_gpu[__I4_W(i,j,k,4)]*vv *dcsiyjdcsi)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)]+((w_gpu[__I4_W(i,j,k,5)]+pp)*uu*dcsixjdcsi+(w_gpu[__I4_W(i,j,k,5)]+pp)*vv*dcsiyjdcsi)*jac_gpu[__I2_JAC(i,j)];

}
}


extern "C"{
void bc_nr_lat_x_c2_kernel_wrapper(hipStream_t stream,int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dcsidx_gpu,real *winf_gpu,real *jac_gpu,real *mcsi_gpu,real *dcsidxnc2_gpu,real *dcsidync2_gpu,real *dcsidxc2_gpu,real *dcsidyc2_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((ny)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bc_nr_lat_x_c2_kernel),grid,block,0,stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dcsidx_gpu,winf_gpu,jac_gpu,mcsi_gpu,dcsidxnc2_gpu,dcsidync2_gpu,dcsidxc2_gpu,dcsidyc2_gpu);
}
}

__device__ void eigenvectors_y_bc_nr_lat_y_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_y_bc_nr_lat_y_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + vv * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu );
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -uu;
el[__LI_EL(2,3)] = 1.0;
el[__LI_EL(3,3)] = 0.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - vv * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu );
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 1.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu;
er[__LI_ER(1,3)] = vv - c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 0.0;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - vv * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = uu;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + vv * c;


}

__device__ real get_gamloc_dev_bc_nr_lat_y_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bc_nr_lat_y_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) bc_nr_lat_y_kernel(int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *detady_gpu,real *winf_gpu){
//Kernel for bc_nr_lat_y_kernel
real w_target;real df;real uu;
real vv;real ww;real h;
real qq;real cc;real c;
real ci;real b2;real b1;
real rho;real tt;real gamloc;
real p_rho;real p_e;real etot;
real b3;
int i;int j;int k;
int l;int m;int mm;
int sgn_dw;
real c_one[3];
#undef __LI_C_ONE
#define __LI_C_ONE(i) (i-(1))
real dw_dn[5];
#undef __LI_DW_DN
#define __LI_DW_DN(i) (i-(1))
real dwc_dn[5];
#undef __LI_DWC_DN
#define __LI_DWC_DN(i) (i-(1))
real ev[5];
#undef __LI_EV
#define __LI_EV(i) (i-(1))
real dw_dn_outer[5];
#undef __LI_DW_DN_OUTER
#define __LI_DW_DN_OUTER(i) (i-(1))
real dwc_dn_outer[5];
#undef __LI_DWC_DN_OUTER
#define __LI_DWC_DN_OUTER(i) (i-(1))
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
c_one[__LI_C_ONE(1)] = -1.50;
c_one[__LI_C_ONE(2)] =  2.0;
c_one[__LI_C_ONE(3)] =  -0.50;
if(start_or_end == 1) {
j = 1;
sgn_dw = 1;
}else if(start_or_end == 2) {
j = ny;
sgn_dw = -1;
}
for(int m=1; m<5+1; m++){
dw_dn[__LI_DW_DN(m)] = 0.0;
for(int l=1; l<3+1; l++){
dw_dn[__LI_DW_DN(m)] = dw_dn[__LI_DW_DN(m)] + sgn_dw * c_one[__LI_C_ONE(l)]*w_gpu[__I4_W(i,j+sgn_dw*(l-1),k,m)];
}
w_target = 0.0;
if (nr_type == 2) {
w_target = w_gpu[__I4_W(i,j-sgn_dw,k,m)];
}else if(nr_type == 3) {
w_target = winf_gpu[__I1_WINF(m)];
}
dw_dn_outer[__LI_DW_DN_OUTER(m)] = sgn_dw * (w_gpu[__I4_W(i,j,k,m)]-w_target);
}
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
gamloc = get_gamloc_dev_bc_nr_lat_y_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
cc = gamloc * tt * rgas0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*(gamloc-1.0);
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);
eigenvectors_y_bc_nr_lat_y_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = 0.0;
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
dwc_dn[__LI_DWC_DN(m)] = dwc_dn[__LI_DWC_DN(m)] + el[__LI_EL(mm,m)] * dw_dn[__LI_DW_DN(mm)];
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)] + el[__LI_EL(mm,m)] * dw_dn_outer[__LI_DW_DN_OUTER(mm)];
}
}
ev[__LI_EV(1)] = vv-c;
ev[__LI_EV(2)] = vv;
ev[__LI_EV(3)] = ev[__LI_EV(2)];
ev[__LI_EV(4)] = ev[__LI_EV(2)];
ev[__LI_EV(5)] = vv+c;
if (nr_type == 1) {
if(start_or_end == 1) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = min(ev[__LI_EV(l)] ,0.0);
}
}else if(start_or_end == 2) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = max(ev[__LI_EV(l)] ,0.0);
}
}
}
if (nr_type == 2 || nr_type == 3) {
for(int m=1; m<5+1; m++){
if(sgn_dw*ev[__LI_EV(m)] > 0.0) {
dwc_dn[__LI_DWC_DN(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)];
}
}
}
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = ev[__LI_EV(m)] * dwc_dn[__LI_DWC_DN(m)];
}
if (nr_type == 6) {
dwc_dn[__LI_DWC_DN(2)] = 0.0;
dwc_dn[__LI_DWC_DN(3)] = 0.0;
dwc_dn[__LI_DWC_DN(4)] = 0.0;
if(start_or_end == 1) {
dwc_dn[__LI_DWC_DN(5)] = dwc_dn[__LI_DWC_DN(1)];
}else if(start_or_end == 2) {
dwc_dn[__LI_DWC_DN(1)] = dwc_dn[__LI_DWC_DN(5)];
}
}
for(int m=1; m<5+1; m++){
df = 0.0;
for(int mm=1; mm<5+1; mm++){
df = df + er[__LI_ER(mm,m)] * dwc_dn[__LI_DWC_DN(mm)];
}
fl_gpu[__I4_FL(i,j,k,m)] = fl_gpu[__I4_FL(i,j,k,m)] + df * detady_gpu[__I1_DETADY(j)];
}

}
}


extern "C"{
void bc_nr_lat_y_kernel_wrapper(hipStream_t stream,int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *detady_gpu,real *winf_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bc_nr_lat_y_kernel),grid,block,0,stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,winf_gpu);
}
}

__device__ void eigenvectors_y_c2_bc_nr_lat_y_c2_kernel_0(real detadxav,real detadyav,real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real ut,real *el,real *er){
//Device kernel for eigenvectors_y_c2_bc_nr_lat_y_c2_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ut * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu + detadxav * ci);
el[__LI_EL(3,1)] = -0.50 * (b2 * vv + detadyav * ci);
el[__LI_EL(4,1)] = -0.50 * (b2 * ww );
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -detadyav * uu + detadxav * vv;
el[__LI_EL(2,3)] = +detadyav;
el[__LI_EL(3,3)] = -detadxav;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -ww;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 0.0;
el[__LI_EL(4,4)] = 1.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ut * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu - detadxav * ci);
el[__LI_EL(3,5)] = -0.50 * (b2 * vv - detadyav * ci);
el[__LI_EL(4,5)] = -0.50 * (b2 * ww );
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu - detadxav * c;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = detadyav;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu + detadxav * c;
er[__LI_ER(1,3)] = vv - detadyav * c;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = - detadxav;
er[__LI_ER(4,3)] = 0.0;
er[__LI_ER(5,3)] = vv + detadyav * c;
er[__LI_ER(1,4)] = ww;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 1.0;
er[__LI_ER(5,4)] = ww;
er[__LI_ER(1,5)] = h - ut * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = detadyav * uu - detadxav * vv;
er[__LI_ER(4,5)] = ww;
er[__LI_ER(5,5)] = h + ut * c;


}

__device__ real get_gamloc_dev_bc_nr_lat_y_c2_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bc_nr_lat_y_c2_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) bc_nr_lat_y_c2_kernel(int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,int *wall_tag_gpu,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *detady_gpu,real *wfar_gpu,real *jac_gpu,real *meta_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *detadxc2_gpu,real *detadyc2_gpu){
//Kernel for bc_nr_lat_y_c2_kernel
real w_target;real df;real uu;
real vv;real ww;real h;
real qq;real cc;real c;
real ci;real b2;real b1;
real rho;real tt;real gamloc;
real p_rho;real p_e;real etot;
real b3;real detaxjdeta;real detayjdeta;
real pp;real ut;
int i;int j;int k;
int l;int m;int mm;
int sgn_dw;
real c_one[3];
#undef __LI_C_ONE
#define __LI_C_ONE(i) (i-(1))
real dw_dn[5];
#undef __LI_DW_DN
#define __LI_DW_DN(i) (i-(1))
real dwc_dn[5];
#undef __LI_DWC_DN
#define __LI_DWC_DN(i) (i-(1))
real ev[5];
#undef __LI_EV
#define __LI_EV(i) (i-(1))
real dw_dn_outer[5];
#undef __LI_DW_DN_OUTER
#define __LI_DW_DN_OUTER(i) (i-(1))
real dwc_dn_outer[5];
#undef __LI_DWC_DN_OUTER
#define __LI_DWC_DN_OUTER(i) (i-(1))
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
c_one[__LI_C_ONE(1)] = -1.50;
c_one[__LI_C_ONE(2)] =  2.0;
c_one[__LI_C_ONE(3)] =  -0.50;
if(start_or_end == 1) {
j = 1;
sgn_dw = 1;
}else if(start_or_end == 2) {
j = ny;
sgn_dw = -1;
}
if(((start_or_end == 1) && (wall_tag_gpu[__I1_WALL_TAG(i)] > 0)) ) {
sgn_dw = 1;
}else {
for(int m=1; m<nv+1; m++){
dw_dn[__LI_DW_DN(m)] = 0.0;
for(int l=1; l<3+1; l++){
dw_dn[__LI_DW_DN(m)] = dw_dn[__LI_DW_DN(m)] + sgn_dw * c_one[__LI_C_ONE(l)]*w_gpu[__I4_W(i,j+sgn_dw*(l-1),k,m)];
}
w_target = 0.0;
if (nr_type == 2) {
w_target = w_gpu[__I4_W(i,j-sgn_dw,k,m)];
}else if(nr_type == 3) {
w_target = wfar_gpu[__I2_WFAR(i,m)];
}
dw_dn_outer[__LI_DW_DN_OUTER(m)] = sgn_dw * (w_gpu[__I4_W(i,j,k,m)]-w_target);
}
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
gamloc = get_gamloc_dev_bc_nr_lat_y_c2_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
cc = gamloc * tt * rgas0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*(gamloc-1.0);
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);
ut = detadxnc2_gpu[__I2_DETADXNC2(i,j)] * uu + detadync2_gpu[__I2_DETADYNC2(i,j)] * vv;
eigenvectors_y_c2_bc_nr_lat_y_c2_kernel_0(detadxnc2_gpu[__I2_DETADXNC2(i,j)],detadync2_gpu[__I2_DETADYNC2(i,j)],b1,b2,b3,uu,vv,ww,c,ci,h,ut,el,er);
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = 0.0;
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
dwc_dn[__LI_DWC_DN(m)] = dwc_dn[__LI_DWC_DN(m)] + el[__LI_EL(mm,m)] * dw_dn[__LI_DW_DN(mm)];
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)] + el[__LI_EL(mm,m)] * dw_dn_outer[__LI_DW_DN_OUTER(mm)];
}
}
ev[__LI_EV(1)] = ut-c;
ev[__LI_EV(2)] = ut;
ev[__LI_EV(3)] = ev[__LI_EV(2)];
ev[__LI_EV(4)] = ev[__LI_EV(2)];
ev[__LI_EV(5)] = ut+c;
if(nr_type == 1) {
if(start_or_end == 1) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = min(ev[__LI_EV(l)] ,0.0);
}
}else if(start_or_end == 2) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = max(ev[__LI_EV(l)] ,0.0);
}
}
}
if (nr_type == 2 || nr_type == 3) {
for(int m=1; m<5+1; m++){
if(sgn_dw*ev[__LI_EV(m)] > 0.0) {
dwc_dn[__LI_DWC_DN(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)];
}
}
}
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = ev[__LI_EV(m)] * dwc_dn[__LI_DWC_DN(m)];
}
if (nr_type == 6) {
dwc_dn[__LI_DWC_DN(2)] = 0.0;
dwc_dn[__LI_DWC_DN(3)] = 0.0;
dwc_dn[__LI_DWC_DN(4)] = 0.0;
if(start_or_end == 1) {
dwc_dn[__LI_DWC_DN(5)] = dwc_dn[__LI_DWC_DN(1)];
}else if(start_or_end == 2) {
dwc_dn[__LI_DWC_DN(1)] = dwc_dn[__LI_DWC_DN(5)];
}
}
for(int m=1; m<5+1; m++){
df = 0.0;
for(int mm=1; mm<5+1; mm++){
df = df + er[__LI_ER(mm,m)] * dwc_dn[__LI_DWC_DN(mm)];
}
fl_gpu[__I4_FL(i,j,k,m)] = fl_gpu[__I4_FL(i,j,k,m)] + df * meta_gpu[__I2_META(i,j)];
}
detaxjdeta=sgn_dw*(c_one[__LI_C_ONE(1)]*detadxc2_gpu[__I2_DETADXC2(i,j)]/jac_gpu[__I2_JAC(i,j)]+c_one[__LI_C_ONE(2)]*detadxc2_gpu[__I2_DETADXC2(i,j+1*sgn_dw)]/jac_gpu[__I2_JAC(i,j+1*sgn_dw)]+ c_one[__LI_C_ONE(3)]*detadxc2_gpu[__I2_DETADXC2(i,j+2*sgn_dw)]/jac_gpu[__I2_JAC(i,j+2*sgn_dw)]);
detayjdeta=sgn_dw*(c_one[__LI_C_ONE(1)]*detadyc2_gpu[__I2_DETADYC2(i,j)]/jac_gpu[__I2_JAC(i,j)]+c_one[__LI_C_ONE(2)]*detadyc2_gpu[__I2_DETADYC2(i,j+1*sgn_dw)]/jac_gpu[__I2_JAC(i,j+1*sgn_dw)]+ c_one[__LI_C_ONE(3)]*detadyc2_gpu[__I2_DETADYC2(i,j+2*sgn_dw)]/jac_gpu[__I2_JAC(i,j+2*sgn_dw)]);
pp = p_rho * rho;
fl_gpu[__I4_FL(i,j,k,1)] = fl_gpu[__I4_FL(i,j,k,1)] +( w_gpu[__I4_W(i,j,k,1)]*uu *detaxjdeta + w_gpu[__I4_W(i,j,k,1)]*vv *detayjdeta)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,2)] = fl_gpu[__I4_FL(i,j,k,2)] +((w_gpu[__I4_W(i,j,k,2)]*uu+pp)*detaxjdeta + w_gpu[__I4_W(i,j,k,2)]*vv *detayjdeta)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,3)] = fl_gpu[__I4_FL(i,j,k,3)] +( w_gpu[__I4_W(i,j,k,3)]*uu *detaxjdeta + (w_gpu[__I4_W(i,j,k,3)]*vv+pp)*detayjdeta)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,4)] = fl_gpu[__I4_FL(i,j,k,4)] +( w_gpu[__I4_W(i,j,k,4)]*uu *detaxjdeta + w_gpu[__I4_W(i,j,k,4)]*vv *detayjdeta)*jac_gpu[__I2_JAC(i,j)];
fl_gpu[__I4_FL(i,j,k,5)] = fl_gpu[__I4_FL(i,j,k,5)] +((w_gpu[__I4_W(i,j,k,5)]+pp)*uu*detaxjdeta + (w_gpu[__I4_W(i,j,k,5)]+pp)*vv*detayjdeta)*jac_gpu[__I2_JAC(i,j)];
}

}
}


extern "C"{
void bc_nr_lat_y_c2_kernel_wrapper(hipStream_t stream,int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,int *wall_tag_gpu,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *detady_gpu,real *wfar_gpu,real *jac_gpu,real *meta_gpu,real *detadxnc2_gpu,real *detadync2_gpu,real *detadxc2_gpu,real *detadyc2_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((bc_nr_lat_y_c2_kernel),grid,block,0,stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,wall_tag_gpu,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,detady_gpu,wfar_gpu,jac_gpu,meta_gpu,detadxnc2_gpu,detadync2_gpu,detadxc2_gpu,detadyc2_gpu);
}
}

__device__ void eigenvectors_z_bc_nr_lat_z_kernel_0(real b1,real b2,real b3,real uu,real vv,real ww,real c,real ci,real h,real *el,real *er){
//Device kernel for eigenvectors_z_bc_nr_lat_z_kernel_0
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))


el[__LI_EL(1,1)] = 0.50 * (b1 + ww * ci);
el[__LI_EL(2,1)] = -0.50 * (b2 * uu );
el[__LI_EL(3,1)] = -0.50 * (b2 * vv );
el[__LI_EL(4,1)] = -0.50 * (b2 * ww + ci);
el[__LI_EL(5,1)] = 0.50 * b2;
el[__LI_EL(1,2)] = 1.0 - b1;
el[__LI_EL(2,2)] = b2 * uu;
el[__LI_EL(3,2)] = b2 * vv;
el[__LI_EL(4,2)] = b2 * ww;
el[__LI_EL(5,2)] = -b2;
el[__LI_EL(1,3)] = -uu;
el[__LI_EL(2,3)] = 1.0;
el[__LI_EL(3,3)] = 0.0;
el[__LI_EL(4,3)] = 0.0;
el[__LI_EL(5,3)] = 0.0;
el[__LI_EL(1,4)] = -vv;
el[__LI_EL(2,4)] = 0.0;
el[__LI_EL(3,4)] = 1.0;
el[__LI_EL(4,4)] = 0.0;
el[__LI_EL(5,4)] = 0.0;
el[__LI_EL(1,5)] = 0.50 * (b1 - ww * ci);
el[__LI_EL(2,5)] = -0.50 * (b2 * uu );
el[__LI_EL(3,5)] = -0.50 * (b2 * vv );
el[__LI_EL(4,5)] = -0.50 * (b2 * ww - ci);
el[__LI_EL(5,5)] = 0.50 * b2;
er[__LI_ER(1,1)] = 1.0;
er[__LI_ER(2,1)] = 1.0;
er[__LI_ER(3,1)] = 0.0;
er[__LI_ER(4,1)] = 0.0;
er[__LI_ER(5,1)] = 1.0;
er[__LI_ER(1,2)] = uu;
er[__LI_ER(2,2)] = uu;
er[__LI_ER(3,2)] = 1.0;
er[__LI_ER(4,2)] = 0.0;
er[__LI_ER(5,2)] = uu;
er[__LI_ER(1,3)] = vv;
er[__LI_ER(2,3)] = vv;
er[__LI_ER(3,3)] = 0.0;
er[__LI_ER(4,3)] = 1.0;
er[__LI_ER(5,3)] = vv;
er[__LI_ER(1,4)] = ww - c;
er[__LI_ER(2,4)] = ww;
er[__LI_ER(3,4)] = 0.0;
er[__LI_ER(4,4)] = 0.0;
er[__LI_ER(5,4)] = ww + c;
er[__LI_ER(1,5)] = h - ww * c;
er[__LI_ER(2,5)] = b3;
er[__LI_ER(3,5)] = uu;
er[__LI_ER(4,5)] = vv;
er[__LI_ER(5,5)] = h + ww * c;


}

__device__ real get_gamloc_dev_bc_nr_lat_z_kernel_0(int indx_cp_l,int indx_cp_r,int calorically_perfect,real t0,real rgas0,real tt,real *cp_coeff_gpu){
//Device kernel for get_gamloc_dev_bc_nr_lat_z_kernel_0
real get_gamloc_dev;real cploc;real gamloc;
int l;


if (calorically_perfect==1) {
cploc = cp_coeff_gpu[__I1_CP_COEFF(0)];
}else {
cploc = 0.0;
for(int l=indx_cp_l; l<indx_cp_r+1; l++){
cploc=cploc+cp_coeff_gpu[__I1_CP_COEFF(l)]*pow((tt/t0),l);
}
}
gamloc = cploc/(cploc-rgas0);
get_gamloc_dev = gamloc;


return get_gamloc_dev;
}




__global__ void __launch_bounds__(256) bc_nr_lat_z_kernel(int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dzitdz_gpu,real *winf_gpu){
//Kernel for bc_nr_lat_z_kernel
real w_target;real df;real uu;
real vv;real ww;real h;
real qq;real cc;real c;
real ci;real b2;real b1;
real rho;real tt;real gamloc;
real p_rho;real p_e;real etot;
real b3;
int i;int j;int k;
int l;int m;int mm;
int sgn_dw;
real c_one[3];
#undef __LI_C_ONE
#define __LI_C_ONE(i) (i-(1))
real dw_dn[5];
#undef __LI_DW_DN
#define __LI_DW_DN(i) (i-(1))
real dwc_dn[5];
#undef __LI_DWC_DN
#define __LI_DWC_DN(i) (i-(1))
real ev[5];
#undef __LI_EV
#define __LI_EV(i) (i-(1))
real dw_dn_outer[5];
#undef __LI_DW_DN_OUTER
#define __LI_DW_DN_OUTER(i) (i-(1))
real dwc_dn_outer[5];
#undef __LI_DWC_DN_OUTER
#define __LI_DWC_DN_OUTER(i) (i-(1))
real el[5*5];
#undef __LI_EL
#define __LI_EL(i,j) ((i-(1))+5*(j-(1)))
real er[5*5];
#undef __LI_ER
#define __LI_ER(i,j) ((i-(1))+5*(j-(1)))

i = __GIDX(x,1);
j = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(j,ny,1)){
c_one[__LI_C_ONE(1)] = -1.50;
c_one[__LI_C_ONE(2)] =  2.0;
c_one[__LI_C_ONE(3)] =  -0.50;
if(start_or_end == 1) {
k = 1;
sgn_dw = 1;
}else if(start_or_end == 2) {
k = nz;
sgn_dw = -1;
}
for(int m=1; m<5+1; m++){
dw_dn[__LI_DW_DN(m)] = 0.0;
for(int l=1; l<3+1; l++){
dw_dn[__LI_DW_DN(m)] = dw_dn[__LI_DW_DN(m)] + sgn_dw * c_one[__LI_C_ONE(l)]*w_gpu[__I4_W(i,j,k+sgn_dw*(l-1),m)];
}
if (nr_type == 2) {
w_target = w_gpu[__I4_W(i,j,k-sgn_dw,m)];
}else if(nr_type == 3) {
w_target = winf_gpu[__I1_WINF(m)];
}
dw_dn_outer[__LI_DW_DN_OUTER(m)] = sgn_dw * (w_gpu[__I4_W(i,j,k,m)]-w_target);
}
rho = w_aux_gpu[__I4_W_AUX(i,j,k,1)];
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
h = w_aux_gpu[__I4_W_AUX(i,j,k,5)];
tt = w_aux_gpu[__I4_W_AUX(i,j,k,6)];
qq = 0.50 * (uu*uu +vv*vv + ww*ww);
gamloc = get_gamloc_dev_bc_nr_lat_z_kernel_0(indx_cp_l,indx_cp_r,calorically_perfect,t0,rgas0,tt,cp_coeff_gpu);
cc = gamloc * tt * rgas0;
c = sqrt(cc);
ci = 1.0/c;
p_rho = tt*rgas0;
p_e = rho*(gamloc-1.0);
etot = h - tt*rgas0;
b3 = etot - rho * p_rho/p_e;
b2 = p_e/(rho*cc);
b1 = p_rho/cc - b2*(etot - 2.0*qq);
eigenvectors_z_bc_nr_lat_z_kernel_0(b1,b2,b3,uu,vv,ww,c,ci,h,el,er);
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = 0.0;
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = 0.0;
for(int mm=1; mm<5+1; mm++){
dwc_dn[__LI_DWC_DN(m)] = dwc_dn[__LI_DWC_DN(m)] + el[__LI_EL(mm,m)] * dw_dn[__LI_DW_DN(mm)];
dwc_dn_outer[__LI_DWC_DN_OUTER(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)] + el[__LI_EL(mm,m)] * dw_dn_outer[__LI_DW_DN_OUTER(mm)];
}
}
ev[__LI_EV(1)] = ww-c;
ev[__LI_EV(2)] = ww;
ev[__LI_EV(3)] = ev[__LI_EV(2)];
ev[__LI_EV(4)] = ev[__LI_EV(2)];
ev[__LI_EV(5)] = ww+c;
if (nr_type == 1) {
if(start_or_end == 1) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = min(ev[__LI_EV(l)] ,0.0);
}
}else if(start_or_end == 2) {
for(int l=1; l<5+1; l++){
ev[__LI_EV(l)] = max(ev[__LI_EV(l)] ,0.0);
}
}
}
if (nr_type == 2 || nr_type == 3) {
for(int m=1; m<5+1; m++){
if(sgn_dw*ev[__LI_EV(m)] > 0.0) {
dwc_dn[__LI_DWC_DN(m)] = dwc_dn_outer[__LI_DWC_DN_OUTER(m)];
}
}
}
for(int m=1; m<5+1; m++){
dwc_dn[__LI_DWC_DN(m)] = ev[__LI_EV(m)] * dwc_dn[__LI_DWC_DN(m)];
}
if (nr_type == 6) {
dwc_dn[__LI_DWC_DN(2)] = 0.0;
dwc_dn[__LI_DWC_DN(3)] = 0.0;
dwc_dn[__LI_DWC_DN(4)] = 0.0;
if(start_or_end == 1) {
dwc_dn[__LI_DWC_DN(5)] = dwc_dn[__LI_DWC_DN(1)];
}else if(start_or_end == 2) {
dwc_dn[__LI_DWC_DN(1)] = dwc_dn[__LI_DWC_DN(5)];
}
}
for(int m=1; m<5+1; m++){
df = 0.0;
for(int mm=1; mm<5+1; mm++){
df = df + er[__LI_ER(mm,m)] * dwc_dn[__LI_DWC_DN(mm)];
}
fl_gpu[__I4_FL(i,j,k,m)] = fl_gpu[__I4_FL(i,j,k,m)] + df * dzitdz_gpu[__I1_DZITDZ(k)];
}

}
}


extern "C"{
void bc_nr_lat_z_kernel_wrapper(hipStream_t stream,int start_or_end,int nr_type,int nx,int ny,int nz,int ng,int nv,int indx_cp_l,int indx_cp_r,int calorically_perfect,real rgas0,real t0,real *cp_coeff_gpu,real *fl_gpu,real *w_aux_gpu,real *w_gpu,real *dzitdz_gpu,real *winf_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((ny)-(1)+1,block.y));

hipLaunchKernelGGL((bc_nr_lat_z_kernel),grid,block,0,stream,start_or_end,nr_type,nx,ny,nz,ng,nv,indx_cp_l,indx_cp_r,calorically_perfect,rgas0,t0,cp_coeff_gpu,fl_gpu,w_aux_gpu,w_gpu,dzitdz_gpu,winf_gpu);
}
}

__device__ void eigs33_insitu_swirling_kernel_0(real *rmat,real *rex,real *rimx){
//Device kernel for eigs33_insitu_swirling_kernel_0
real ddbl;real pdbl;real qdbl;
real otrd;real pi;real a;
real temps;real tempw;real b;
real somma;real c;real p;
real q;real sqdel;real sqp;
real teta;
int ii;int jj;int k;
#undef __LI_REX
#define __LI_REX(i) (i-(1))
#undef __LI_RIMX
#define __LI_RIMX(i) (i-(1))
real rey[3];
#undef __LI_REY
#define __LI_REY(i) (i-(1))
real rimy[3];
#undef __LI_RIMY
#define __LI_RIMY(i) (i-(1))
real reu[3];
#undef __LI_REU
#define __LI_REU(i) (i-(1))
real rimu[3];
#undef __LI_RIMU
#define __LI_RIMU(i) (i-(1))
real rev[3];
#undef __LI_REV
#define __LI_REV(i) (i-(1))
real rimv[3];
#undef __LI_RIMV
#define __LI_RIMV(i) (i-(1))
#undef __LI_RMAT
#define __LI_RMAT(i,j) ((i-(1))+3*(j-(1)))
real s[3*3];
#undef __LI_S
#define __LI_S(i,j) ((i-(1))+3*(j-(1)))
real ww[3*3];
#undef __LI_WW
#define __LI_WW(i,j) ((i-(1))+3*(j-(1)))
real temp[3*3];
#undef __LI_TEMP
#define __LI_TEMP(i,j) ((i-(1))+3*(j-(1)))


otrd = 1./3.0;
pi = acos(-1.0);
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
s[__LI_S(ii,jj)] = 0.5*(rmat[__LI_RMAT(ii,jj)]+rmat[__LI_RMAT(jj,ii)]);
ww[__LI_WW(ii,jj)] = 0.5*(rmat[__LI_RMAT(ii,jj)]-rmat[__LI_RMAT(jj,ii)]);
}
}
a = -(s[__LI_S(1,1)]+s[__LI_S(2,2)]+s[__LI_S(3,3)]);
temps = 0.;
tempw = 0.;
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
temps = temps + s[__LI_S(ii,jj)]*s[__LI_S(ii,jj)];
tempw = tempw + ww[__LI_WW(ii,jj)]*ww[__LI_WW(ii,jj)];
}
}
b=0.5*(((a)*(a))-temps+tempw);
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
temp[__LI_TEMP(ii,jj)]=0.;
for(int k=1; k<3+1; k++){
temp[__LI_TEMP(ii,jj)] = temp[__LI_TEMP(ii,jj)]+s[__LI_S(ii,k)]*s[__LI_S(k,jj)]+3.*ww[__LI_WW(ii,k)]*ww[__LI_WW(k,jj)];
}
}
}
somma=0.;
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
somma= somma+temp[__LI_TEMP(ii,jj)]*s[__LI_S(jj,ii)];
}
}
c=-1./3.*(pow(a,3)-3.*a*b+somma);
pdbl=b-((a)*(a))/3.;
qdbl=c-a*b/3.+2*pow(a,3)/27.;
ddbl=(((qdbl)*(qdbl)))/4.E0+(pow(pdbl,3))/27.E0;
p = pdbl;
q = qdbl;
if(ddbl>0.E0) {
sqdel = sqrt(ddbl);
reu[__LI_REU(1)] =-0.5*q+sqdel;
rev[__LI_REV(1)] =-0.5*q-sqdel;
reu[__LI_REU(1)]=SIGN(1.0,reu[__LI_REU(1)])*pow((abs(reu[__LI_REU(1)])),otrd);
rev[__LI_REV(1)]=SIGN(1.0,rev[__LI_REV(1)])*pow((abs(rev[__LI_REV(1)])),otrd);
reu[__LI_REU(2)] = -0.5*reu[__LI_REU(1)];
rev[__LI_REV(2)] = -0.5*rev[__LI_REV(1)];
reu[__LI_REU(3)] = reu[__LI_REU(2)];
rev[__LI_REV(3)] = rev[__LI_REV(2)];
rimu[__LI_RIMU(1)] = 0.;
rimv[__LI_RIMV(1)] = 0.;
rimu[__LI_RIMU(2)] = sqrt(3.0)/2.*reu[__LI_REU(1)];
rimv[__LI_RIMV(2)] = sqrt(3.0)/2.*rev[__LI_REV(1)];
rimu[__LI_RIMU(3)] = -rimu[__LI_RIMU(2)];
rimv[__LI_RIMV(3)] = -rimv[__LI_RIMV(2)];
rey[__LI_REY(1)] = reu[__LI_REU(1)]+rev[__LI_REV(1)];
rimy[__LI_RIMY(1)] = rimu[__LI_RIMU(1)]+rimv[__LI_RIMV(1)];
rey[__LI_REY(2)] = reu[__LI_REU(2)]+rev[__LI_REV(3)];
rimy[__LI_RIMY(2)] = rimu[__LI_RIMU(2)]+rimv[__LI_RIMV(3)];
rey[__LI_REY(3)] = reu[__LI_REU(3)]+rev[__LI_REV(2)];
rimy[__LI_RIMY(3)] = rimu[__LI_RIMU(3)]+rimv[__LI_RIMV(2)];
}else {
if (q==0.) {
rey[__LI_REY(1)] = 0.;
rey[__LI_REY(2)] = sqrt(-p);
rey[__LI_REY(3)] = -rey[__LI_REY(2)];
}else {
sqp = 2.*sqrt(-p/3.);
sqdel = sqrt(-ddbl);
if (q<0.) {
teta = atan(-2.*sqdel/q);
}else {
teta = pi+atan(-2.*sqdel/q);
}
rey[__LI_REY(1)] = sqp*cos(teta/3.);
rey[__LI_REY(2)] = sqp*cos((teta+2*pi)/3.);
rey[__LI_REY(3)] = sqp*cos((teta+4*pi)/3.);
}
rimy[__LI_RIMY(1)] = 0.;
rimy[__LI_RIMY(2)] = 0.;
rimy[__LI_RIMY(3)] = 0.;
}
rex[__LI_REX(1)] = rey[__LI_REY(1)]-(a/3.);
rimx[__LI_RIMX(1)] = rimy[__LI_RIMY(1)];
rex[__LI_REX(2)] = rey[__LI_REY(2)]-(a/3.);
rimx[__LI_RIMX(2)] = rimy[__LI_RIMY(2)];
rex[__LI_REX(3)] = rey[__LI_REY(3)]-(a/3.);
rimx[__LI_RIMX(3)] = rimy[__LI_RIMY(3)];
if (rimy[__LI_RIMY(2)]<0.) {
rex[__LI_REX(2)] = rey[__LI_REY(3)]-(a/3.);
rimx[__LI_RIMX(2)] = rimy[__LI_RIMY(3)];
rex[__LI_REX(3)] = rey[__LI_REY(2)]-(a/3.);
rimx[__LI_RIMX(3)] = rimy[__LI_RIMY(2)];
}


}




__global__ void __launch_bounds__(256) insitu_swirling_kernel(int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
//Kernel for insitu_swirling_kernel
real uu;real vv;real ww;
real ux;real uy;real uz;
real vx;real vy;real vz;
real wx;real wy;real wz;
real ccl;real div3l;real epsi2;
real omx;real omy;real omz;
real omod2;real div2;real div;
int i;int j;int k;
int l;
real eigr_a[3];
#undef __LI_EIGR_A
#define __LI_EIGR_A(i) (i-(1))
real eigi_a[3];
#undef __LI_EIGI_A
#define __LI_EIGI_A(i) (i-(1))
real astar[3*3];
#undef __LI_ASTAR
#define __LI_ASTAR(i,j) ((i-(1))+3*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=1; j<ny+1; j++){
uu = w_aux_gpu[__I4_W_AUX(i,j,k,2)];
vv = w_aux_gpu[__I4_W_AUX(i,j,k,3)];
ww = w_aux_gpu[__I4_W_AUX(i,j,k,4)];
ux = 0.0;
vx = 0.0;
wx = 0.0;
uy = 0.0;
vy = 0.0;
wy = 0.0;
uz = 0.0;
vz = 0.0;
wz = 0.0;
for(int l=1; l<visc_order/2+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,visc_order/2)];
ux = ux+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vx = vx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wx = wx+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
uy = uy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
vy = vy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
wy = wy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
uz = uz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vz = vz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wz = wz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
}
ux = ux*dcsidx_gpu[__I1_DCSIDX(i)];
vx = vx*dcsidx_gpu[__I1_DCSIDX(i)];
wx = wx*dcsidx_gpu[__I1_DCSIDX(i)];
uy = uy*detady_gpu[__I1_DETADY(j)];
vy = vy*detady_gpu[__I1_DETADY(j)];
wy = wy*detady_gpu[__I1_DETADY(j)];
uz = uz*dzitdz_gpu[__I1_DZITDZ(k)];
vz = vz*dzitdz_gpu[__I1_DZITDZ(k)];
wz = wz*dzitdz_gpu[__I1_DZITDZ(k)];
div = ux+vy+wz;
div3l = div/3.0;
epsi2=((u0)*(u0));
omz = vx-uy;
omx = wy-vz;
omy = uz-wx;
omod2 = omx*omx+omy*omy+omz*omz;
div2 = div*div;
astar[__LI_ASTAR(1,1)] = ux-div3l;
astar[__LI_ASTAR(1,2)] = uy;
astar[__LI_ASTAR(1,3)] = uz;
astar[__LI_ASTAR(2,1)] = vx;
astar[__LI_ASTAR(2,2)] = vy-div3l;
astar[__LI_ASTAR(2,3)] = vz;
astar[__LI_ASTAR(3,1)] = wx;
astar[__LI_ASTAR(3,2)] = wy;
astar[__LI_ASTAR(3,3)] = wz-div3l;
eigs33_insitu_swirling_kernel_0(astar,eigr_a,eigi_a);
psi_gpu[__I4_PSI(i,j,k,mpsi)] = 2.0*max(0.0,eigi_a[__LI_EIGI_A(2)]);
}

}
}


extern "C"{
void insitu_swirling_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((insitu_swirling_kernel),grid,block,0,stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu);
}
}

__device__ void eigs33_insitu_swirling_c2_kernel_0(real *rmat,real *rex,real *rimx){
//Device kernel for eigs33_insitu_swirling_c2_kernel_0
real ddbl;real pdbl;real qdbl;
real otrd;real pi;real a;
real temps;real tempw;real b;
real somma;real c;real p;
real q;real sqdel;real sqp;
real teta;
int ii;int jj;int k;
#undef __LI_REX
#define __LI_REX(i) (i-(1))
#undef __LI_RIMX
#define __LI_RIMX(i) (i-(1))
real rey[3];
#undef __LI_REY
#define __LI_REY(i) (i-(1))
real rimy[3];
#undef __LI_RIMY
#define __LI_RIMY(i) (i-(1))
real reu[3];
#undef __LI_REU
#define __LI_REU(i) (i-(1))
real rimu[3];
#undef __LI_RIMU
#define __LI_RIMU(i) (i-(1))
real rev[3];
#undef __LI_REV
#define __LI_REV(i) (i-(1))
real rimv[3];
#undef __LI_RIMV
#define __LI_RIMV(i) (i-(1))
#undef __LI_RMAT
#define __LI_RMAT(i,j) ((i-(1))+3*(j-(1)))
real s[3*3];
#undef __LI_S
#define __LI_S(i,j) ((i-(1))+3*(j-(1)))
real ww[3*3];
#undef __LI_WW
#define __LI_WW(i,j) ((i-(1))+3*(j-(1)))
real temp[3*3];
#undef __LI_TEMP
#define __LI_TEMP(i,j) ((i-(1))+3*(j-(1)))


otrd = 1./3.0;
pi = acos(-1.0);
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
s[__LI_S(ii,jj)] = 0.5*(rmat[__LI_RMAT(ii,jj)]+rmat[__LI_RMAT(jj,ii)]);
ww[__LI_WW(ii,jj)] = 0.5*(rmat[__LI_RMAT(ii,jj)]-rmat[__LI_RMAT(jj,ii)]);
}
}
a = -(s[__LI_S(1,1)]+s[__LI_S(2,2)]+s[__LI_S(3,3)]);
temps = 0.;
tempw = 0.;
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
temps = temps + s[__LI_S(ii,jj)]*s[__LI_S(ii,jj)];
tempw = tempw + ww[__LI_WW(ii,jj)]*ww[__LI_WW(ii,jj)];
}
}
b=0.5*(((a)*(a))-temps+tempw);
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
temp[__LI_TEMP(ii,jj)]=0.;
for(int k=1; k<3+1; k++){
temp[__LI_TEMP(ii,jj)] = temp[__LI_TEMP(ii,jj)]+s[__LI_S(ii,k)]*s[__LI_S(k,jj)]+3.*ww[__LI_WW(ii,k)]*ww[__LI_WW(k,jj)];
}
}
}
somma=0.;
for(int ii=1; ii<3+1; ii++){
for(int jj=1; jj<3+1; jj++){
somma= somma+temp[__LI_TEMP(ii,jj)]*s[__LI_S(jj,ii)];
}
}
c=-1./3.*(pow(a,3)-3.*a*b+somma);
pdbl=b-((a)*(a))/3.;
qdbl=c-a*b/3.+2*pow(a,3)/27.;
ddbl=(((qdbl)*(qdbl)))/4.E0+(pow(pdbl,3))/27.E0;
p = pdbl;
q = qdbl;
if(ddbl>0.E0) {
sqdel = sqrt(ddbl);
reu[__LI_REU(1)] =-0.5*q+sqdel;
rev[__LI_REV(1)] =-0.5*q-sqdel;
reu[__LI_REU(1)]=SIGN(1.0,reu[__LI_REU(1)])*pow((abs(reu[__LI_REU(1)])),otrd);
rev[__LI_REV(1)]=SIGN(1.0,rev[__LI_REV(1)])*pow((abs(rev[__LI_REV(1)])),otrd);
reu[__LI_REU(2)] = -0.5*reu[__LI_REU(1)];
rev[__LI_REV(2)] = -0.5*rev[__LI_REV(1)];
reu[__LI_REU(3)] = reu[__LI_REU(2)];
rev[__LI_REV(3)] = rev[__LI_REV(2)];
rimu[__LI_RIMU(1)] = 0.;
rimv[__LI_RIMV(1)] = 0.;
rimu[__LI_RIMU(2)] = sqrt(3.0)/2.*reu[__LI_REU(1)];
rimv[__LI_RIMV(2)] = sqrt(3.0)/2.*rev[__LI_REV(1)];
rimu[__LI_RIMU(3)] = -rimu[__LI_RIMU(2)];
rimv[__LI_RIMV(3)] = -rimv[__LI_RIMV(2)];
rey[__LI_REY(1)] = reu[__LI_REU(1)]+rev[__LI_REV(1)];
rimy[__LI_RIMY(1)] = rimu[__LI_RIMU(1)]+rimv[__LI_RIMV(1)];
rey[__LI_REY(2)] = reu[__LI_REU(2)]+rev[__LI_REV(3)];
rimy[__LI_RIMY(2)] = rimu[__LI_RIMU(2)]+rimv[__LI_RIMV(3)];
rey[__LI_REY(3)] = reu[__LI_REU(3)]+rev[__LI_REV(2)];
rimy[__LI_RIMY(3)] = rimu[__LI_RIMU(3)]+rimv[__LI_RIMV(2)];
}else {
if (q==0.) {
rey[__LI_REY(1)] = 0.;
rey[__LI_REY(2)] = sqrt(-p);
rey[__LI_REY(3)] = -rey[__LI_REY(2)];
}else {
sqp = 2.*sqrt(-p/3.);
sqdel = sqrt(-ddbl);
if (q<0.) {
teta = atan(-2.*sqdel/q);
}else {
teta = pi+atan(-2.*sqdel/q);
}
rey[__LI_REY(1)] = sqp*cos(teta/3.);
rey[__LI_REY(2)] = sqp*cos((teta+2*pi)/3.);
rey[__LI_REY(3)] = sqp*cos((teta+4*pi)/3.);
}
rimy[__LI_RIMY(1)] = 0.;
rimy[__LI_RIMY(2)] = 0.;
rimy[__LI_RIMY(3)] = 0.;
}
rex[__LI_REX(1)] = rey[__LI_REY(1)]-(a/3.);
rimx[__LI_RIMX(1)] = rimy[__LI_RIMY(1)];
rex[__LI_REX(2)] = rey[__LI_REY(2)]-(a/3.);
rimx[__LI_RIMX(2)] = rimy[__LI_RIMY(2)];
rex[__LI_REX(3)] = rey[__LI_REY(3)]-(a/3.);
rimx[__LI_RIMX(3)] = rimy[__LI_RIMY(3)];
if (rimy[__LI_RIMY(2)]<0.) {
rex[__LI_REX(2)] = rey[__LI_REY(3)]-(a/3.);
rimx[__LI_RIMX(2)] = rimy[__LI_RIMY(3)];
rex[__LI_REX(3)] = rey[__LI_REY(2)]-(a/3.);
rimx[__LI_RIMX(3)] = rimy[__LI_RIMY(2)];
}


}




__global__ void __launch_bounds__(256) insitu_swirling_c2_kernel(int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,int *vis_tag_gpu,int *wall_tag_gpu,real *dzitdz_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
//Kernel for insitu_swirling_c2_kernel
real uu;real vv;real ww;
real ux;real uy;real uz;
real vx;real vy;real vz;
real wx;real wy;real wz;
real ucsi;real ueta;real uzit;
real vcsi;real veta;real vzit;
real wcsi;real weta;real wzit;
real ccl;real div3l;real cli;
real clj;real clk;real epsi2;
real omx;real omy;real omz;
real omod2;real div2;real div;
int i;int j;int k;
int l;int lmaxi;int lmaxj;
real eigr_a[3];
#undef __LI_EIGR_A
#define __LI_EIGR_A(i) (i-(1))
real eigi_a[3];
#undef __LI_EIGI_A
#define __LI_EIGI_A(i) (i-(1))
real astar[3*3];
#undef __LI_ASTAR
#define __LI_ASTAR(i,j) ((i-(1))+3*(j-(1)))

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=1; j<ny+1; j++){
ucsi = 0.0;
vcsi = 0.0;
wcsi = 0.0;
ueta = 0.0;
veta = 0.0;
weta = 0.0;
uzit = 0.0;
vzit = 0.0;
wzit = 0.0;
lmaxi = visc_order/2;
lmaxj = visc_order/2;
if (j == 1) lmaxi = vis_tag_gpu[__I1_VIS_TAG(i)];
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) lmaxj = min(j,visc_order/2);
for(int l=1; l<visc_order/2+1; l++){
cli = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxi)];
clj = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxj)];
clk = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,visc_order/2 )];
ucsi = ucsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,2)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,2)]);
vcsi = vcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,3)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,3)]);
wcsi = wcsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,4)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,4)]);
ueta = ueta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,2)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,2)]);
veta = veta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,3)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,3)]);
weta = weta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,4)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,4)]);
uzit = uzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,2)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,2)]);
vzit = vzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,3)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,3)]);
wzit = wzit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,4)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,4)]);
}
ux = ucsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + ueta *detadxc2_gpu[__I2_DETADXC2(i,j)];
vx = vcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + veta *detadxc2_gpu[__I2_DETADXC2(i,j)];
wx = wcsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + weta *detadxc2_gpu[__I2_DETADXC2(i,j)];
uy = ucsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + ueta *detadyc2_gpu[__I2_DETADYC2(i,j)];
vy = vcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + veta *detadyc2_gpu[__I2_DETADYC2(i,j)];
wy = wcsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + weta *detadyc2_gpu[__I2_DETADYC2(i,j)];
uz = uzit *dzitdz_gpu[__I1_DZITDZ(k)];
vz = vzit *dzitdz_gpu[__I1_DZITDZ(k)];
wz = wzit *dzitdz_gpu[__I1_DZITDZ(k)];
div = ux+vy+wz;
div3l = div/3.0;
epsi2=((u0)*(u0));
omz = vx-uy;
omx = wy-vz;
omy = uz-wx;
omod2 = omx*omx+omy*omy+omz*omz;
div2 = div*div;
astar[__LI_ASTAR(1,1)] = ux-div3l;
astar[__LI_ASTAR(1,2)] = uy;
astar[__LI_ASTAR(1,3)] = uz;
astar[__LI_ASTAR(2,1)] = vx;
astar[__LI_ASTAR(2,2)] = vy-div3l;
astar[__LI_ASTAR(2,3)] = vz;
astar[__LI_ASTAR(3,1)] = wx;
astar[__LI_ASTAR(3,2)] = wy;
astar[__LI_ASTAR(3,3)] = wz-div3l;
eigs33_insitu_swirling_c2_kernel_0(astar,eigr_a,eigi_a);
psi_gpu[__I4_PSI(i,j,k,mpsi)] = 2.0*max(0.0,eigi_a[__LI_EIGI_A(2)]);
}

}
}


extern "C"{
void insitu_swirling_c2_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,int *vis_tag_gpu,int *wall_tag_gpu,real *dzitdz_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((insitu_swirling_c2_kernel),grid,block,0,stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,vis_tag_gpu,wall_tag_gpu,dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu);
}
}



__global__ void __launch_bounds__(256) insitu_schlieren_kernel(int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
//Kernel for insitu_schlieren_kernel
real rhox;real rhoy;real rhoz;
real ccl;
int i;int j;int k;
int l;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=1; j<ny+1; j++){
rhox = 0.0;
rhoy = 0.0;
rhoz = 0.0;
for(int l=1; l<visc_order/2+1; l++){
ccl = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,visc_order/2)];
rhox = rhox+ccl*(w_aux_gpu[__I4_W_AUX(i+l,j,k,1)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,1)]);
rhoy = rhoy+ccl*(w_aux_gpu[__I4_W_AUX(i,j+l,k,1)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,1)]);
rhoz = rhoz+ccl*(w_aux_gpu[__I4_W_AUX(i,j,k+l,1)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,1)]);
}
rhox = rhox*dcsidx_gpu[__I1_DCSIDX(i)];
rhoy = rhoy*detady_gpu[__I1_DETADY(j)];
rhoz = rhoz*dzitdz_gpu[__I1_DZITDZ(k)];
psi_gpu[__I4_PSI(i,j,k,mpsi)]=exp(-sqrt((((rhox))*((rhox)))+(((rhoy))*((rhoy)))+(((rhoz))*((rhoz)))));
}

}
}


extern "C"{
void insitu_schlieren_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,real *dcsidx_gpu,real *detady_gpu,real *dzitdz_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((insitu_schlieren_kernel),grid,block,0,stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,dcsidx_gpu,detady_gpu,dzitdz_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu);
}
}



__global__ void __launch_bounds__(256) insitu_schlieren_c2_kernel(int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,int *vis_tag_gpu,int *wall_tag_gpu,real *dzitdz_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
//Kernel for insitu_schlieren_c2_kernel
real rhocsi;real rhoeta;real rhozit;
real rhox;real rhoy;real rhoz;
real cli;real clj;real clk;
int i;int j;int k;
int l;int lmaxi;int lmaxj;

i = __GIDX(x,1);
k = __GIDX(y,1);


if(loop_cond(i,nx,1)&&loop_cond(k,nz,1)){
for(int j=1; j<ny+1; j++){
rhocsi = 0.0;
rhoeta = 0.0;
rhozit = 0.0;
lmaxi = visc_order/2;
lmaxj = visc_order/2;
if (j == 1) lmaxi = vis_tag_gpu[__I1_VIS_TAG(i)];
if (wall_tag_gpu[__I1_WALL_TAG(i)] < 1) lmaxj = min(j,visc_order/2);
for(int l=1; l<visc_order/2+1; l++){
cli = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxi)];
clj = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,lmaxj)];
clk = coeff_deriv1_gpu[__I2_COEFF_DERIV1(l,visc_order/2 )];
rhocsi = rhocsi +cli*(w_aux_gpu[__I4_W_AUX(i+l,j,k,1)]-w_aux_gpu[__I4_W_AUX(i-l,j,k,1)]);
rhoeta = rhoeta +clj*(w_aux_gpu[__I4_W_AUX(i,j+l,k,1)]-w_aux_gpu[__I4_W_AUX(i,j-l,k,1)]);
rhozit = rhozit +clk*(w_aux_gpu[__I4_W_AUX(i,j,k+l,1)]-w_aux_gpu[__I4_W_AUX(i,j,k-l,1)]);
}
rhox = rhocsi *dcsidxc2_gpu[__I2_DCSIDXC2(i,j)] + rhoeta *detadxc2_gpu[__I2_DETADXC2(i,j)];
rhoy = rhocsi *dcsidyc2_gpu[__I2_DCSIDYC2(i,j)] + rhoeta *detadyc2_gpu[__I2_DETADYC2(i,j)];
rhoz = rhozit *dzitdz_gpu[__I1_DZITDZ(k)];
psi_gpu[__I4_PSI(i,j,k,mpsi)]=exp(-sqrt((((rhox))*((rhox)))+(((rhoy))*((rhoy)))+(((rhoz))*((rhoz)))));
}

}
}


extern "C"{
void insitu_schlieren_c2_kernel_wrapper(hipStream_t stream,int nv,int nx,int ny,int nz,int visc_order,int ng,int npsi,int mpsi,real u0,int *vis_tag_gpu,int *wall_tag_gpu,real *dzitdz_gpu,real *dcsidxc2_gpu,real *detadyc2_gpu,real *detadxc2_gpu,real *dcsidyc2_gpu,real *psi_gpu,real *w_aux_gpu,real *coeff_deriv1_gpu,real *x_gpu){
dim3 block(EULERWENO_THREADS_X,EULERWENO_THREADS_Y);
dim3 grid(divideAndRoundUp((nx)-(1)+1,block.x),divideAndRoundUp((nz)-(1)+1,block.y));

hipLaunchKernelGGL((insitu_schlieren_c2_kernel),grid,block,0,stream,nv,nx,ny,nz,visc_order,ng,npsi,mpsi,u0,vis_tag_gpu,wall_tag_gpu,dzitdz_gpu,dcsidxc2_gpu,detadyc2_gpu,detadxc2_gpu,dcsidyc2_gpu,psi_gpu,w_aux_gpu,coeff_deriv1_gpu,x_gpu);
}
}


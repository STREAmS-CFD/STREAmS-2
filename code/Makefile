#----------------------------------------------------------------------
# Edit makefile.inc starting from avilable templates to configure Make process
#----------------------------------------------------------------------
include makefile.inc
#----------------------------------------------------------------------

#----------------------------------------------------------------------
ifeq (${PREC},"single") 
    PRECFLAGS = -DSINGLE_PRECISION 
endif

#----------------------------------------------------------------------
# Default C settings
#----------------------------------------------------------------------
CC       = gcc
CFLAGS   = -O2
CXX      = g++
CXXFLAGS = -std=c++17 -frounding-math -O2
# original INSITU was using c++11
#CXXFLAGS = -std=c++11 -O2
#----------------------------------------------------------------------
# Trick to create a space
E := 
S := $E $E
#----------------------------------------------------------------------

ifeq (${COMPILE},"cray-cuda") # mpif90 shipped with PGI compiler
    # Pitz-Daint instructions (2022)
    #module swap PrgEnv-cray PrgEnv-nvidia
    #module load cudatoolkit/10.2.89_3.28-2.1__g52c0314
    #export CPATH=`echo $CPATH | sed "s@$CRAY_NVIDIA_PREFIX/compilers/include@@"`
    #export CPATH=$CRAY_NVIDIA_PREFIX/cuda/10.2/targets/x86_64-linux/include:$CPATH
    #export MPICH_RDMA_ENABLED_CUDA=1 #export this variable in the submission script
    MODULE_FLAG = -module
    CUDA_HOME = $(CUDATOOLKIT_HOME) 
    INCLUDE   = -I $(MPICH_DIR)/include 
    FCLIBS    = -L$(MPICH_DIR)/lib -lmpich 
    FC        = ftn$S
    ifeq (${MODE},"opt") 
        OLEVEL = CUDA_HOME=$(CUDATOOLKIT_HOME) -O3 -Mcuda=ccall 
    else
        OLEVEL = -Mcuda=ccall -g 
    endif
endif

ifeq (${COMPILE},"nv-cuda") # mpif90 shipped with PGI compiler
    FC = mpif90 
    MODULE_FLAG = -module$S
    ifeq (${MODE},"opt") 
        #OLEVEL = -O3 -Mcuda=cc70,cuda10.1,ptxinfo -DUSE_CUDA -DCUDA_ASYNC #-DNOCUDAAWAREMPI
        OLEVEL = -O3 -Mcuda=ccall,cuda$(CUDA_VER)
    else
        OLEVEL = -Mcuda=ccall,cuda$(CUDA_VER) -g -Minfo
        #OLEVEL = -Mcuda=ccall,cuda11.0,ptxinfo -g -Minfo
    endif
endif

ifeq (${COMPILE},"nv") # mpif90 shipped with PGI compiler or other MPI
    FC = mpif90 
    MODULE_FLAG = -module$S
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 
    else
        OLEVEL = -O2 -Mbounds -Ktrap=fp -g
    endif
endif

ifeq (${COMPILE},"intel") # mpiifort from IntelMPI
    FC = mpiifort
    MODULE_FLAG = -module$S
    ifeq (${MODE},"opt") 
        OLEVEL = -O3
    else
        OLEVEL = -O2 -check bounds,uninit -g -fpe0 -traceback
    endif
endif

ifeq (${COMPILE},"gnu") # mpif90 from MPI library
    FC = mpif90 
    MODULE_FLAG = -J
    ifeq (${MODE},"opt") 
        OLEVEL = -O3 #-ffree-line-length-256
    else
        OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core -ffpe-trap=zero,invalid,overflow -ffree-line-length-256
        #OLEVEL = -O2 -fbounds-check -g -fbacktrace -fdump-core
    endif
endif

ifeq (${COMPILE},"ibmxl") # mpif90 from MPI library
    FC = mpif90
    MODULE_FLAG = -qmoddir=
    ifeq (${MODE},"opt") 
        OLEVEL = -O4
    else
        OLEVEL = -O2 -g
    endif
endif

ifeq (${COMPILE},"hip")
    FC  = $(HIPFC)
    CC  = $(HIPCC)
    CXX = $(HIPCC)
    MODULE_FLAG := -J
    ifeq (${MODE},"opt")
        OLEVEL = -std=f2008 -O3 -fallow-argument-mismatch
    else
	      CXXFLAGS := $(CXXFLAGS) -save-temps
        OLEVEL = -std=f2008 -O2 -fallow-argument-mismatch -fbounds-check -g -fbacktrace -fdump-core -ffpe-trap=zero,invalid,overflow -save-temps
    endif
endif

ifdef MPI_CMD
    FC = ${MPI_CMD}
endif

RM = rm -f

EQUATION := $(patsubst "%",%,$(EQUATION))

#----------------------------------------------------------------------
# Base code
#----------------------------------------------------------------------
TP_MODS  = third-party/cfgio/string_conv_mod.o third-party/cfgio/cfgio_mod.o third-party/crandom/crandom.o third-party/crandom/crandom_f.o third-party/cpumem/getmeminterface.o third-party/cpumem/getmemory.o
ifeq (${CGAL},"yes")
    TP_MODS := $(TP_MODS) third-party/cgal/cgal_c_wrappers.o third-party/cgal/cgal_f90_wrappers.o
    INCLUDE := $(INCLUDE) $(CGAL_INCLUDE)
    LDFLAGS := $(LDFLAGS) $(CGAL_LDFLAGS)
else
    TP_MODS := $(TP_MODS) third-party/cgal/cgal_c_wrappers_dummy.o third-party/cgal/cgal_f90_wrappers.o
endif
ifeq (${CATALYST},"yes")
    TP_MODS := $(TP_MODS) third-party/catalyst/catalyst_adaptor_c.o third-party/catalyst/catalyst_adaptor.o
    INCLUDE := $(INCLUDE) $(CATALYST_INCLUDE)
    LDFLAGS := $(LDFLAGS) $(CATALYST_LDFLAGS)
else
    TP_MODS := $(TP_MODS) third-party/catalyst/catalyst_adaptor_c_dummy.o third-party/catalyst/catalyst_adaptor_dummy.o
endif
LIB_MODS = src/parameters.o src/grid.o src/field.o
EQU_MODS = src/$(EQUATION)/$(EQUATION).o 

ifeq (${BACKEND},"cpu") 
    MAIN     = src/$(EQUATION)/main_cpu.o
    LIB_MODS := $(LIB_MODS) src/base_cpu.o
    EQU_MODS := $(EQU_MODS) src/$(EQUATION)/kernels_cpu.o src/$(EQUATION)/$(EQUATION)_cpu.o
else ifeq (${BACKEND},"gpu") 
    MAIN     = src/$(EQUATION)/main_gpu.o
    LIB_MODS := $(LIB_MODS) src/base_gpu.o
    EQU_MODS := $(EQU_MODS) src/$(EQUATION)/kernels_gpu.o src/$(EQUATION)/$(EQUATION)_gpu.o
else ifeq (${BACKEND},"amd") 
    MAIN     = src/$(EQUATION)/main_amd.o
    LIB_MODS := $(LIB_MODS) src/base_amd.o src/base_amd_cpp.o
    EQU_MODS := $(EQU_MODS) src/$(EQUATION)/kernels_amd.o src/$(EQUATION)/kernels_amd_cpp.o src/$(EQUATION)/$(EQUATION)_amd.o 
endif

OBJS = $(TP_MODS) $(LIB_MODS) $(EQU_MODS) $(MAIN)
OBJS_BUILD = $(OBJS:%.o=build/%.o)
EXEC = streams_2.exe

$(EXEC): $(OBJS_BUILD) Makefile makefile.inc
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -o $@ $(OBJS_BUILD) $(FCLIBS) $(LDFLAGS)

build/%.o:%.F90 Makefile makefile.inc
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -c $< -o $@ $(MODULE_FLAG)build/mods

build/%.o:%.c Makefile makefile.inc
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

build/%.o:%.cpp Makefile makefile.inc
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

build/%.o:%.cxx Makefile makefile.inc
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# Module dependencies
build/src/grid.o: build/mods/streams_parameters.mod
build/src/field.o: build/mods/streams_parameters.mod build/mods/streams_grid_object.mod
build/src/base_cpu.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod
build/src/base_gpu.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod
build/src/base_amd.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod
build/src/$(EQUATION)/$(EQUATION).o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod build/mods/streams_grid_object.mod
build/src/$(EQUATION)/$(EQUATION)_cpu.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod build/mods/streams_grid_object.mod build/mods/streams_base_cpu_object.mod build/mods/streams_equation_$(EQUATION)_object.mod build/mods/streams_kernels_cpu.mod
build/src/$(EQUATION)/$(EQUATION)_gpu.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod build/mods/streams_grid_object.mod build/mods/streams_base_gpu_object.mod build/mods/streams_equation_$(EQUATION)_object.mod build/mods/streams_kernels_gpu.mod
build/src/$(EQUATION)/$(EQUATION)_amd.o: build/mods/streams_parameters.mod build/mods/streams_field_object.mod build/mods/streams_grid_object.mod build/mods/streams_base_amd_object.mod build/mods/streams_equation_$(EQUATION)_object.mod build/mods/streams_kernels_amd.mod
build/src/$(EQUATION)/kernels_cpu.o: build/mods/streams_parameters.mod
build/src/$(EQUATION)/kernels_gpu.o: build/mods/streams_parameters.mod
build/src/$(EQUATION)/kernels_amd.o: build/mods/streams_parameters.mod
build/src/$(EQUATION)/main_cpu.o: build/mods/streams_equation_$(EQUATION)_cpu_object.mod
build/src/$(EQUATION)/main_gpu.o: build/mods/streams_equation_$(EQUATION)_gpu_object.mod
build/src/$(EQUATION)/main_amd.o: build/mods/streams_equation_$(EQUATION)_amd_object.mod

.PHONY: clean
clean:
	$(RM) $(EXEC) $(MAIN) $(OBJS_BUILD) build/mods/*.mod

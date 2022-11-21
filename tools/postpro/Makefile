#----------------------------------------------------------------------
# SUPPORTED COMPILE: gnu
# SUPPORTED PRECISIONS (PREC) : single, double [default=double]
#----------------------------------------------------------------------
COMPILE = "gnu"
PREC    = "double"
#----------------------------------------------------------------------

ifeq (${PREC},"single") 
    PRECFLAGS = -DSINGLE_PRECISION 
endif

ifeq (${COMPILE},"gnu")
    FC = gfortran
    OLEVEL = -O3
    MODULE_FLAG = -J
endif

RM = rm -f

#----------------------------------------------------------------------
# Base code
#----------------------------------------------------------------------
MAIN = src/main.o
OBJ_FILES = src/parameters.o src/global_variables.o src/uty.o src/reader.o src/derivatives.o src/comp_transform.o src/postpro_cha.o src/postpro_bl.o 

OBJ_MODS = third-party/cfgio/string_conv_mod.o third-party/cfgio/cfgio_mod.o
OBJS = $(OBJ_MODS) $(OBJ_FILES) $(MAIN)
OBJS_BUILD = $(OBJS:%.o=build/%.o)
EXEC = postpro.exe
#----------------------------------------------------------------------

$(EXEC): $(OBJS_BUILD) Makefile 
	$(FC) $(PRECFLAGS) $(OLEVEL) -o $@ $(OBJS_BUILD) $(FCLIBS) $(LDFLAGS)

build/%.o:%.F90 Makefile
	$(FC) $(PRECFLAGS) $(INCLUDE) $(OLEVEL) -c $< -o $@ $(MODULE_FLAG) build/mods

.PHONY: clean
clean:
	$(RM) $(EXEC) $(MAIN) $(OBJS_BUILD) build/mods/*.mod

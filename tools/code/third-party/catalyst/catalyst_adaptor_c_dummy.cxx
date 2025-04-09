#include <cstdlib>
#include <cstdio>

// Adaptor for getting Fortran simulation code into ParaView CoProcessor.
#ifdef SINGLE_PRECISION
#define REAL float
#else
#define REAL double
#endif



// These will be called from the Fortran "glue" code"
// Completely dependent on data layout, structured vs. unstructured, etc.
// since VTK/ParaView uses different internal layouts for each.

// Creates the data container for the CoProcessor.
extern "C" void createcpimagedata(int* nxstart, int* nxend, int* nzstart, int* nzend,
 int* nx, int* nz,
 REAL* dx, REAL* dz,
 REAL* origin, char* zone, int* nrank)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtoimage(char* zone, REAL* scalars, char* name, int* nrank)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// Creates the data container for the CoProcessor.
extern "C" void createcprectilineardata(int* nxstart, int* nxend, int* nystart, int* nyend, int* nzstart, int* nzend,
 int* nxstartg, int* nystartg, int* nzstartg,
 int* nxendg, int* nyendg, int* nzendg,
 REAL* x, REAL* y, REAL* z,
 int* nrank, int* nproc,
 char* zone)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtorectilinear(char* zone, REAL* scalars, char* name, int* nrank)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// Creates the data container for the CoProcessor.
extern "C" void createcpstructureddata(int* nxp2, int* nyp2, int* nzp2,
 int* nxstart, int* nystart, int* nzstart, int* nxend, int* nyend, int* nzend,
 int* nxmaxp2, int* nymaxp2, int* nzmaxp2,
 REAL* x, REAL* y, REAL* z,
 int* nrank, int* nproc,
 char* zone,
 REAL* xyzc)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// Add field(s) to the data container.
// Separate from above because this will be dynamic, grid is static.
// By hand name mangling for fortran.
extern "C" void addfieldtostructured(char* zone, REAL* scalars, char* name, int* nrank)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// wrapper to empty paraview stuff
extern "C" void coprocessorinitializewithpython(char* vtkpipeline, int* length)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// wrapper to empty paraview stuff
extern "C" void coprocessorfinalize()
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// wrapper to empty paraview stuff
extern "C" void requestdatadescription(int* icyc, REAL* time, int* flag)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// wrapper to empty paraview stuff
extern "C" void needtocreategrid(int* flag)
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}

// wrapper to empty paraview stuff
extern "C" void coprocess()
{
        printf("CATALYST compilation not enabled! Exiting...\n");
        exit(0);
}


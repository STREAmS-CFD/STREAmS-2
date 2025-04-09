#include <stdlib.h>                                                             
#include <stdio.h>                                                              

#ifdef SINGLE_PRECISION
#define REAL float
#else
#define REAL double
#endif

REAL rescale;
                                                                                
void init_crandom(int seed) {
    //printf("seed %d :",seed);
    if(seed == 0) {
        rescale = 0.;
    } else { 
        rescale = 1.;
    }
    srand(seed);
};

void get_crandom(REAL *a) 
{                                                                               
    int i;
    i = rand();
    // printf("rand int %d :",i);
    *a = rescale * (REAL)i/RAND_MAX;
    // printf("rescale %lf :",rescale);
    // printf("rand float %lf :",*a);
};       

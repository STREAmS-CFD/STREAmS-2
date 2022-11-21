#include <stdlib.h>                                                             
#include <stdio.h>                                                              
                                                                                
void init_crandom(int seed) {
    //printf("seed %d :",seed);
    srand(seed);
};

void get_crandom(double *a) 
{                                                                               
    int i;
    i = rand();
    //printf("rand int %d :",i);
    *a = (double)i/RAND_MAX;
    //printf("rand float %lf :",*a);
};       

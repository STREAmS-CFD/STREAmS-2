gcc crandom.c -c
gfortran crandom_f.F90 -c
gfortran test_f.F90 -c
gfortran test_f.o crandom_f.o crandom.o 

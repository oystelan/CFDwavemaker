#!/bin/bash
# Creates DLL shared library from input file
module rm group_users
export LIBRARY_PATH=/usr/lib64

gcc -fPIC -fopenmp -std=c99 -O2 -c  Stokes5.c
g++ -fPIC -fopenmp -std=c++11 -O2 -c  CFDwavemaker.cpp

#g++ -L/usr/lib64 -shared -pthread -o comflow_wavemaker.so wavemaker.o
g++ -L/usr/lib64 -shared -fopenmp -o CFDwavemaker_omp.so CFDwavemaker.o Stokes5.o
module add group_users

#g++ -fPIC wavemaker.cpp -shared -std=c++11 -o comflow_wavemaker.so -Wl,--whole-archive -Wl,--no-whole-archive


cp CFDwavemaker_omp.so /cm/shared/apps/comflow/4.2.0-lib/.
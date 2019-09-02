#!/bin/bash
# Creates DLL shared library from input file

export LIBRARY_PATH=/usr/lib64

#gcc -pthread -std=c99 -O2 -c  Stokes5.c
#g++ -pthread -std=c++11 -O2 -c  CFDwavemaker.cpp

gcc -std=c99 -fopenmp -O2 -c  Stokes5.c
gcc -std=c++11 -fopenmp -lstdc++ -O2 -c  CFDwavemaker.cpp

#g++ -L/usr/lib64 -shared -pthread -o comflow_wavemaker.so wavemaker.o
#g++ -L/usr/lib64 -shared -pthread -o CFDwavemaker.so CFDwavemaker.o Stokes5.o
#g++ -L/usr/lib64  -pthread -o CFDwavemaker.a CFDwavemaker.o Stokes5.o
ar rvs libcfdwavemaker.a CFDwavemaker.o Stokes5.o

#g++ -fPIC wavemaker.cpp -shared -std=c++11 -o comflow_wavemaker.so -Wl,--whole-archive -Wl,--no-whole-archive

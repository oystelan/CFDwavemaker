# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
import os

# give location of dll
#mydll = cdll.comflow
temp = os.path.abspath(__file__)
temp = os.path.realpath(temp)
temp = os.path.dirname(temp)

# we load the shared library 
mydll = cdll.LoadLibrary("../../../builds/linux64/libCFDwavemaker_openmp.so")     

# and define the functions which CFDwavemaker.h tells us are in the library

def wave_VeloX(mydll,t,x,y,z):
    aa = mydll.wave_VeloX
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_VeloY(mydll,t,x,y,z):
    aa = mydll.wave_VeloY
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_VeloZ(mydll,t,x,y,z):
    aa = mydll.wave_VeloZ
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t))

def wave_SurfElev(mydll,t,x,y):
    aa = mydll.wave_SurfElev
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(t))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0))

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa()


# All done. lets test:

# we start by initializing (reads data from waveinput.dat and computes kinematics at t=0 since lsgrid is used) 
init_dll(mydll)


# Now that the initialization is done, we could extract interpolated results at any position within the bounds of the lsgrid. Linear interpolation by default.
print('such as wave elevation: ' + wave_SurfElev(mydll,0,1.12,3.223))
print('or kinematics: ' + wave_VeloX(mydll,0,1.1345,1.234,-0.45))
# you see what i mean...., but more importantly the [vtk output] is used in the waveinput.dat file, which means that kinematics for the entire lsgrid is dumped into a vts file every time the timestep of the lsgrid is updated. (hint: look in the ./vtk folder after running python3 case3.py)


# all done. remember to clean up after us.
clean_up(mydll)






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
#path = os.path.join(temp, "libCFDwavemaker_openmp.so")
print(temp)
#lib = CDLL("./comflow_wavemaker.so")

mydll = cdll.LoadLibrary("../../../builds/linux64/libCFDwavemaker_swd_openmp.so")     

# and define the functions which CFDwavemaker.h tells us are in the library
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

# we start by initializing (reads data from waveinput.dat) 
init_dll(mydll)

x = 4.8
y = 0.23

time = np.arange(0,20,0.2)

# this time we only extract elevation at every 0.2 sec, but have turned on [lsgrid] and [vtk output] (see. waveinput.dat) to write complete kinematics to vtk files which may be viewed in paraview
elev = []
for t in time:
    elev.append(wave_SurfElev(mydll,t,x,y))
    
# now, lets view the resulting surface elevation
plt.plot(time, elev, label="surface elevation")
plt.legend()
plt.xlabel('Time [sec]')
plt.ylabel('Surface elevation [m]')
plt.grid(True)
plt.savefig("./result_eta.png")

# all done. remember to clean up after us.
clean_up(mydll)






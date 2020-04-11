# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np

# give location of dll
mydll = cdll.CFDwavemaker
def tic():
    import time
    #Homemade version of matlab tic and toc functions    
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time    
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")

def velocityx(mydll,t,x,y,z):
    aa = mydll.VelocityX
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))    


    
def velocityy(mydll,t,x,y,z):
    aa = mydll.VelocityY
    aa.restype = c_double
    aa.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def velocityz(mydll,t,x,y,z):
    aa = mydll.VelocityZ
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def waveelev(mydll,t,x,y):
    aa = mydll.SurfaceElevation
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))

def volfrac(mydll,x,y,z,t,delta):
    aa = mydll.VolumeFraction
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(delta))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0));

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa();


# We start of by calling the init function. This will read the waveinput.dat file and perform necessary initialization
print(init_dll(mydll))

# We create a time vector from 0 to 100 sec, with 0.5sec spacing
time = np.arange(0, 100, 0.5)

# lets extract the surface elevation at x=y=0
wave_elev = np.zeros(len(time))
for wave,t in zip(wave_elev,time):
    wave = waveelev(mydll,t,0.0,0.0) # x and y position set to 0.

# or perhaps the horizontal particle velocity at position x = 0, y=0, z = -10.

ux = np.zeros(len(time))
for u,t in zip(ux,time):
    u = waveelev(mydll,t,0.0,0.0,-11.0) # x and y position set to 0.

# You get the idea.

# Plotting some kinematics data
plt.plot(time,wave)
plt.plot(time,ux)
plt.show()

# now that we are done, we close down the link to the library before exiting
clean_up(mydll)

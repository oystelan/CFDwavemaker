# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 14:58:16 2016

@author: oland
"""

from ctypes import*
import matplotlib.pylab as plt
import numpy as np
import shutil
import os.path

# give location of dll
#mydll = cdll.CFDwavemaker
dll_name = "CFDwavemaker.dll"
dllpath = os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../builds/win64",dll_name)
print(dllpath)
mydll = CDLL(dllpath)


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

def velocityX(mydll,t,x,y,z):
    aa = mydll.VelocityX
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))
    
def velocityY(mydll,t,x,y,z):
    aa = mydll.VelocityY
    aa.restype = c_double
    aa.argtypes = [c_int, c_int, c_int, c_double, c_double, c_double, c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def velocityZ(mydll,t,x,y,z):
    aa = mydll.VelocityZ
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_int,c_double,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_int(0),c_double(x),c_double(y),c_double(z),c_double(t))

def waveelev(mydll,t,x,y):
    aa = mydll.SurfaceElevation
    aa.restype = c_double
    aa.argtypes = [c_int,c_int,c_double,c_double,c_double]
    return aa(c_int(0),c_int(0),c_double(x),c_double(y),c_double(t))

def update_probes(mydll,t):
    aa = mydll.update_probes
    aa.restype = c_double
    aa.argtypes = [c_double]
    return aa(c_double(t))

def volfrac(mydll,x,y,z,t,delta):
    aa = mydll.VolumeFraction
    aa.restype = c_double
    aa.argtypes = [c_double,c_double,c_double,c_double,c_double]
    return aa(c_double(x),c_double(y),c_double(z),c_double(t),c_double(delta))

def init_dll(mydll):
    aa = mydll.Init
    aa.restype = c_int
    aa.argtypes = [POINTER(c_double),POINTER(c_double)]
    return aa(c_double(0),c_double(0))

def clean_up(mydll):
    aa = mydll.Cleanup
    aa.restype = c_int
    return aa()


# We start of by calling the init function. This will read the waveinput.dat file and perform necessary initialization


# We create a time vector from 0 to 100 sec, with 0.5sec spacing
# time = np.arange(0, 100, 0.5)

# # lets extract the surface elevation at x=y=0
# wave_elev = np.zeros(len(time))
# for wave,t in zip(wave_elev,time):
#     wave = waveelev(mydll,t,0.0,0.0) # x and y position set to 0.

# or perhaps the horizontal particle velocity at position x = 0, y=0, z = -10.




# t = np.arange(0, 300,0.87)
# x = 0.
# #shutil.copy2('./waveinput1.dat','./waveinput.dat')
# print(init_dll(mydll))

# eta = []
# for tt in t:
#     dd = waveelev(mydll,tt,x,0.0)
#     eta.append(dd) # x and y position set to 0.
# plt.plot(t, eta)
# plt.show()
# clean_up(mydll)

#exit()

#
time = np.arange(0, 10., 0.06)
x = 0.
#shutil.copy2('./waveinput3.dat','./waveinput.dat')
tic()
print(init_dll(mydll))
toc()

for t in time:
    #update_probes(mydll, t)
    wave = waveelev(mydll,t,x,0.0)
print("Alle done")
# z = np.arange(-76.4,wave,2.)
# # wave = velocityX(mydll,t,x,0.0,-76.3)
# # print(wave)
# exit()
# ux = np.zeros(len(z))
# ux = []
# for zz in z:
#     dd = velocityZ(mydll,t,x,0.0,zz)
#     print(zz, dd)
#     ux.append(dd) # x and y position set to 0.
# plt.plot(ux, z)
clean_up(mydll)

#
#
# shutil.copy2('./waveinput2.dat','./waveinput.dat')
# print(init_dll(mydll))
# wave = waveelev(mydll,t,x,0.0)
# print(wave)
# tic()
# z = np.arange(-76.4,wave,0.1)
# print(len(z))
# ux = np.zeros(len(z))
# ux = []
# for zz in z:
#     ux.append(velocityZ(mydll,t,x,0.0,zz)) # x and y position set to 0.
#
#
# toc()
# plt.plot(ux, z,'r')
# plt.show()
#
# clean_up(mydll)
